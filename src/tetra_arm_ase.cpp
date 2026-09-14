// Molecule-aware, donor-directed chromosome-arm ASE evidence from CellBouncer
// interindividual pileup sidecars. No BAM rescan and no species-panel proxy.

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <htslib/vcf.h>
#include <htswrapper/bc.h>

namespace fs = std::filesystem;

namespace {

constexpr const char* PROGRAM = "tetra_arm_ase";
constexpr const char* VERSION = "2.4.0";
constexpr const char* SCHEMA = "tetra_arm_ase_evidence_v2";
constexpr const char* QC_SCHEMA = "tetra_arm_ase_qc_v2";
constexpr const char* CELL_MANIFEST_SCHEMA = "tetra_arm_cell_manifest_v1";
constexpr const char* AMBIENT_SOURCES_SCHEMA = "tetra_arm_ambient_sources_v1";
constexpr size_t LINE_BUFFER = 1U << 20;
constexpr size_t BARCODE_BASES = 16;
#ifdef BC_LENX2
static_assert(BC_LENX2 == 2 * BARCODE_BASES,
              "tetra_arm_ase requires the 16-base CellBouncer barcode encoding");
#endif

class UserError : public std::runtime_error {
public:
    explicit UserError(const std::string& message) : std::runtime_error(message) {}
};

struct Options {
    std::string samples;
    std::string panel_bcf;
    std::string pileup_sites;
    std::string pileup_molecules;
    std::string pileup_observations;
    std::string cells;
    std::string ambient_sources;
    std::string arms;
    std::string output;
    std::string qc;
    std::string library;
    int threads = 8;
    // The producer supplies allele-support weights, not a calibrated
    // likelihood. --hard-posterior is accepted only as a historical alias.
    double hard_allele_threshold = 0.80;
    bool include_sex_chromosomes = false;
    bool self_test = false;
    bool version = false;
};

struct Interval {
    int64_t start = 0;
    int64_t end = 0;
    int arm_index = -1;
};

struct Arm {
    std::string chrom;
    std::string name;
    int64_t start = 0;
    int64_t end = 0;
};

struct Site {
    int tid = -1;
    int pos = -1;
    int arm = -1;
    char ref = 'N';
    char alt = 'N';
    std::vector<int8_t> genotype;
};

struct AmbientSource {
    int sample = -1;
    double mass = 0.0;
};

struct Cell {
    std::string barcode;
    std::string donor_a;
    std::string donor_b;
    std::string donor_pair;
    int donor_a_index = -1;
    int donor_b_index = -1;
    double ambient_c = 0.0;
    double ambient_c_se = std::numeric_limits<double>::quiet_NaN();
    bool model_eligible = false;
    std::vector<AmbientSource> ambient;
    double ambient_unmapped_mass = 0.0;
};

enum Orientation : uint8_t { A_IS_REF = 0, A_IS_ALT = 1 };

struct Observation {
    uint64_t molecule = 0;
    uint64_t site_key = 0;
    int arm = -1;
    double a = 0.0;
    double b = 0.0;
    double ambient_a = 0.5;
    double ambient_genotyped_mass = 0.0;
    Orientation orientation = A_IS_REF;
    bool qname_fallback = false;
};

struct ArmEvidence {
    uint64_t molecules = 0;
    uint64_t informative = 0;
    uint64_t molecules_ref = 0;
    uint64_t molecules_alt = 0;
    uint64_t molecules_mixed = 0;
    uint64_t a_ref = 0;
    uint64_t b_ref = 0;
    uint64_t a_alt = 0;
    uint64_t b_alt = 0;
    uint64_t a_mixed = 0;
    uint64_t b_mixed = 0;
    uint64_t ambiguous = 0;
    uint64_t qname = 0;
    double soft_a = 0.0;
    double soft_b = 0.0;
    double soft_a_ref = 0.0;
    double soft_b_ref = 0.0;
    double soft_a_alt = 0.0;
    double soft_b_alt = 0.0;
    double soft_a_mixed = 0.0;
    double soft_b_mixed = 0.0;
    double soft_a_sumsq_ref = 0.0;
    double soft_a_sumsq_alt = 0.0;
    double soft_a_sumsq_mixed = 0.0;
    double effective_a_ref = 0.0;
    double effective_b_ref = 0.0;
    double effective_a_alt = 0.0;
    double effective_b_alt = 0.0;
    double effective_a_mixed = 0.0;
    double effective_b_mixed = 0.0;
    double effective_weight_ref = 0.0;
    double effective_weight_alt = 0.0;
    double effective_weight_mixed = 0.0;
    double ambient_a_ref_sum = 0.0;
    double ambient_a_alt_sum = 0.0;
    double ambient_a_mixed_sum = 0.0;
    double ambient_genotyped_mass_sum = 0.0;
    double ambient_ref_weight = 0.0;
    double ambient_alt_weight = 0.0;
    double ambient_mixed_weight = 0.0;
    double ambient_genotyped_weight = 0.0;
    uint64_t molecule_site_units = 0;
    std::unordered_set<uint64_t> sites;
};

struct CellOutput {
    unsigned long encoded = 0;
    std::vector<ArmEvidence> arms;
};

std::string trim(const std::string& value) {
    const size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    size_t begin = 0;
    while (true) {
        const size_t tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(begin));
            return fields;
        }
        fields.push_back(line.substr(begin, tab - begin));
        begin = tab + 1;
    }
}

bool parse_double(const std::string& text, double& value) {
    errno = 0;
    char* end = nullptr;
    value = std::strtod(text.c_str(), &end);
    return errno == 0 && end != text.c_str() && *end == '\0' &&
           std::isfinite(value);
}

bool missing_value(const std::string& value) {
    std::string text = trim(value);
    std::transform(text.begin(), text.end(), text.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::toupper(character));
        });
    return text.empty() || text == "NA" || text == "N/A" || text == "NONE" ||
           text == "NULL" || text == ".";
}

long long parse_ll(const std::string& text, const std::string& context) {
    errno = 0;
    char* end = nullptr;
    const long long result = std::strtoll(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0') {
        throw UserError(context + ": expected integer, saw '" + text + "'");
    }
    return result;
}

uint64_t parse_u64(const std::string& text, const std::string& context) {
    if (text.empty() || !std::all_of(text.begin(), text.end(),
            [](unsigned char character) { return std::isdigit(character) != 0; })) {
        throw UserError(context + ": expected unsigned integer, saw '" + text + "'");
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long result = std::strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0') {
        throw UserError(context + ": expected unsigned integer, saw '" + text + "'");
    }
    return static_cast<uint64_t>(result);
}

int parse_nonnegative_int(const std::string& text, const std::string& context) {
    const long long value = parse_ll(text, context);
    if (value < 0 || value > std::numeric_limits<int>::max()) {
        throw UserError(context + ": integer is outside [0, INT_MAX]: '" + text + "'");
    }
    return static_cast<int>(value);
}

bool parse_bool(const std::string& value, const std::string& context) {
    std::string text = trim(value);
    std::transform(text.begin(), text.end(), text.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::toupper(character));
        });
    if (text == "1" || text == "TRUE" || text == "T" || text == "YES" ||
            text == "Y" || text == "PASS") return true;
    if (text == "0" || text == "FALSE" || text == "F" || text == "NO" ||
            text == "N" || text == "FAIL") return false;
    throw UserError(context + ": expected boolean, saw '" + value + "'");
}

unsigned long encode_barcode(const std::string& value,
                             const std::string& context) {
    std::string barcode = trim(value);
    if (barcode.size() != BARCODE_BASES ||
            !std::all_of(barcode.begin(), barcode.end(), [](unsigned char character) {
                return character == 'A' || character == 'C' ||
                       character == 'G' || character == 'T';
            })) {
        throw UserError(context + ": expected exactly 16 uppercase A/C/G/T bases, saw '" +
                        value + "'");
    }
    return bc_ul(barcode);
}

unsigned long parse_encoded_barcode(const std::string& value,
                                    const std::string& context) {
    const uint64_t encoded = parse_u64(value, context);
    if (encoded > static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
        throw UserError(context + ": encoded 16-base barcode exceeds UINT32_MAX");
    }
    return static_cast<unsigned long>(encoded);
}

bool gz_getline_checked(gzFile input, std::vector<char>& buffer,
                        std::string& line, const std::string& path) {
    line.clear();
    while (true) {
        char* result = gzgets(input, buffer.data(), static_cast<int>(buffer.size()));
        if (result == nullptr) {
            if (gzeof(input)) return !line.empty();
            int zlib_error = Z_OK;
            const char* detail = gzerror(input, &zlib_error);
            throw UserError("gzip read failure for " + path + ": " +
                            (detail == nullptr ? "unknown zlib error" : detail));
        }
        line.append(result);
        if (!line.empty() && line.back() == '\n') return true;
        if (gzeof(input)) return true;
    }
}

uint64_t site_key(int tid, int pos) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(tid)) << 32) |
           static_cast<uint32_t>(pos);
}

std::string make_temp_path(const std::string& destination) {
    // Reserve a collision-proof sibling path. PIDs and hostnames alone do not
    // provide O_EXCL semantics across requeues, containers, or stale files.
    std::string pattern = destination + ".tmp.XXXXXX";
    std::vector<char> mutable_pattern(pattern.begin(), pattern.end());
    mutable_pattern.push_back('\0');
    const int descriptor = mkstemp(mutable_pattern.data());
    if (descriptor < 0) {
        const int error = errno;
        throw UserError("cannot reserve temporary for " + destination + ": " +
                        std::strerror(error));
    }
    const std::string result(mutable_pattern.data());
    const mode_t mask = umask(0);
    (void)umask(mask);
    const mode_t requested = static_cast<mode_t>(0666) & ~mask;
    if (fchmod(descriptor, requested) != 0) {
        const int error = errno;
        (void)::close(descriptor);
        (void)::unlink(result.c_str());
        throw UserError("cannot set permissions on temporary " + result + ": " +
                        std::strerror(error));
    }
    if (::close(descriptor) != 0) {
        const int error = errno;
        (void)::unlink(result.c_str());
        throw UserError("cannot close reserved temporary " + result + ": " +
                        std::strerror(error));
    }
    return result;
}

bool has_temporary_sibling(const std::string& destination) {
    const fs::path path(destination);
    const fs::path parent = path.parent_path();
    const std::string prefix = path.filename().string() + ".tmp.";
    for (const auto& entry : fs::directory_iterator(parent)) {
        if (entry.path().filename().string().rfind(prefix, 0) == 0) return true;
    }
    return false;
}

class TemporaryFileGuard {
public:
    explicit TemporaryFileGuard(std::string path) : path_(std::move(path)) {}
    ~TemporaryFileGuard() {
        if (!path_.empty()) std::remove(path_.c_str());
    }
    void release() { path_.clear(); }

private:
    std::string path_;
};

class GzipOutput {
public:
    explicit GzipOutput(const std::string& path) : path_(path) {
        handle_ = gzopen(path.c_str(), "wb");
        if (handle_ == nullptr) throw UserError("cannot open gzip output: " + path);
    }
    ~GzipOutput() {
        if (handle_ != nullptr) gzclose(handle_);
    }
    gzFile get() const { return handle_; }
    void close_checked() {
        gzFile handle = handle_;
        handle_ = nullptr;
        if (gzclose(handle) != Z_OK) {
            throw UserError("failed closing gzip output: " + path_);
        }
    }

private:
    std::string path_;
    gzFile handle_ = nullptr;
};

void require_gzprintf(gzFile output, int written, const std::string& path) {
    if (written > 0) return;
    int zlib_error = Z_OK;
    const char* detail = gzerror(output, &zlib_error);
    throw UserError("gzip write failure for " + path + ": " +
                    (detail == nullptr ? "unknown zlib error" : detail));
}

void require_input(const std::string& label, const std::string& path) {
    if (path.empty() || !fs::is_regular_file(path) || fs::file_size(path) == 0) {
        throw UserError(label + " is missing or empty: " + path);
    }
}

void require_output_parent(const std::string& path) {
    fs::path parent = fs::path(path).parent_path();
    if (parent.empty() || !fs::is_directory(parent)) {
        throw UserError("output parent directory is missing: " + parent.string());
    }
}

std::unordered_map<std::string, size_t> header_index(
        const std::vector<std::string>& fields, const std::string& path) {
    std::unordered_map<std::string, size_t> result;
    for (size_t index = 0; index < fields.size(); ++index) {
        if (fields[index].empty() || !result.emplace(fields[index], index).second) {
            throw UserError("invalid/duplicate header column in " + path);
        }
    }
    return result;
}

size_t required_column(const std::unordered_map<std::string, size_t>& index,
                       const std::string& name, const std::string& path) {
    const auto found = index.find(name);
    if (found == index.end()) throw UserError(path + ": missing column " + name);
    return found->second;
}

std::vector<std::string> load_samples(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw UserError("cannot open samples: " + path);
    std::vector<std::string> samples;
    std::unordered_set<std::string> seen;
    std::string line;
    while (std::getline(input, line)) {
        const std::string value = trim(line);
        if (value.empty()) continue;
        if (!seen.insert(value).second) throw UserError("duplicate sample: " + value);
        samples.push_back(value);
    }
    if (input.bad()) throw UserError("I/O error while reading samples: " + path);
    if (samples.empty()) throw UserError("samples file is empty: " + path);
    return samples;
}

void validate_samples_against_panel(
        const std::string& path, const std::vector<std::string>& samples) {
    htsFile* input = bcf_open(path.c_str(), "r");
    if (input == nullptr) throw UserError("cannot open main panel BCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(input);
    if (header == nullptr) {
        bcf_close(input);
        throw UserError("cannot read main panel BCF header: " + path);
    }
    std::string problem;
    const int count = bcf_hdr_nsamples(header);
    if (count != static_cast<int>(samples.size())) {
        problem = "sample count differs between .samples and main panel BCF";
    } else {
        for (int index = 0; index < count; ++index) {
            const std::string observed = header->samples[index] == nullptr
                ? "" : header->samples[index];
            if (observed != samples[static_cast<size_t>(index)]) {
                problem = "sample order differs at column " +
                    std::to_string(index + 1) + ": .samples=" +
                    samples[static_cast<size_t>(index)] + ", panel=" + observed;
                break;
            }
        }
    }
    bcf_hdr_destroy(header);
    if (bcf_close(input) != 0 && problem.empty()) {
        problem = "failed closing main panel BCF";
    }
    if (!problem.empty()) throw UserError(problem + ": " + path);
}

std::string canonical_contig(const std::string& value) {
    std::string contig = trim(value);
    std::transform(contig.begin(), contig.end(), contig.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::tolower(character));
        });
    if (contig.rfind("chr", 0) == 0) contig.erase(0, 3);
    if (contig == "m" || contig == "mt") return "mt";
    return contig;
}

bool is_sex_contig(const std::string& value) {
    const std::string contig = canonical_contig(value);
    return contig == "x" || contig == "y";
}

std::string chromosome_from_arm(const std::string& arm) {
    std::string chromosome = trim(arm);
    if (!chromosome.empty() &&
            (chromosome.back() == 'p' || chromosome.back() == 'q')) {
        chromosome.pop_back();
    }
    return chromosome;
}

void load_arms(const std::string& path, std::vector<Arm>& arms,
               std::unordered_map<std::string, std::vector<Interval>>& by_chrom,
               bool include_sex_chromosomes) {
    std::ifstream input(path);
    if (!input) throw UserError("cannot open arms BED: " + path);
    std::string line;
    long line_number = 0;
    std::unordered_map<std::string, int> name_to_index;
    std::unordered_map<int, std::string> first_interval_contig;
    std::unordered_set<int> fragmented_arms;
    std::unordered_map<std::string, std::string> raw_contigs;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> fields = split_tabs(line);
        if (fields.size() < 4) {
            throw UserError(path + ":" + std::to_string(line_number) +
                            ": expected BED4");
        }
        const std::string chrom = trim(fields[0]);
        const std::string contig_key = canonical_contig(chrom);
        const std::string name = trim(fields[3]);
        const long long start = parse_ll(trim(fields[1]), path);
        const long long end = parse_ll(trim(fields[2]), path);
        if (start < 0 || end <= start || chrom.empty() || contig_key.empty() ||
                name.empty()) {
            throw UserError(path + ": invalid BED interval on line " +
                            std::to_string(line_number));
        }
        const auto raw_inserted = raw_contigs.emplace(contig_key, chrom);
        if (!raw_inserted.second && raw_inserted.first->second != chrom) {
            throw UserError("ambiguous chromosome aliases in arms BED: " +
                            raw_inserted.first->second + " and " + chrom);
        }
        const std::string biological_chromosome = chromosome_from_arm(name);
        if (!include_sex_chromosomes && is_sex_contig(biological_chromosome)) continue;
        int arm_index = -1;
        const auto existing = name_to_index.find(name);
        if (existing == name_to_index.end()) {
            arm_index = static_cast<int>(arms.size());
            name_to_index.emplace(name, arm_index);
            first_interval_contig.emplace(arm_index, contig_key);
            arms.push_back(Arm{biological_chromosome, name, start, end});
        } else {
            arm_index = existing->second;
            Arm& arm = arms[arm_index];
            if (first_interval_contig.at(arm_index) == contig_key &&
                    fragmented_arms.count(arm_index) == 0) {
                arm.start = std::min<int64_t>(arm.start, start);
                arm.end = std::max<int64_t>(arm.end, end);
            } else {
                // A biological arm projected onto multiple ancestral contigs has
                // no single meaningful coordinate span.  Preserve the logical
                // chromosome/arm identity and mark the bounds as fragmented.
                arm.start = 0;
                arm.end = 0;
                fragmented_arms.insert(arm_index);
            }
        }
        by_chrom[contig_key].push_back(Interval{start, end, arm_index});
    }
    if (input.bad()) throw UserError("I/O error while reading arms BED: " + path);
    if (arms.empty()) {
        throw UserError("arms BED contains no included intervals: " + path +
                        (include_sex_chromosomes
                            ? "" : " (sex chromosomes are excluded by default)"));
    }
    for (auto& item : by_chrom) {
        auto& intervals = item.second;
        std::sort(intervals.begin(), intervals.end(),
                  [](const Interval& a, const Interval& b) { return a.start < b.start; });
        std::vector<Interval> merged;
        merged.reserve(intervals.size());
        for (const Interval& interval : intervals) {
            if (!merged.empty() && interval.start < merged.back().end) {
                if (interval.arm_index != merged.back().arm_index) {
                    throw UserError("overlapping intervals assigned to different arms on " +
                                    item.first);
                }
                merged.back().end = std::max(merged.back().end, interval.end);
            } else if (!merged.empty() && interval.start == merged.back().end &&
                       interval.arm_index == merged.back().arm_index) {
                merged.back().end = interval.end;
            } else {
                merged.push_back(interval);
            }
        }
        intervals.swap(merged);
    }
}

int find_arm(const std::unordered_map<std::string, std::vector<Interval>>& by_chrom,
             const std::string& chrom, int pos) {
    const auto found = by_chrom.find(canonical_contig(chrom));
    if (found == by_chrom.end()) return -1;
    const auto& values = found->second;
    auto candidate = std::upper_bound(
        values.begin(), values.end(), static_cast<int64_t>(pos),
        [](int64_t value, const Interval& interval) { return value < interval.start; });
    if (candidate == values.begin()) return -1;
    --candidate;
    return pos >= candidate->start && pos < candidate->end
        ? candidate->arm_index : -1;
}

void load_cells(const Options& options, const std::vector<std::string>& samples,
                std::unordered_map<unsigned long, Cell>& cells,
                std::unordered_map<unsigned long, std::string>& manifest_barcodes) {
    std::unordered_map<std::string, int> sample_index;
    for (size_t index = 0; index < samples.size(); ++index) {
        sample_index[samples[index]] = static_cast<int>(index);
    }
    gzFile input = gzopen(options.cells.c_str(), "rb");
    if (!input) throw UserError("cannot open cell manifest: " + options.cells);
    std::vector<char> buffer(LINE_BUFFER);
    std::string header_line;
    if (!gz_getline_checked(input, buffer, header_line, options.cells)) {
        gzclose(input);
        throw UserError("empty cell manifest: " + options.cells);
    }
    while (!header_line.empty() && (header_line.back() == '\n' || header_line.back() == '\r'))
        header_line.pop_back();
    const std::vector<std::string> header = split_tabs(header_line);
    const auto index = header_index(header, options.cells);
    const size_t library_col = required_column(index, "library", options.cells);
    const size_t barcode_col = required_column(index, "barcode", options.cells);
    const size_t donor_a_col = required_column(index, "donor_a", options.cells);
    const size_t donor_b_col = required_column(index, "donor_b", options.cells);
    const size_t pair_col = required_column(index, "donor_pair", options.cells);
    const size_t c_col = required_column(index, "ambient_c", options.cells);
    const size_t cse_col = required_column(index, "ambient_c_se", options.cells);
    const size_t eligible_col = required_column(index, "model_eligible", options.cells);
    const size_t schema_col = required_column(index, "schema_version", options.cells);
    long line_number = 1;
    std::string line;
    while (gz_getline_checked(input, buffer, line, options.cells)) {
        ++line_number;
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_tabs(line);
        if (fields.size() != header.size()) {
            gzclose(input);
            throw UserError(options.cells + ": malformed line " + std::to_string(line_number));
        }
        const std::string row_context = options.cells + ":" +
            std::to_string(line_number);
        if (trim(fields[library_col]) != options.library) {
            gzclose(input);
            throw UserError(row_context + ": library does not match --library " +
                            options.library);
        }
        if (trim(fields[schema_col]) != CELL_MANIFEST_SCHEMA) {
            gzclose(input);
            throw UserError(row_context + ": expected schema_version " +
                            CELL_MANIFEST_SCHEMA);
        }
        Cell cell;
        cell.barcode = trim(fields[barcode_col]);
        const unsigned long encoded = encode_barcode(cell.barcode, row_context);
        const auto barcode_inserted = manifest_barcodes.emplace(encoded, cell.barcode);
        if (!barcode_inserted.second) {
            gzclose(input);
            if (barcode_inserted.first->second == cell.barcode) {
                throw UserError("duplicate barcode in cell manifest: " + cell.barcode);
            }
            throw UserError("encoded barcode collision in cell manifest: " +
                            barcode_inserted.first->second + " and " + cell.barcode);
        }
        cell.donor_a = trim(fields[donor_a_col]);
        cell.donor_b = trim(fields[donor_b_col]);
        cell.donor_pair = trim(fields[pair_col]);
        cell.model_eligible = parse_bool(fields[eligible_col], row_context);
        // Validate all manifest rows, including rows that are not heterotypic
        // targets. A header-only target set is a valid terminal state, but it
        // must not turn malformed upstream data into a successful run.
        if (!parse_double(fields[c_col], cell.ambient_c) ||
                cell.ambient_c < 0.0 || cell.ambient_c >= 1.0) {
            gzclose(input);
            throw UserError("invalid ambient_c for " + cell.barcode);
        }
        if (missing_value(fields[cse_col])) {
            cell.ambient_c_se = std::numeric_limits<double>::quiet_NaN();
        } else if (!parse_double(fields[cse_col], cell.ambient_c_se) ||
                   cell.ambient_c_se < 0.0) {
            gzclose(input);
            throw UserError("invalid ambient_c_se for " + cell.barcode);
        }
        const bool donor_a_missing = missing_value(cell.donor_a);
        const bool donor_b_missing = missing_value(cell.donor_b);
        const auto a = donor_a_missing ? sample_index.end()
                                       : sample_index.find(cell.donor_a);
        const auto b = donor_b_missing ? sample_index.end()
                                       : sample_index.find(cell.donor_b);
        if ((!donor_a_missing && a == sample_index.end()) ||
                (!donor_b_missing && b == sample_index.end())) {
            gzclose(input);
            throw UserError("cell donor is absent from .samples for " + cell.barcode +
                            ": " + cell.donor_a + "+" + cell.donor_b);
        }
        if (donor_a_missing || donor_b_missing) continue;
        const std::string pair_ab = cell.donor_a + "+" + cell.donor_b;
        const std::string pair_ba = cell.donor_b + "+" + cell.donor_a;
        if (cell.donor_pair != pair_ab && cell.donor_pair != pair_ba) {
            gzclose(input);
            throw UserError("donor_pair does not match donor_a/donor_b for " +
                            cell.barcode);
        }
        if (cell.donor_a == cell.donor_b) continue;
        cell.donor_a_index = a->second;
        cell.donor_b_index = b->second;
        const auto inserted = cells.emplace(encoded, cell);
        if (!inserted.second) {
            gzclose(input);
            throw UserError("encoded barcode collision in cell manifest: " + cell.barcode);
        }
    }
    if (gzclose(input) != Z_OK) throw UserError("failed closing cell manifest");
    if (manifest_barcodes.empty()) {
        throw UserError("cell manifest contains no data rows: " + options.cells);
    }
}

void load_ambient(const Options& options, const std::vector<std::string>& samples,
                  const std::unordered_map<unsigned long, std::string>& manifest_barcodes,
                  std::unordered_map<unsigned long, Cell>& cells) {
    std::unordered_map<std::string, int> sample_index;
    for (size_t index = 0; index < samples.size(); ++index)
        sample_index[samples[index]] = static_cast<int>(index);
    gzFile input = gzopen(options.ambient_sources.c_str(), "rb");
    if (!input) throw UserError("cannot open ambient sources: " + options.ambient_sources);
    std::vector<char> buffer(LINE_BUFFER);
    std::string header_line;
    if (!gz_getline_checked(input, buffer, header_line, options.ambient_sources)) {
        gzclose(input);
        throw UserError("empty ambient source table");
    }
    while (!header_line.empty() && (header_line.back() == '\n' || header_line.back() == '\r'))
        header_line.pop_back();
    const auto header = split_tabs(header_line);
    const auto index = header_index(header, options.ambient_sources);
    const size_t library_col = required_column(index, "library", options.ambient_sources);
    const size_t barcode_col = required_column(index, "barcode", options.ambient_sources);
    const size_t source_col = required_column(index, "source_label", options.ambient_sources);
    const size_t mass_col = required_column(index, "scoring_profile_mass", options.ambient_sources);
    const size_t schema_col = required_column(index, "schema_version", options.ambient_sources);
    std::unordered_map<unsigned long, std::unordered_set<std::string>> seen;
    std::unordered_map<unsigned long, std::string> ambient_barcodes;
    std::unordered_map<unsigned long, double> total_mass;
    long line_number = 1;
    std::string line;
    while (gz_getline_checked(input, buffer, line, options.ambient_sources)) {
        ++line_number;
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != header.size()) {
            gzclose(input);
            throw UserError(options.ambient_sources + ": malformed line " +
                            std::to_string(line_number));
        }
        const std::string row_context = options.ambient_sources + ":" +
            std::to_string(line_number);
        if (trim(fields[library_col]) != options.library) {
            gzclose(input);
            throw UserError(row_context + ": library does not match --library " +
                            options.library);
        }
        if (trim(fields[schema_col]) != AMBIENT_SOURCES_SCHEMA) {
            gzclose(input);
            throw UserError(row_context + ": expected schema_version " +
                            AMBIENT_SOURCES_SCHEMA);
        }
        const std::string barcode = trim(fields[barcode_col]);
        const unsigned long encoded = encode_barcode(
            barcode, row_context);
        const auto manifest = manifest_barcodes.find(encoded);
        if (manifest == manifest_barcodes.end()) {
            gzclose(input);
            throw UserError(row_context + ": ambient barcode is absent from cell manifest: " +
                            barcode);
        }
        if (manifest->second != barcode) {
            gzclose(input);
            throw UserError("encoded barcode collision between cell manifest and ambient "
                            "source table: " + manifest->second + " and " + barcode);
        }
        const auto barcode_inserted = ambient_barcodes.emplace(encoded, barcode);
        if (!barcode_inserted.second && barcode_inserted.first->second != barcode) {
            gzclose(input);
            throw UserError("encoded barcode collision in ambient source table: " +
                            barcode_inserted.first->second + " and " + barcode);
        }
        double mass = 0.0;
        if (!parse_double(fields[mass_col], mass) || mass < 0.0) {
            gzclose(input);
            throw UserError("invalid ambient source mass on line " +
                            std::to_string(line_number));
        }
        const std::string source_label = trim(fields[source_col]);
        if (source_label.empty()) {
            gzclose(input);
            throw UserError("empty ambient source label on line " +
                            std::to_string(line_number));
        }
        if (!seen[encoded].insert(source_label).second) {
            gzclose(input);
            throw UserError("duplicate cell/source ambient row for " + barcode +
                            "/" + source_label);
        }
        total_mass[encoded] += mass;
        if (!std::isfinite(total_mass[encoded])) {
            gzclose(input);
            throw UserError("ambient source mass overflow for " + barcode);
        }
        auto cell = cells.find(encoded);
        if (cell == cells.end()) continue;
        if (cell->second.barcode != barcode) {
            gzclose(input);
            throw UserError("encoded barcode collision between cell manifest and ambient "
                            "source table: " + cell->second.barcode + " and " + barcode);
        }
        const auto source = sample_index.find(source_label);
        if (source == sample_index.end()) {
            cell->second.ambient_unmapped_mass += mass;
            continue;
        }
        cell->second.ambient.push_back(AmbientSource{source->second, mass});
    }
    if (gzclose(input) != Z_OK) throw UserError("failed closing ambient source table");
    for (const auto& item : manifest_barcodes) {
        const auto total = total_mass.find(item.first);
        if (total == total_mass.end()) {
            throw UserError("cell manifest barcode has no ambient source rows: " +
                            item.second);
        }
        if (std::fabs(total->second - 1.0) > 1e-6) {
            throw UserError("ambient source masses do not sum to one for " +
                            item.second + ": " + std::to_string(total->second));
        }
    }
    for (auto& item : cells) {
        double total = item.second.ambient_unmapped_mass;
        for (const auto& source : item.second.ambient) total += source.mass;
        if (std::fabs(total - 1.0) > 1e-6) {
            throw UserError("ambient source masses do not sum to one for " +
                            item.second.barcode + ": " + std::to_string(total));
        }
        // Normalize away harmless decimal serialization drift after enforcing
        // the simplex contract, so the ungenotyped remainder is well-defined.
        for (auto& source : item.second.ambient) source.mass /= total;
        item.second.ambient_unmapped_mass /= total;
    }
}

struct EvidenceLayout {
    bool molecules = true;
    std::string path;
    size_t barcode_col = 0;
    size_t molecule_col = 1;
    size_t basis_col = 2;
    size_t tid_col = 3;
    size_t pos_col = 4;
    size_t ref_col = 5;
    size_t alt_col = 6;
};

bool evidence_file_has_nonblank_row(const std::string& path) {
    if (path.empty() || !fs::is_regular_file(path) || fs::file_size(path) == 0) {
        return false;
    }
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input) throw UserError("cannot open pileup evidence: " + path);
    std::vector<char> buffer(LINE_BUFFER);
    std::string line;
    bool found = false;
    try {
        while (gz_getline_checked(input, buffer, line, path)) {
            if (!trim(line).empty()) {
                found = true;
                break;
            }
        }
    } catch (...) {
        (void)gzclose(input);
        throw;
    }
    if (gzclose(input) != Z_OK) {
        throw UserError("failed closing pileup evidence: " + path);
    }
    return found;
}

EvidenceLayout select_evidence(const Options& options) {
    // A valid empty gzip stream still has a nonzero byte size. Inspect its
    // logical rows so an empty molecule sidecar does not mask a usable
    // site-level fallback.
    if (evidence_file_has_nonblank_row(options.pileup_molecules)) {
        EvidenceLayout result;
        result.path = options.pileup_molecules;
        return result;
    }
    if (evidence_file_has_nonblank_row(options.pileup_observations)) {
        EvidenceLayout result;
        result.molecules = false;
        result.path = options.pileup_observations;
        result.molecule_col = 0;
        result.basis_col = 0;
        result.tid_col = 1;
        result.pos_col = 2;
        result.ref_col = 3;
        result.alt_col = 4;
        return result;
    }
    throw UserError("neither molecule-aware nor site-level pileup evidence is available");
}

std::unordered_set<uint64_t> discover_observed_sites(
        const EvidenceLayout& layout,
        const std::unordered_map<unsigned long, Cell>& cells,
        uint64_t& rows_scanned, uint64_t& rows_selected) {
    gzFile input = gzopen(layout.path.c_str(), "rb");
    if (!input) throw UserError("cannot open pileup evidence: " + layout.path);
    std::vector<char> buffer(LINE_BUFFER);
    std::unordered_set<uint64_t> sites;
    std::string line;
    const size_t expected_fields = layout.molecules ? 7U : 5U;
    while (gz_getline_checked(input, buffer, line, layout.path)) {
        ++rows_scanned;
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != expected_fields) {
            gzclose(input);
            throw UserError(layout.path + ": expected " +
                            std::to_string(expected_fields) + " tab-separated fields");
        }
        const unsigned long barcode = parse_encoded_barcode(
            fields[layout.barcode_col], layout.path);
        if (cells.find(barcode) == cells.end()) continue;
        const int tid = parse_nonnegative_int(fields[layout.tid_col], layout.path);
        const int pos = parse_nonnegative_int(fields[layout.pos_col], layout.path);
        double ref = 0.0, alt = 0.0;
        if (!parse_double(fields[layout.ref_col], ref) ||
                !parse_double(fields[layout.alt_col], alt) || ref < 0.0 || alt < 0.0 ||
                !std::isfinite(ref + alt) || ref + alt <= 0.0) {
            gzclose(input);
            throw UserError(layout.path + ": invalid ref/alt evidence values");
        }
        sites.insert(site_key(tid, pos));
        ++rows_selected;
    }
    if (gzclose(input) != Z_OK) throw UserError("failed closing pileup evidence");
    return sites;
}

std::unordered_map<uint64_t, Site> load_sites(
        const Options& options, size_t n_samples,
        const std::unordered_set<uint64_t>& observed,
        const std::unordered_map<std::string, std::vector<Interval>>& arm_intervals,
        uint64_t& site_rows_scanned, uint64_t& observed_site_keys_found,
        uint64_t& observed_sites_in_arms) {
    gzFile input = gzopen(options.pileup_sites.c_str(), "rb");
    if (!input) throw UserError("cannot open pileup sites: " + options.pileup_sites);
    std::vector<char> buffer(LINE_BUFFER);
    std::unordered_map<uint64_t, Site> result;
    std::unordered_map<int, std::string> contig_by_tid;
    std::unordered_map<std::string, int> tid_by_contig;
    std::unordered_set<uint64_t> found_observed;
    result.reserve(observed.size());
    found_observed.reserve(observed.size());
    std::string line;
    const size_t expected_fields = 5U + n_samples;
    while (gz_getline_checked(input, buffer, line, options.pileup_sites)) {
        ++site_rows_scanned;
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != expected_fields) {
            gzclose(input);
            throw UserError(options.pileup_sites + ": expected " +
                            std::to_string(expected_fields) + " tab-separated fields");
        }
        const int tid = parse_nonnegative_int(fields[0], options.pileup_sites);
        const int pos = parse_nonnegative_int(fields[2], options.pileup_sites);
        const std::string contig_key = canonical_contig(fields[1]);
        if (contig_key.empty()) {
            gzclose(input);
            throw UserError(options.pileup_sites + ": empty chromosome name");
        }
        const auto tid_contig = contig_by_tid.emplace(tid, contig_key);
        if (!tid_contig.second && tid_contig.first->second != contig_key) {
            gzclose(input);
            throw UserError(options.pileup_sites +
                            ": one TID maps to multiple chromosome aliases");
        }
        const auto contig_tid = tid_by_contig.emplace(contig_key, tid);
        if (!contig_tid.second && contig_tid.first->second != tid) {
            gzclose(input);
            throw UserError(options.pileup_sites +
                            ": ambiguous duplicate chromosome alias maps to multiple TIDs");
        }
        const uint64_t key = site_key(tid, pos);
        if (observed.find(key) == observed.end()) continue;
        if (!found_observed.insert(key).second) {
            gzclose(input);
            throw UserError("duplicate observed pileup site key in " +
                            options.pileup_sites + ": tid=" + std::to_string(tid) +
                            " pos=" + std::to_string(pos));
        }
        // vcf_hts.cpp writes bcf1_t::pos, which is zero-based. It therefore
        // enters the BED half-open lookup directly with no +/-1 conversion.
        const int arm = find_arm(arm_intervals, trim(fields[1]), pos);
        if (arm < 0) continue;
        const std::string ref = trim(fields[3]);
        const std::string alt = trim(fields[4]);
        const auto is_base = [](const std::string& allele) {
            return allele.size() == 1U && (allele[0] == 'A' || allele[0] == 'C' ||
                allele[0] == 'G' || allele[0] == 'T');
        };
        if (!is_base(ref) || !is_base(alt) || ref == alt) {
            gzclose(input);
            throw UserError(options.pileup_sites +
                            ": pileup site is not a distinct biallelic A/C/G/T SNP");
        }
        Site site;
        site.tid = tid;
        site.pos = pos;
        site.arm = arm;
        site.ref = ref[0];
        site.alt = alt[0];
        site.genotype.resize(n_samples, -1);
        for (size_t sample = 0; sample < n_samples; ++sample) {
            const long long genotype = parse_ll(fields[5 + sample], options.pileup_sites);
            if (genotype < -1 || genotype > 2) {
                gzclose(input);
                throw UserError(options.pileup_sites +
                                ": genotype must be -1, 0, 1, or 2");
            }
            site.genotype[sample] = static_cast<int8_t>(genotype);
        }
        if (!result.emplace(key, std::move(site)).second) {
            gzclose(input);
            throw UserError("duplicate pileup site key in " + options.pileup_sites);
        }
        ++observed_sites_in_arms;
    }
    if (gzclose(input) != Z_OK) throw UserError("failed closing pileup sites");
    observed_site_keys_found = static_cast<uint64_t>(found_observed.size());
    if (found_observed.size() != observed.size()) {
        uint64_t missing_key = 0;
        for (const uint64_t key : observed) {
            if (found_observed.find(key) == found_observed.end()) {
                missing_key = key;
                break;
            }
        }
        const uint32_t missing_tid = static_cast<uint32_t>(missing_key >> 32);
        const uint32_t missing_pos = static_cast<uint32_t>(missing_key);
        throw UserError(options.pileup_sites + ": missing " +
                        std::to_string(observed.size() - found_observed.size()) +
                        " site key(s) referenced by pileup evidence; first missing tid=" +
                        std::to_string(missing_tid) + " pos=" +
                        std::to_string(missing_pos));
    }
    return result;
}

double ambient_a_probability(const Cell& cell, const Site& site,
                             Orientation orientation, double& genotyped_mass) {
    double numerator = 0.0;
    genotyped_mass = 0.0;
    for (const auto& source : cell.ambient) {
        if (source.sample < 0 || source.sample >= static_cast<int>(site.genotype.size())) continue;
        const int genotype = site.genotype[static_cast<size_t>(source.sample)];
        if (genotype < 0 || genotype > 2) continue;
        const double alt_fraction = genotype / 2.0;
        const double a_like = orientation == A_IS_REF ? 1.0 - alt_fraction : alt_fraction;
        numerator += source.mass * a_like;
        genotyped_mass += source.mass;
    }
    // Preserve the whole ambient simplex. Sources absent from the genotype
    // panel, plus panel sources missing at this site, are neutral rather than
    // being discarded and renormalizing the typed subset toward an extreme.
    genotyped_mass = std::clamp(genotyped_mass, 0.0, 1.0);
    const double ungenotyped_mass = 1.0 - genotyped_mass;
    return std::clamp(numerator + 0.5 * ungenotyped_mass, 0.0, 1.0);
}

std::unordered_map<unsigned long, std::vector<Observation>> load_observations(
        const EvidenceLayout& layout,
        const std::unordered_map<unsigned long, Cell>& cells,
        const std::unordered_map<uint64_t, Site>& sites,
        uint64_t& rows_scanned, uint64_t& rows_retained,
        uint64_t& nondiscriminating_rows) {
    gzFile input = gzopen(layout.path.c_str(), "rb");
    if (!input) throw UserError("cannot open pileup evidence: " + layout.path);
    std::vector<char> buffer(LINE_BUFFER);
    std::unordered_map<unsigned long, std::vector<Observation>> result;
    result.reserve(cells.size());
    std::string line;
    const size_t expected_fields = layout.molecules ? 7U : 5U;
    while (gz_getline_checked(input, buffer, line, layout.path)) {
        ++rows_scanned;
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != expected_fields) {
            gzclose(input);
            throw UserError(layout.path + ": expected " +
                            std::to_string(expected_fields) + " tab-separated fields");
        }
        const unsigned long barcode = parse_encoded_barcode(
            fields[layout.barcode_col], layout.path);
        const auto cell = cells.find(barcode);
        if (cell == cells.end()) continue;
        const int tid = parse_nonnegative_int(fields[layout.tid_col], layout.path);
        const int pos = parse_nonnegative_int(fields[layout.pos_col], layout.path);
        const uint64_t key = site_key(tid, pos);
        const auto site_it = sites.find(key);
        if (site_it == sites.end()) continue;
        const Site& site = site_it->second;
        const int ga = site.genotype[
            static_cast<size_t>(cell->second.donor_a_index)];
        const int gb = site.genotype[
            static_cast<size_t>(cell->second.donor_b_index)];
        if (!((ga == 0 && gb == 2) || (ga == 2 && gb == 0))) {
            ++nondiscriminating_rows;
            continue;
        }
        double ref = 0.0, alt = 0.0;
        if (!parse_double(fields[layout.ref_col], ref) ||
                !parse_double(fields[layout.alt_col], alt) || ref < 0.0 || alt < 0.0 ||
                !std::isfinite(ref + alt) || ref + alt <= 0.0) {
            gzclose(input);
            throw UserError(layout.path + ": invalid ref/alt evidence values");
        }
        Observation observation;
        observation.site_key = key;
        observation.arm = site.arm;
        observation.orientation = ga == 0 ? A_IS_REF : A_IS_ALT;
        observation.a = observation.orientation == A_IS_REF ? ref : alt;
        observation.b = observation.orientation == A_IS_REF ? alt : ref;
        if (layout.molecules) {
            observation.molecule = parse_u64(fields[layout.molecule_col], layout.path);
            const std::string basis = trim(fields[layout.basis_col]);
            if (basis != "UB_GX" && basis != "UB_GN" && basis != "QNAME_FALLBACK") {
                gzclose(input);
                throw UserError(layout.path + ": unknown molecule basis '" + basis + "'");
            }
            observation.qname_fallback = basis == "QNAME_FALLBACK";
        } else {
            // Site-level fallback: one conservative independent unit per site.
            observation.molecule = key;
            observation.qname_fallback = true;
        }
        observation.ambient_a = ambient_a_probability(
            cell->second, site, observation.orientation,
            observation.ambient_genotyped_mass);
        result[barcode].push_back(observation);
        ++rows_retained;
    }
    if (gzclose(input) != Z_OK) throw UserError("failed closing pileup evidence");
    return result;
}

CellOutput summarize_cell(unsigned long encoded, const Cell& cell,
                          std::vector<Observation>& observations,
                          size_t n_arms, double hard_allele_threshold) {
    // Preserve producer row order among duplicate keys so floating-point
    // accumulation is reproducible for a fixed sidecar.
    std::stable_sort(observations.begin(), observations.end(),
        [](const Observation& a, const Observation& b) {
            return std::tie(a.molecule, a.arm, a.site_key) <
                   std::tie(b.molecule, b.arm, b.site_key);
        });

    // Merge duplicate molecule/site rows left by independent producer threads.
    size_t write = 0;
    for (size_t read = 0; read < observations.size(); ++read) {
        if (write > 0 && observations[write - 1].molecule == observations[read].molecule &&
                observations[write - 1].arm == observations[read].arm &&
                observations[write - 1].site_key == observations[read].site_key) {
            Observation& target = observations[write - 1];
            const double target_depth = target.a + target.b;
            const double source_depth = observations[read].a + observations[read].b;
            const double combined = target_depth + source_depth;
            if (combined > 0.0) {
                target.ambient_a =
                    (target.ambient_a * target_depth + observations[read].ambient_a * source_depth) /
                    combined;
                target.ambient_genotyped_mass =
                    (target.ambient_genotyped_mass * target_depth +
                     observations[read].ambient_genotyped_mass * source_depth) / combined;
            }
            target.a += observations[read].a;
            target.b += observations[read].b;
            // A merged unit is conservatively fallback-derived if any of its
            // producer-thread fragments required the read-name key.
            target.qname_fallback = target.qname_fallback ||
                observations[read].qname_fallback;
        } else {
            if (write != read) observations[write] = observations[read];
            ++write;
        }
    }
    observations.resize(write);

    CellOutput output;
    output.encoded = encoded;
    output.arms.resize(n_arms);
    size_t begin = 0;
    while (begin < observations.size()) {
        size_t end = begin + 1;
        while (end < observations.size() &&
                observations[end].molecule == observations[begin].molecule &&
                observations[end].arm == observations[begin].arm) ++end;
        ArmEvidence& arm = output.arms[
            static_cast<size_t>(observations[begin].arm)];
        ++arm.molecules;
        double a = 0.0, b = 0.0;
        double ref_orientation_depth = 0.0, alt_orientation_depth = 0.0;
        double ambient_weighted = 0.0, ambient_genotyped_weighted = 0.0;
        double total_depth = 0.0;
        bool any_qname = false;
        for (size_t index = begin; index < end; ++index) {
            const Observation& observation = observations[index];
            const double depth = observation.a + observation.b;
            a += observation.a;
            b += observation.b;
            total_depth += depth;
            ambient_weighted += depth * observation.ambient_a;
            ambient_genotyped_weighted += depth * observation.ambient_genotyped_mass;
            if (observation.orientation == A_IS_REF) ref_orientation_depth += depth;
            else alt_orientation_depth += depth;
            any_qname = any_qname || observation.qname_fallback;
            arm.sites.insert(observation.site_key);
            ++arm.molecule_site_units;
        }
        if (any_qname) ++arm.qname;
        // These are within-molecule allele-support fractions. They are not a
        // calibrated posterior because the sidecar carries support weights,
        // not an observation likelihood model.
        const double total_support = a + b;
        const double support_a = total_support > 0.0 ? a / total_support : 0.5;
        const double support_b = total_support > 0.0 ? b / total_support : 0.5;
        arm.soft_a += support_a;
        arm.soft_b += support_b;
        const double ambient_a = total_depth > 0.0 ? ambient_weighted / total_depth : 0.5;
        const double genotyped_mass = total_depth > 0.0
            ? ambient_genotyped_weighted / total_depth : 0.0;
        const double directional_score = 2.0 * support_a - 1.0;
        const double effective_weight = directional_score * directional_score;
        const double effective_a = effective_weight * support_a;
        const double effective_b = effective_weight - effective_a;
        // Raw q is weighted over exactly the same directional pseudo-depth as
        // the primary effective evidence. It is never an expectation for the
        // threshold-selected hard calls retained below as QC.
        arm.ambient_genotyped_mass_sum += effective_weight * genotyped_mass;
        arm.ambient_genotyped_weight += effective_weight;
        enum OrientationClass : uint8_t { REF_CLASS, ALT_CLASS, MIXED_CLASS };
        OrientationClass orientation_class = MIXED_CLASS;
        if (ref_orientation_depth > 1.5 * alt_orientation_depth) {
            orientation_class = REF_CLASS;
            ++arm.molecules_ref;
            arm.soft_a_ref += support_a;
            arm.soft_b_ref += support_b;
            arm.soft_a_sumsq_ref += support_a * support_a;
            arm.effective_a_ref += effective_a;
            arm.effective_b_ref += effective_b;
            arm.effective_weight_ref += effective_weight;
            arm.ambient_a_ref_sum += effective_weight * ambient_a;
            arm.ambient_ref_weight += effective_weight;
        } else if (alt_orientation_depth > 1.5 * ref_orientation_depth) {
            orientation_class = ALT_CLASS;
            ++arm.molecules_alt;
            arm.soft_a_alt += support_a;
            arm.soft_b_alt += support_b;
            arm.soft_a_sumsq_alt += support_a * support_a;
            arm.effective_a_alt += effective_a;
            arm.effective_b_alt += effective_b;
            arm.effective_weight_alt += effective_weight;
            arm.ambient_a_alt_sum += effective_weight * ambient_a;
            arm.ambient_alt_weight += effective_weight;
        } else {
            ++arm.molecules_mixed;
            arm.soft_a_mixed += support_a;
            arm.soft_b_mixed += support_b;
            arm.soft_a_sumsq_mixed += support_a * support_a;
            arm.effective_a_mixed += effective_a;
            arm.effective_b_mixed += effective_b;
            arm.effective_weight_mixed += effective_weight;
            arm.ambient_a_mixed_sum += effective_weight * ambient_a;
            arm.ambient_mixed_weight += effective_weight;
        }
        const bool hard_a = support_a >= hard_allele_threshold;
        const bool hard_b = support_b >= hard_allele_threshold;
        if (!hard_a && !hard_b) {
            ++arm.ambiguous;
            begin = end;
            continue;
        }
        ++arm.informative;
        if (orientation_class == REF_CLASS) {
            if (hard_a) ++arm.a_ref; else ++arm.b_ref;
        } else if (orientation_class == ALT_CLASS) {
            if (hard_a) ++arm.a_alt; else ++arm.b_alt;
        } else {
            if (hard_a) ++arm.a_mixed; else ++arm.b_mixed;
        }
        begin = end;
    }
    (void)cell;
    return output;
}

std::string format_double(double value) {
    if (!std::isfinite(value)) return "NA";
    std::ostringstream output;
    output << std::setprecision(17) << value;
    return output.str();
}

void publish(const std::string& temporary, const std::string& destination) {
    // Both paths are in the destination directory. link(2) gives us an
    // atomic no-replace publication primitive; rename(2) would silently
    // overwrite a result created by a concurrent/requeued array task.
    if (::link(temporary.c_str(), destination.c_str()) != 0) {
        const int error = errno;
        throw UserError("cannot publish " + destination + ": " + std::strerror(error));
    }
    if (::unlink(temporary.c_str()) != 0) {
        const int error = errno;
        const int rollback = ::unlink(destination.c_str());
        throw UserError("published " + destination + " but could not remove temporary " +
                        temporary + ": " + std::strerror(error) +
                        (rollback == 0 ? "; publication rolled back"
                                       : "; publication rollback also failed"));
    }
}

void write_outputs(const Options& options, const EvidenceLayout& layout,
                   const std::vector<Arm>& arms,
                   const std::unordered_map<unsigned long, Cell>& cells,
                   const std::vector<CellOutput>& outputs,
                   uint64_t manifest_cells,
                   uint64_t discovery_rows, uint64_t discovery_selected,
                   uint64_t observed_site_keys_discovered,
                   uint64_t observed_site_keys_found,
                   uint64_t site_rows, uint64_t sites_in_arms,
                   uint64_t evidence_rows, uint64_t evidence_retained,
                   uint64_t nondiscriminating_rows) {
    const std::string out_tmp = make_temp_path(options.output);
    TemporaryFileGuard out_guard(out_tmp);
    const std::string qc_tmp = make_temp_path(options.qc);
    TemporaryFileGuard qc_guard(qc_tmp);
    GzipOutput compressed_output(out_tmp);
    gzFile output = compressed_output.get();
    require_gzprintf(output, gzprintf(output,
        "library\tbarcode\tdonor_a\tdonor_b\tdonor_pair\tarm\tchromosome"
        "\tarm_start\tarm_end\tambient_c\tambient_c_se\tn_sites"
        "\tn_molecules\tn_informative_molecules"
        "\tn_molecules_ref\tn_molecules_alt\tn_molecules_mixed"
        "\ta_ref\tb_ref\ta_alt\tb_alt\ta_mixed\tb_mixed\tn_ambiguous"
        "\tsoft_a\tsoft_b"
        "\tsoft_a_ref\tsoft_b_ref\tsoft_a_alt\tsoft_b_alt"
        "\tsoft_a_mixed\tsoft_b_mixed"
        "\tsoft_a_sumsq_ref\tsoft_a_sumsq_alt\tsoft_a_sumsq_mixed"
        "\teffective_a_ref\teffective_b_ref\teffective_a_alt\teffective_b_alt"
        "\teffective_a_mixed\teffective_b_mixed"
        "\teffective_weight_ref\teffective_weight_alt\teffective_weight_mixed"
        "\tambient_a_ref\tambient_a_alt\tambient_a_mixed"
        "\tambient_genotyped_mass\tqname_fallback_fraction"
        "\tmean_sites_per_molecule\tevidence_basis\tmodel_eligible"
        "\tevidence_status\tschema_version\n"), out_tmp);
    uint64_t output_rows = 0;
    uint64_t passing_rows = 0;
    uint64_t soft_molecule_units = 0;
    uint64_t hard_assigned_molecule_units = 0;
    uint64_t ambiguous_molecule_units = 0;
    double effective_directional_depth = 0.0;
    for (const auto& cell_output : outputs) {
        const auto cell_it = cells.find(cell_output.encoded);
        if (cell_it == cells.end()) continue;
        const Cell& cell = cell_it->second;
        for (size_t arm_index = 0; arm_index < arms.size(); ++arm_index) {
            const ArmEvidence& evidence = cell_output.arms[arm_index];
            if (evidence.molecules == 0) continue;
            const Arm& arm = arms[arm_index];
            const double ambient_ref = evidence.ambient_ref_weight > 0.0
                ? evidence.ambient_a_ref_sum /
                    evidence.ambient_ref_weight : NAN;
            const double ambient_alt = evidence.ambient_alt_weight > 0.0
                ? evidence.ambient_a_alt_sum /
                    evidence.ambient_alt_weight : NAN;
            const double ambient_mixed = evidence.ambient_mixed_weight > 0.0
                ? evidence.ambient_a_mixed_sum /
                    evidence.ambient_mixed_weight : NAN;
            const double genotyped_mass = evidence.ambient_genotyped_weight > 0.0
                ? evidence.ambient_genotyped_mass_sum /
                    evidence.ambient_genotyped_weight : NAN;
            const double qname_fraction = evidence.molecules
                ? static_cast<double>(evidence.qname) /
                    static_cast<double>(evidence.molecules) : NAN;
            const double sites_per_molecule = evidence.molecules
                ? static_cast<double>(evidence.molecule_site_units) /
                    static_cast<double>(evidence.molecules) : NAN;
            const double row_effective_depth =
                evidence.effective_weight_ref + evidence.effective_weight_alt +
                evidence.effective_weight_mixed;
            const std::string status = row_effective_depth > 0.0
                ? (layout.molecules ? "PASS" : "PASS_SITE_FALLBACK")
                : "NO_DIRECTIONAL_SOFT_EVIDENCE";
            if (row_effective_depth > 0.0) ++passing_rows;
            soft_molecule_units += evidence.molecules;
            hard_assigned_molecule_units += evidence.informative;
            ambiguous_molecule_units += evidence.ambiguous;
            effective_directional_depth += row_effective_depth;
            require_gzprintf(output, gzprintf(output,
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%lld\t%lld\t%s\t%s",
                options.library.c_str(), cell.barcode.c_str(), cell.donor_a.c_str(),
                cell.donor_b.c_str(), cell.donor_pair.c_str(), arm.name.c_str(),
                arm.chrom.c_str(), static_cast<long long>(arm.start),
                static_cast<long long>(arm.end), format_double(cell.ambient_c).c_str(),
                format_double(cell.ambient_c_se).c_str()), out_tmp);
            require_gzprintf(output, gzprintf(output,
                "\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu"
                "\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu",
                static_cast<unsigned long long>(evidence.sites.size()),
                static_cast<unsigned long long>(evidence.molecules),
                static_cast<unsigned long long>(evidence.informative),
                static_cast<unsigned long long>(evidence.molecules_ref),
                static_cast<unsigned long long>(evidence.molecules_alt),
                static_cast<unsigned long long>(evidence.molecules_mixed),
                static_cast<unsigned long long>(evidence.a_ref),
                static_cast<unsigned long long>(evidence.b_ref),
                static_cast<unsigned long long>(evidence.a_alt),
                static_cast<unsigned long long>(evidence.b_alt),
                static_cast<unsigned long long>(evidence.a_mixed),
                static_cast<unsigned long long>(evidence.b_mixed),
                static_cast<unsigned long long>(evidence.ambiguous)), out_tmp);
            require_gzprintf(output, gzprintf(output,
                "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                format_double(evidence.soft_a).c_str(),
                format_double(evidence.soft_b).c_str(),
                format_double(evidence.soft_a_ref).c_str(),
                format_double(evidence.soft_b_ref).c_str(),
                format_double(evidence.soft_a_alt).c_str(),
                format_double(evidence.soft_b_alt).c_str(),
                format_double(evidence.soft_a_mixed).c_str(),
                format_double(evidence.soft_b_mixed).c_str(),
                format_double(evidence.soft_a_sumsq_ref).c_str(),
                format_double(evidence.soft_a_sumsq_alt).c_str(),
                format_double(evidence.soft_a_sumsq_mixed).c_str()), out_tmp);
            require_gzprintf(output, gzprintf(output,
                "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                format_double(evidence.effective_a_ref).c_str(),
                format_double(evidence.effective_b_ref).c_str(),
                format_double(evidence.effective_a_alt).c_str(),
                format_double(evidence.effective_b_alt).c_str(),
                format_double(evidence.effective_a_mixed).c_str(),
                format_double(evidence.effective_b_mixed).c_str(),
                format_double(evidence.effective_weight_ref).c_str(),
                format_double(evidence.effective_weight_alt).c_str(),
                format_double(evidence.effective_weight_mixed).c_str()), out_tmp);
            require_gzprintf(output, gzprintf(output,
                "\t%s\t%s\t%s\t%s\t%s\t%s",
                format_double(ambient_ref).c_str(), format_double(ambient_alt).c_str(),
                format_double(ambient_mixed).c_str(), format_double(genotyped_mass).c_str(),
                format_double(qname_fraction).c_str(),
                format_double(sites_per_molecule).c_str()), out_tmp);
            require_gzprintf(output, gzprintf(output, "\t%s\t%d\t%s\t%s\n",
                layout.molecules ? "MOLECULE_UB_GENE_OR_QNAME" : "SITE_LEVEL_FALLBACK",
                cell.model_eligible ? 1 : 0, status.c_str(), SCHEMA), out_tmp);
            ++output_rows;
        }
    }
    compressed_output.close_checked();
    const bool no_heterotypic_targets = cells.empty();
    const bool no_informative_evidence = !no_heterotypic_targets && passing_rows == 0;

    std::ofstream qc(qc_tmp);
    if (!qc) throw UserError("cannot open QC output: " + qc_tmp);
    qc << "metric\tvalue\n"
       << "schema_version\t" << QC_SCHEMA << "\n"
       << "tool_version\t" << VERSION << "\n"
       << "library\t" << options.library << "\n"
       << "panel_bcf\t" << options.panel_bcf << "\n"
       << "evidence_basis\t" << (layout.molecules ? "MOLECULE" : "SITE_FALLBACK") << "\n"
       << "sex_chromosomes\t"
       << (options.include_sex_chromosomes ? "INCLUDED_BY_OPT_IN" : "EXCLUDED_DEFAULT")
       << "\n"
       << "manifest_cells_total\t" << manifest_cells << "\n"
       << "heterotypic_targets\t" << cells.size() << "\n"
       << "cells_loaded\t" << cells.size() << "\n"
       << "arms_loaded\t" << arms.size() << "\n"
       << "discovery_rows_scanned\t" << discovery_rows << "\n"
       << "discovery_rows_selected\t" << discovery_selected << "\n"
       << "observed_site_keys_discovered\t" << observed_site_keys_discovered << "\n"
       << "observed_site_keys_found\t" << observed_site_keys_found << "\n"
       << "pileup_site_rows_scanned\t" << site_rows << "\n"
       << "observed_sites_in_arms\t" << sites_in_arms << "\n"
       << "evidence_rows_scanned\t" << evidence_rows << "\n"
       << "homozygous_discordant_rows_retained\t" << evidence_retained << "\n"
       << "non_homozygous_discordant_rows\t" << nondiscriminating_rows << "\n"
       << "output_rows\t" << output_rows << "\n"
       << "passing_rows\t" << passing_rows << "\n"
       << "soft_molecule_arm_units\t" << soft_molecule_units << "\n"
       << "effective_directional_depth\t" << std::setprecision(17)
       << effective_directional_depth << "\n"
       << "hard_assigned_molecule_arm_units\t" << hard_assigned_molecule_units << "\n"
       << "ambiguous_molecule_arm_units\t" << ambiguous_molecule_units << "\n"
       << "hard_allele_threshold\t" << std::setprecision(17)
       << options.hard_allele_threshold << "\n"
       << "hard_allele_threshold_legacy_alias\t--hard-posterior\n"
       << "allele_evidence_contract\tORIENTATION_CONFIDENCE_WEIGHTED_SOFT_PRIMARY_HARD_CALLS_QC_ONLY\n"
       << "effective_weight_contract\tW_EQUALS_SQUARED_TWO_P_MINUS_ONE\n"
       << "soft_likelihood_contract\tEFFECTIVE_FRACTIONAL_QUASI_LIKELIHOOD_NOT_EXACT_BINOMIAL_PMF\n"
       << "ambient_a_contract\tRAW_FULL_SIMPLEX_Q_WEIGHTED_BY_EFFECTIVE_W_UNTYPED_AND_UNMAPPED_NEUTRAL_0.5\n"
       << "ambient_genotyped_mass_role\tQC_ONLY_NOT_RENORMALIZATION_DENOMINATOR\n"
       << "status\t"
       << (no_heterotypic_targets
               ? "PASS_NO_HETEROTYPIC_TARGETS"
               : (no_informative_evidence
                      ? "PASS_NO_INFORMATIVE_EVIDENCE" : "PASS"))
       << "\n";
    qc.close();
    if (!qc) throw UserError("failed closing QC output");
    // Publish the primary data file last. A visible data file therefore never
    // points at a run whose QC sidecar failed to publish.
    publish(qc_tmp, options.qc);
    qc_guard.release();
    try {
        publish(out_tmp, options.output);
        out_guard.release();
    } catch (...) {
        std::remove(options.qc.c_str());
        throw;
    }
}

std::string usage() {
    return std::string(PROGRAM) + " " + VERSION + "\n\n"
        "Molecule-aware donor-directed chromosome-arm ASE from existing\n"
        "CellBouncer interindividual pileup sidecars.\n\n"
        "Required:\n"
        "  --samples FILE              demux .samples (genotype column order)\n"
        "  --panel-bcf FILE            main interindividual BCF; header order is verified\n"
        "  --pileup-sites FILE         demux .pileup_sites.tsv.gz\n"
        "  --cells FILE                prepared cell_manifest.tsv.gz\n"
        "  --ambient-sources FILE      prepared ambient_sources.tsv.gz\n"
        "  --arms FILE                 BED4, 0-based half-open chromosome arms\n"
        "  --library LABEL             library label\n"
        "  --output FILE               output .tsv.gz\n"
        "  --qc FILE                   output QC TSV\n\n"
        "Evidence (at least one):\n"
        "  --pileup-molecules FILE     preferred molecule-aware sidecar\n"
        "  --pileup-observations FILE  conservative site-level fallback\n\n"
        "Options:\n"
        "  --threads INT               OpenMP cell summarization threads [8]\n"
        "  --hard-allele-threshold F   molecule A/B support-fraction threshold [0.80]\n"
        "  --hard-posterior FLOAT      deprecated alias for --hard-allele-threshold\n"
        "  --include-sex-chromosomes   opt in to chrX/chrY evidence (default: exclude)\n"
        "  --self-test                 run internal invariants and exit\n"
        "  --version                   show version\n";
}

Options parse_options(int argc, char** argv) {
    Options options;
    bool hard_threshold_set = false;
    auto need = [&](int& index, const std::string& option) -> std::string {
        if (index + 1 >= argc) throw UserError("missing value after " + option);
        return argv[++index];
    };
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--samples") options.samples = need(index, arg);
        else if (arg == "--panel-bcf") options.panel_bcf = need(index, arg);
        else if (arg == "--pileup-sites") options.pileup_sites = need(index, arg);
        else if (arg == "--pileup-molecules") options.pileup_molecules = need(index, arg);
        else if (arg == "--pileup-observations") options.pileup_observations = need(index, arg);
        else if (arg == "--cells") options.cells = need(index, arg);
        else if (arg == "--ambient-sources") options.ambient_sources = need(index, arg);
        else if (arg == "--arms") options.arms = need(index, arg);
        else if (arg == "--output") options.output = need(index, arg);
        else if (arg == "--qc") options.qc = need(index, arg);
        else if (arg == "--library") options.library = need(index, arg);
        else if (arg == "--threads") {
            options.threads = parse_nonnegative_int(need(index, arg), arg);
        } else if (arg == "--hard-allele-threshold" || arg == "--hard-posterior") {
            const std::string value = need(index, arg);
            double parsed = 0.0;
            if (!parse_double(value, parsed)) {
                throw UserError(arg + ": expected a finite floating-point value, saw '" +
                                value + "'");
            }
            if (hard_threshold_set && parsed != options.hard_allele_threshold) {
                throw UserError("conflicting --hard-allele-threshold/--hard-posterior values");
            }
            options.hard_allele_threshold = parsed;
            hard_threshold_set = true;
        }
        else if (arg == "--include-sex-chromosomes") {
            options.include_sex_chromosomes = true;
        }
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--version") options.version = true;
        else if (arg == "--help" || arg == "-h") {
            std::cout << usage();
            std::exit(0);
        } else throw UserError("unknown option: " + arg);
    }
    return options;
}

void internal_self_test() {
    std::unordered_map<std::string, std::vector<Interval>> intervals;
    intervals["1"] = {Interval{0, 100, 0}, Interval{120, 200, 1}};
    intervals["mt"] = {Interval{0, 50, 2}};
    if (find_arm(intervals, "chr1", 0) != 0 ||
            find_arm(intervals, "1", 0) != 0 ||
            find_arm(intervals, "chr1", 99) != 0 ||
            find_arm(intervals, "chr1", 100) != -1 ||
            find_arm(intervals, "chr1", 120) != 1 ||
            find_arm(intervals, "chr1", 200) != -1 ||
            find_arm(intervals, "chrM", 1) != 2 ||
            find_arm(intervals, "MT", 1) != 2 ||
            !is_sex_contig("chrX") || !is_sex_contig("Y") ||
            is_sex_contig("chrM")) {
        throw UserError("self-test failed: BED half-open/contig-alias policy");
    }
    if (site_key(1, 2) == site_key(2, 1)) {
        throw UserError("self-test failed: site key collision");
    }

    bool rejected_negative = false;
    try {
        (void)parse_u64("-1", "self-test");
    } catch (const UserError&) {
        rejected_negative = true;
    }
    if (!rejected_negative) {
        throw UserError("self-test failed: negative unsigned integer accepted");
    }

    std::string barcode_a = "AAAAAAAAAAAAAAAA";
    std::string barcode_c = "CAAAAAAAAAAAAAAA";
    const unsigned long encoded_a = encode_barcode(barcode_a, "self-test");
    const unsigned long encoded_c = encode_barcode(barcode_c, "self-test");
    if (encoded_a == encoded_c) {
        throw UserError("self-test failed: distinct valid barcodes collided");
    }
    bool rejected_barcode = false;
    try {
        (void)encode_barcode("AAAAAAAAAAAAAAAA-1", "self-test");
    } catch (const UserError&) {
        rejected_barcode = true;
    }
    if (!rejected_barcode) {
        throw UserError("self-test failed: noncanonical barcode accepted");
    }

    Cell cell;
    cell.barcode = barcode_a;
    cell.donor_a = "donorA";
    cell.donor_b = "donorB";
    cell.donor_pair = "donorA+donorB";
    cell.ambient_c = 0.1;
    cell.ambient_c_se = 0.01;
    cell.model_eligible = true;
    cell.ambient.push_back(AmbientSource{0, 0.2});
    cell.ambient_unmapped_mass = 0.8;
    Site ambient_site;
    ambient_site.genotype = {0};
    double genotyped_mass = 0.0;
    const double ambient_ref = ambient_a_probability(
        cell, ambient_site, A_IS_REF, genotyped_mass);
    const double ambient_alt = ambient_a_probability(
        cell, ambient_site, A_IS_ALT, genotyped_mass);
    if (std::fabs(ambient_ref - 0.6) > 1e-12 ||
            std::fabs(ambient_alt - 0.4) > 1e-12 ||
            std::fabs(genotyped_mass - 0.2) > 1e-12) {
        throw UserError("self-test failed: ambient ungenotyped-mass handling");
    }

    Observation first;
    first.molecule = 7;
    first.site_key = site_key(0, 10);
    first.arm = 0;
    first.a = 0.9;
    first.b = 0.1;
    first.ambient_a = 0.6;
    first.ambient_genotyped_mass = 0.2;
    first.orientation = A_IS_REF;
    Observation duplicate = first;
    duplicate.qname_fallback = true;
    Observation second;
    second.molecule = 7;
    second.site_key = site_key(0, 20);
    second.arm = 0;
    second.a = 0.8;
    second.b = 0.2;
    second.ambient_a = 0.4;
    second.ambient_genotyped_mass = 0.3;
    second.orientation = A_IS_ALT;
    std::vector<Observation> observations{first, duplicate, second};
    CellOutput summarized = summarize_cell(
        encoded_a, cell, observations, 1, 0.8);
    const ArmEvidence& evidence = summarized.arms[0];
    const double expected_soft_a = 2.6 / 3.0;
    const double expected_weight =
        (2.0 * expected_soft_a - 1.0) * (2.0 * expected_soft_a - 1.0);
    if (evidence.molecules != 1 || evidence.informative != 1 ||
            evidence.molecules_ref != 1 || evidence.molecules_alt != 0 ||
            evidence.molecules_mixed != 0 || evidence.a_ref != 1 ||
            evidence.qname != 1 ||
            evidence.sites.size() != 2 || evidence.molecule_site_units != 2 ||
            std::fabs(evidence.soft_a - expected_soft_a) > 1e-12 ||
            std::fabs(evidence.soft_a_ref - evidence.soft_a) > 1e-12 ||
            std::fabs(evidence.soft_b_ref - evidence.soft_b) > 1e-12 ||
            std::fabs(evidence.soft_a_sumsq_ref -
                      evidence.soft_a * evidence.soft_a) > 1e-12 ||
            std::fabs(evidence.effective_weight_ref - expected_weight) > 1e-12 ||
            std::fabs(evidence.effective_a_ref -
                      expected_weight * expected_soft_a) > 1e-12 ||
            std::fabs(evidence.effective_b_ref -
                      expected_weight * (1.0 - expected_soft_a)) > 1e-12 ||
            std::fabs(evidence.ambient_ref_weight - expected_weight) > 1e-12 ||
            std::fabs(evidence.ambient_a_ref_sum / evidence.ambient_ref_weight -
                      (1.6 / 3.0)) > 1e-12) {
        throw UserError("self-test failed: duplicate molecule/orientation summary");
    }

    Observation boundary_a = first;
    boundary_a.molecule = 8;
    boundary_a.site_key = site_key(0, 30);
    boundary_a.a = 4.0;
    boundary_a.b = 1.0;
    boundary_a.qname_fallback = false;
    Observation boundary_b = boundary_a;
    boundary_b.molecule = 9;
    boundary_b.site_key = site_key(0, 40);
    boundary_b.a = 1.0;
    boundary_b.b = 4.0;
    std::vector<Observation> boundary_observations{boundary_a, boundary_b};
    const CellOutput boundary_output = summarize_cell(
        encoded_a, cell, boundary_observations, 1, 0.8);
    const ArmEvidence& boundary_evidence = boundary_output.arms[0];
    if (boundary_evidence.informative != 2 || boundary_evidence.a_ref != 1 ||
            boundary_evidence.b_ref != 1) {
        throw UserError("self-test failed: symmetric hard-threshold boundary");
    }

    Observation tied = first;
    tied.molecule = 10;
    tied.site_key = site_key(0, 50);
    tied.a = 1.0;
    tied.b = 1.0;
    tied.qname_fallback = false;
    std::vector<Observation> tied_observations{tied};
    const CellOutput tied_output = summarize_cell(
        encoded_a, cell, tied_observations, 1, 0.8);
    const ArmEvidence& tied_evidence = tied_output.arms[0];
    if (tied_evidence.molecules != 1 || tied_evidence.informative != 0 ||
            tied_evidence.ambiguous != 1 || tied_evidence.effective_weight_ref != 0.0 ||
            tied_evidence.effective_a_ref != 0.0 ||
            tied_evidence.effective_b_ref != 0.0) {
        throw UserError("self-test failed: zero-weight tied molecule");
    }

    const fs::path test_root = fs::temp_directory_path() /
        (std::string("tetra_arm_ase.selftest.") +
         std::to_string(static_cast<long long>(getpid())));
    fs::remove_all(test_root);
    fs::create_directory(test_root);
    try {
        const fs::path projected_bed = test_root / "projected_arms.bed";
        {
            std::ofstream output(projected_bed);
            output << "ancestor_contig_a\t10\t50\tchr1p\n"
                   << "ancestor_contig_b\t20\t80\tchr1p\n"
                   << "ancestor_contig_b\t100\t140\tchr1q\n";
            output.close();
            if (!output) {
                throw UserError("self-test failed: cannot create projected arms BED");
            }
        }
        std::vector<Arm> projected_arms;
        std::unordered_map<std::string, std::vector<Interval>> projected_intervals;
        load_arms(projected_bed.string(), projected_arms, projected_intervals, false);
        const int contig_a_arm = find_arm(
            projected_intervals, "ancestor_contig_a", 25);
        const int contig_b_arm = find_arm(
            projected_intervals, "ancestor_contig_b", 25);
        if (projected_arms.size() != 2 || contig_a_arm < 0 ||
                contig_a_arm != contig_b_arm ||
                projected_arms.at(static_cast<size_t>(contig_a_arm)).name != "chr1p" ||
                projected_arms.at(static_cast<size_t>(contig_a_arm)).chrom != "chr1" ||
                projected_arms.at(static_cast<size_t>(contig_a_arm)).start != 0 ||
                projected_arms.at(static_cast<size_t>(contig_a_arm)).end != 0 ||
                find_arm(projected_intervals, "ancestor_contig_b", 110) == contig_a_arm) {
            throw UserError(
                "self-test failed: fragmented gene-synteny arm aggregation");
        }

        Options options;
        options.library = "selftest";
        options.output = (test_root / "evidence.tsv.gz").string();
        options.qc = (test_root / "evidence.qc.tsv").string();
        EvidenceLayout layout;
        std::vector<Arm> arms{Arm{"chr1", "chr1p", 0, 100}};
        std::unordered_map<unsigned long, Cell> cells{{encoded_a, cell}};
        std::vector<CellOutput> outputs{summarized};
        write_outputs(options, layout, arms, cells, outputs,
                      1, 3, 3, 2, 2, 2, 2, 3, 3, 0);
        gzFile input = gzopen(options.output.c_str(), "rb");
        if (input == nullptr) {
            throw UserError("self-test failed: cannot reopen gzip output");
        }
        std::vector<char> buffer(LINE_BUFFER);
        std::string header;
        std::string row;
        const bool has_header = gz_getline_checked(
            input, buffer, header, options.output);
        const bool has_row = gz_getline_checked(input, buffer, row, options.output);
        const int close_status = gzclose(input);
        while (!header.empty() && (header.back() == '\n' || header.back() == '\r'))
            header.pop_back();
        while (!row.empty() && (row.back() == '\n' || row.back() == '\r'))
            row.pop_back();
        const auto header_fields = split_tabs(header);
        const auto row_fields = split_tabs(row);
        const auto output_columns = header_index(header_fields, "self-test output");
        if (!has_header || !has_row || close_status != Z_OK ||
                header_fields.size() != 54 || row_fields.size() != header_fields.size() ||
                output_columns.find("n_molecules_ref") == output_columns.end() ||
                output_columns.find("soft_a_ref") == output_columns.end() ||
                output_columns.find("soft_b_mixed") == output_columns.end() ||
                output_columns.find("soft_a_sumsq_alt") == output_columns.end() ||
                output_columns.find("effective_a_ref") == output_columns.end() ||
                output_columns.find("effective_weight_mixed") == output_columns.end() ||
                row_fields.back() != SCHEMA || !fs::is_regular_file(options.qc) ||
                has_temporary_sibling(options.output) ||
                has_temporary_sibling(options.qc)) {
            throw UserError("self-test failed: gzip/atomic output contract");
        }

        Options tied_options;
        tied_options.library = "selftest-tied";
        tied_options.output = (test_root / "tied.tsv.gz").string();
        tied_options.qc = (test_root / "tied.qc.tsv").string();
        const std::vector<CellOutput> tied_outputs{tied_output};
        write_outputs(tied_options, layout, arms, cells, tied_outputs,
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 0);
        input = gzopen(tied_options.output.c_str(), "rb");
        if (input == nullptr) {
            throw UserError("self-test failed: cannot reopen tied output");
        }
        std::string tied_header;
        std::string tied_row;
        const bool tied_has_header = gz_getline_checked(
            input, buffer, tied_header, tied_options.output);
        const bool tied_has_row = gz_getline_checked(
            input, buffer, tied_row, tied_options.output);
        const int tied_close_status = gzclose(input);
        while (!tied_header.empty() &&
                (tied_header.back() == '\n' || tied_header.back() == '\r')) {
            tied_header.pop_back();
        }
        while (!tied_row.empty() &&
                (tied_row.back() == '\n' || tied_row.back() == '\r')) {
            tied_row.pop_back();
        }
        const auto tied_columns = header_index(
            split_tabs(tied_header), "self-test tied output");
        const auto tied_fields = split_tabs(tied_row);
        if (!tied_has_header || !tied_has_row || tied_close_status != Z_OK ||
                tied_fields.at(tied_columns.at("evidence_status")) !=
                    "NO_DIRECTIONAL_SOFT_EVIDENCE" ||
                tied_fields.at(tied_columns.at("effective_weight_ref")) != "0" ||
                tied_fields.at(tied_columns.at("ambient_a_ref")) != "NA") {
            throw UserError("self-test failed: tied evidence output contract");
        }

        bool rejected_clobber = false;
        try {
            write_outputs(options, layout, arms, cells, outputs,
                          1, 3, 3, 2, 2, 2, 2, 3, 3, 0);
        } catch (const UserError&) {
            rejected_clobber = true;
        }
        if (!rejected_clobber || !fs::is_regular_file(options.output) ||
                !fs::is_regular_file(options.qc)) {
            throw UserError("self-test failed: no-clobber publication contract");
        }

        Options blocked_options;
        blocked_options.library = "selftest-blocked";
        blocked_options.output = (test_root / "blocked.tsv.gz").string();
        blocked_options.qc = (test_root / "blocked.qc.tsv").string();
        {
            std::ofstream blocked(blocked_options.output);
            blocked << "sentinel\n";
            blocked.close();
            if (!blocked) {
                throw UserError("self-test failed: cannot create blocked destination");
            }
        }
        bool rejected_primary_clobber = false;
        try {
            write_outputs(blocked_options, layout, arms, cells, outputs,
                          1, 3, 3, 2, 2, 2, 2, 3, 3, 0);
        } catch (const UserError&) {
            rejected_primary_clobber = true;
        }
        std::ifstream blocked(blocked_options.output);
        std::string sentinel;
        std::getline(blocked, sentinel);
        if (!rejected_primary_clobber || sentinel != "sentinel" ||
                fs::exists(blocked_options.qc) ||
                has_temporary_sibling(blocked_options.output) ||
                has_temporary_sibling(blocked_options.qc)) {
            throw UserError("self-test failed: primary no-clobber/QC rollback contract");
        }

        Options empty_options;
        empty_options.library = "selftest-empty";
        empty_options.output = (test_root / "empty.tsv.gz").string();
        empty_options.qc = (test_root / "empty.qc.tsv").string();
        const std::unordered_map<unsigned long, Cell> empty_cells;
        const std::vector<CellOutput> empty_outputs;
        write_outputs(empty_options, layout, arms, empty_cells, empty_outputs,
                      1, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        input = gzopen(empty_options.output.c_str(), "rb");
        if (input == nullptr) {
            throw UserError("self-test failed: cannot reopen header-only gzip output");
        }
        std::string empty_header;
        std::string unexpected_row;
        const bool empty_has_header = gz_getline_checked(
            input, buffer, empty_header, empty_options.output);
        const bool empty_has_row = gz_getline_checked(
            input, buffer, unexpected_row, empty_options.output);
        const int empty_close_status = gzclose(input);
        std::ifstream empty_qc(empty_options.qc);
        std::ostringstream empty_qc_text;
        empty_qc_text << empty_qc.rdbuf();
        if (!empty_has_header || empty_has_row || empty_close_status != Z_OK ||
                !empty_qc || empty_qc_text.str().find(
                    "status\tPASS_NO_HETEROTYPIC_TARGETS\n") == std::string::npos ||
                has_temporary_sibling(empty_options.output) ||
                has_temporary_sibling(empty_options.qc)) {
            throw UserError("self-test failed: zero-target terminal state");
        }

        Options no_data_options;
        no_data_options.library = "selftest-no-data";
        no_data_options.output = (test_root / "no_data.tsv.gz").string();
        no_data_options.qc = (test_root / "no_data.qc.tsv").string();
        write_outputs(no_data_options, layout, arms, cells, empty_outputs,
                      1, 5, 0, 0, 0, 0, 0, 0, 0, 0);
        input = gzopen(no_data_options.output.c_str(), "rb");
        if (input == nullptr) {
            throw UserError("self-test failed: cannot reopen no-data gzip output");
        }
        std::string no_data_header;
        std::string no_data_row;
        const bool no_data_has_header = gz_getline_checked(
            input, buffer, no_data_header, no_data_options.output);
        const bool no_data_has_row = gz_getline_checked(
            input, buffer, no_data_row, no_data_options.output);
        const int no_data_close_status = gzclose(input);
        std::ifstream no_data_qc(no_data_options.qc);
        std::ostringstream no_data_qc_text;
        no_data_qc_text << no_data_qc.rdbuf();
        if (!no_data_has_header || no_data_has_row || no_data_close_status != Z_OK ||
                !no_data_qc || no_data_qc_text.str().find(
                    "status\tPASS_NO_INFORMATIVE_EVIDENCE\n") == std::string::npos ||
                has_temporary_sibling(no_data_options.output) ||
                has_temporary_sibling(no_data_options.qc)) {
            throw UserError("self-test failed: no-evidence terminal state");
        }
    } catch (...) {
        fs::remove_all(test_root);
        throw;
    }
    fs::remove_all(test_root);
    std::cout << "tetra_arm_ase self-test PASS\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        if (options.version) {
            std::cout << PROGRAM << " " << VERSION << "\n";
            return 0;
        }
        if (options.self_test) {
            internal_self_test();
            return 0;
        }
        if (options.samples.empty() || options.panel_bcf.empty() ||
                options.pileup_sites.empty() ||
                options.cells.empty() || options.ambient_sources.empty() ||
                options.arms.empty() || options.output.empty() || options.qc.empty() ||
                options.library.empty()) {
            throw UserError("missing required option\n" + usage());
        }
        if (options.threads < 1 || options.threads > 256)
            throw UserError("--threads must be between 1 and 256");
        if (!(options.hard_allele_threshold > 0.5 &&
                options.hard_allele_threshold < 1.0))
            throw UserError("--hard-allele-threshold must be between 0.5 and 1");
        require_input("samples", options.samples);
        require_input("main panel BCF", options.panel_bcf);
        require_input("pileup sites", options.pileup_sites);
        require_input("cell manifest", options.cells);
        require_input("ambient sources", options.ambient_sources);
        require_input("arms BED", options.arms);
        require_output_parent(options.output);
        require_output_parent(options.qc);
        if (!fs::path(options.output).is_absolute() || !fs::path(options.qc).is_absolute())
            throw UserError("output paths must be absolute");
        if (fs::path(options.output).lexically_normal() ==
                fs::path(options.qc).lexically_normal()) {
            throw UserError("--output and --qc must be different paths");
        }
        if (fs::path(options.output).extension() != ".gz") {
            throw UserError("--output must end in .gz because ASE evidence is gzip-compressed");
        }
        if (fs::exists(options.output) || fs::exists(options.qc))
            throw UserError("output already exists; select a new run or remove it explicitly");

        const EvidenceLayout evidence = select_evidence(options);
        const auto samples = load_samples(options.samples);
        validate_samples_against_panel(options.panel_bcf, samples);
        std::vector<Arm> arms;
        std::unordered_map<std::string, std::vector<Interval>> arm_intervals;
        load_arms(options.arms, arms, arm_intervals,
                  options.include_sex_chromosomes);
        std::unordered_map<unsigned long, Cell> cells;
        std::unordered_map<unsigned long, std::string> manifest_barcodes;
        load_cells(options, samples, cells, manifest_barcodes);
        load_ambient(options, samples, manifest_barcodes, cells);

        if (cells.empty()) {
            // A library without heterotypic two-donor targets is expected in
            // sparse arrays. Publish the complete schema plus an explicit QC
            // terminal state; do not inspect irrelevant evidence rows.
            const std::vector<CellOutput> outputs;
            write_outputs(options, evidence, arms, cells, outputs,
                          static_cast<uint64_t>(manifest_barcodes.size()),
                          0, 0, 0, 0, 0, 0, 0, 0, 0);
            return 0;
        }

        uint64_t discovery_rows = 0, discovery_selected = 0;
        auto observed_sites = discover_observed_sites(
            evidence, cells, discovery_rows, discovery_selected);
        if (observed_sites.empty()) {
            const std::vector<CellOutput> outputs;
            write_outputs(options, evidence, arms, cells, outputs,
                          static_cast<uint64_t>(manifest_barcodes.size()),
                          discovery_rows, discovery_selected,
                          0, 0, 0, 0, 0, 0, 0);
            return 0;
        }
        const uint64_t observed_site_keys_discovered =
            static_cast<uint64_t>(observed_sites.size());
        uint64_t site_rows = 0, observed_site_keys_found = 0, sites_in_arms = 0;
        auto sites = load_sites(options, samples.size(), observed_sites, arm_intervals,
                                site_rows, observed_site_keys_found, sites_in_arms);
        observed_sites.clear();
        observed_sites.rehash(0);
        if (sites.empty()) {
            const std::vector<CellOutput> outputs;
            write_outputs(options, evidence, arms, cells, outputs,
                          static_cast<uint64_t>(manifest_barcodes.size()),
                          discovery_rows, discovery_selected,
                          observed_site_keys_discovered, observed_site_keys_found,
                          site_rows, sites_in_arms, 0, 0, 0);
            return 0;
        }
        uint64_t evidence_rows = 0, evidence_retained = 0, nondiscriminating_rows = 0;
        auto by_cell = load_observations(evidence, cells, sites, evidence_rows,
                                         evidence_retained, nondiscriminating_rows);

        std::vector<unsigned long> order;
        order.reserve(cells.size());
        for (const auto& item : cells) order.push_back(item.first);
        std::sort(order.begin(), order.end());
        std::vector<CellOutput> outputs(order.size());
        omp_set_num_threads(options.threads);
        #pragma omp parallel for schedule(dynamic, 1)
        for (long index = 0; index < static_cast<long>(order.size()); ++index) {
            const unsigned long encoded = order[static_cast<size_t>(index)];
            auto found = by_cell.find(encoded);
            std::vector<Observation> empty;
            std::vector<Observation>& observations = found == by_cell.end()
                ? empty : found->second;
            outputs[static_cast<size_t>(index)] = summarize_cell(
                encoded, cells.at(encoded), observations, arms.size(),
                options.hard_allele_threshold);
        }
        write_outputs(options, evidence, arms, cells, outputs,
                      static_cast<uint64_t>(manifest_barcodes.size()),
                      discovery_rows, discovery_selected,
                      observed_site_keys_discovered, observed_site_keys_found,
                      site_rows, sites_in_arms,
                      evidence_rows, evidence_retained, nondiscriminating_rows);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << "\n";
        return 2;
    }
}
