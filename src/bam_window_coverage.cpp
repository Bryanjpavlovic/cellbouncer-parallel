#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <unistd.h>

#ifndef CELLBOUNCER_SOURCE_REVISION
#define CELLBOUNCER_SOURCE_REVISION "archive"
#endif

namespace fs = std::filesystem;

namespace {

constexpr const char* kProgram = "bam_window_coverage";
constexpr const char* kVersion = "1.0.0";
constexpr const char* kSchema = "bam_window_coverage_v1";
constexpr uint64_t kDefaultWindowSize = 100;
constexpr uint32_t kDefaultExcludeFlags = 0x0F04U;
constexpr int kDefaultThreads = 4;
constexpr std::array<const char*, 5> kPopulations = {
    "bam_all", "all_barcoded", "cells", "noncell", "empty"};

class UserError : public std::runtime_error {
public:
    explicit UserError(const std::string& message) : std::runtime_error(message) {}
};

std::string version_string() {
    return std::string(kProgram) + " " + kVersion + " (source " +
           CELLBOUNCER_SOURCE_REVISION + ")";
}

std::string build_help() {
    return
        "Usage:\n"
        "  bam_window_coverage build --bam ABS --cell-barcodes ABS\n"
        "      --output ABS --metadata ABS [options]\n\n"
        "Required:\n"
        "  --bam PATH               Coordinate-sorted BAM input (absolute path)\n"
        "  --cell-barcodes PATH     Called-cell barcode roster (absolute path)\n"
        "  --output PATH            Sparse BGZF count TSV (absolute path)\n"
        "  --metadata PATH          Key/value metadata TSV (absolute path)\n\n"
        "Options:\n"
        "  --empty-barcodes PATH    Explicit empty-droplet roster (absolute path)\n"
        "  --window-size INT        Fixed window size [100]\n"
        "  --min-mapq INT           Minimum mapping quality, 0..255 [0]\n"
        "  --exclude-flags INT      SAM flag mask, decimal or 0x-prefixed [0xF04]\n"
        "  --tag TAG                Two-character barcode tag [CB]\n"
        "  --threads INT            HTSlib BGZF reader threads [4]\n"
        "  --no-normalize-10x       Do not strip terminal -<digits> suffixes\n"
        "  --normalize-10x          Explicitly retain the default normalization\n"
        "  --force                  Atomically replace existing outputs\n"
        "  -h, --help               Show this help\n\n"
        "Counts are aligned reference bases from CIGAR M, =, and X operations.\n"
        "CIGAR D and N advance reference position but never contribute coverage.\n"
        "No base-quality filter or mate-overlap collapse is applied; overlapping\n"
        "mates contribute independently.\n"
        "cells, noncell, and empty are disjoint: noncell excludes explicit empties.\n";
}

std::string merge_help() {
    return
        "Usage:\n"
        "  bam_window_coverage merge --input-list PATH --output-prefix PATH\n"
        "      --metadata PATH [options]\n\n"
        "Required:\n"
        "  --input-list PATH        TSV rows: library_id<TAB>coverage_path\n"
        "  --output-prefix PATH     Output prefix\n"
        "  --metadata PATH          Key/value merge metadata TSV\n\n"
        "Options:\n"
        "  --populations LIST       Comma-separated population columns\n"
        "                           [bam_all,all_barcoded,cells,noncell,empty]\n"
        "  --threads INT            Shared BGZF writer threads [4]\n"
        "  --force                  Atomically replace existing outputs\n"
        "  -h, --help               Show this help\n\n"
        "The input list may have the header library_id<TAB>coverage_path and may\n"
        "contain blank lines or lines beginning with '#'. Each output is\n"
        "<prefix>.<population>.bedGraph.gz and has four bedGraph columns. All\n"
        "merge paths, including coverage_path values, must be absolute.\n";
}

std::string main_help() {
    return version_string() + "\n\n" +
        "Build sparse per-BAM aligned-base window counts and stream-merge them.\n\n"
        "Usage:\n"
        "  bam_window_coverage build [options]\n"
        "  bam_window_coverage merge [options]\n"
        "  bam_window_coverage --version\n\n"
        "Run 'bam_window_coverage <subcommand> --help' for full options.\n";
}

bool has_forbidden_tsv_char(const std::string& value) {
    return value.find('\t') != std::string::npos ||
           value.find('\n') != std::string::npos ||
           value.find('\r') != std::string::npos;
}

std::string trim(const std::string& value) {
    size_t first = 0;
    while (first < value.size() &&
           std::isspace(static_cast<unsigned char>(value[first]))) {
        ++first;
    }
    size_t last = value.size();
    while (last > first &&
           std::isspace(static_cast<unsigned char>(value[last - 1]))) {
        --last;
    }
    return value.substr(first, last - first);
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

std::vector<std::string> split_commas(const std::string& value) {
    std::vector<std::string> fields;
    size_t begin = 0;
    while (true) {
        const size_t comma = value.find(',', begin);
        fields.push_back(trim(value.substr(
            begin, comma == std::string::npos ? std::string::npos : comma - begin)));
        if (comma == std::string::npos) return fields;
        begin = comma + 1;
    }
}

uint64_t parse_u64(const std::string& text, const std::string& option,
                   int base = 10) {
    if (text.empty() || text[0] == '-') {
        throw UserError(option + " requires a non-negative integer, got '" + text + "'");
    }
    size_t used = 0;
    unsigned long long parsed = 0;
    try {
        parsed = std::stoull(text, &used, base);
    } catch (const std::exception&) {
        throw UserError(option + " requires a valid integer, got '" + text + "'");
    }
    if (used != text.size()) {
        throw UserError(option + " requires a valid integer, got '" + text + "'");
    }
    return static_cast<uint64_t>(parsed);
}

uint64_t parse_u64_view(std::string_view text, const char* field,
                        const std::string& path) {
    uint64_t value = 0;
    if (text.empty() || text.front() == '-') {
        throw UserError(std::string("invalid ") + field + " in " + path);
    }
    const char* first = text.data();
    const char* last = first + text.size();
    const std::from_chars_result parsed = std::from_chars(first, last, value);
    if (parsed.ec != std::errc() || parsed.ptr != last) {
        throw UserError(std::string("invalid ") + field + " in " + path);
    }
    return value;
}

bool split_eight_tabs(std::string_view line,
                      std::array<std::string_view, 8>& fields) {
    size_t begin = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
        const size_t tab = line.find('\t', begin);
        if (i + 1 == fields.size()) {
            if (tab != std::string_view::npos) return false;
            fields[i] = line.substr(begin);
            return true;
        }
        if (tab == std::string_view::npos) return false;
        fields[i] = line.substr(begin, tab - begin);
        begin = tab + 1;
    }
    return false;
}

int parse_threads(const std::string& text) {
    const uint64_t value = parse_u64(text, "--threads");
    if (value < 1 || value > 256) {
        throw UserError("--threads must be between 1 and 256");
    }
    return static_cast<int>(value);
}

void require_absolute(const std::string& option, const std::string& path) {
    if (!fs::path(path).is_absolute()) {
        throw UserError(option + " must be an absolute path: " + path);
    }
}

void require_input_file(const std::string& option, const std::string& path) {
    std::error_code ec;
    if (!fs::exists(path, ec) || ec || !fs::is_regular_file(path, ec) || ec) {
        throw UserError(option + " is not a readable regular file: " + path);
    }
}

void require_output_parent(const std::string& option, const std::string& path) {
    const fs::path parent = fs::path(path).parent_path();
    std::error_code ec;
    if (parent.empty() || !fs::exists(parent, ec) || ec ||
        !fs::is_directory(parent, ec) || ec) {
        throw UserError(option + " parent directory does not exist: " +
                        (parent.empty() ? std::string(".") : parent.string()));
    }
}

void require_distinct_paths(const std::vector<std::pair<std::string, std::string>>& paths) {
    std::unordered_map<std::string, std::string> seen;
    for (const auto& item : paths) {
        std::error_code ec;
        const fs::path absolute = fs::absolute(fs::path(item.second), ec);
        if (ec) {
            throw UserError("cannot resolve path for " + item.first + ": " + item.second +
                            ": " + ec.message());
        }
        const std::string normalized = absolute.lexically_normal().string();
        const auto inserted = seen.emplace(normalized, item.first);
        if (!inserted.second) {
            throw UserError(item.first + " and " + inserted.first->second +
                            " resolve to the same path: " + normalized);
        }
    }
}

void require_replaceable(const std::string& path, bool force) {
    std::error_code ec;
    if (!fs::exists(path, ec)) {
        if (ec) throw UserError("cannot inspect output path " + path + ": " + ec.message());
        return;
    }
    if (ec) throw UserError("cannot inspect output path " + path + ": " + ec.message());
    if (fs::is_directory(path, ec)) {
        throw UserError("output path is a directory: " + path);
    }
    if (!force) {
        throw UserError("output already exists (use --force to replace): " + path);
    }
}

std::string temporary_path_for(const std::string& destination) {
    static uint64_t counter = 0;
    for (int attempt = 0; attempt < 100; ++attempt) {
        const std::string candidate = destination + ".tmp." +
            std::to_string(static_cast<unsigned long long>(getpid())) + "." +
            std::to_string(++counter);
        std::error_code ec;
        if (!fs::exists(candidate, ec) && !ec) return candidate;
    }
    throw UserError("could not allocate a temporary path beside " + destination);
}

void atomic_publish(const std::string& temporary, const std::string& destination) {
    if (std::rename(temporary.c_str(), destination.c_str()) != 0) {
        const int error_number = errno;
        throw UserError("cannot atomically rename " + temporary + " to " +
                        destination + ": " + std::strerror(error_number));
    }
}

void invalidate_old_metadata(const std::string& metadata, bool force) {
    if (!force) return;
    if (std::remove(metadata.c_str()) != 0 && errno != ENOENT) {
        const int error_number = errno;
        throw UserError("cannot invalidate old completion metadata " + metadata +
                        ": " + std::strerror(error_number));
    }
}

class TempFiles {
public:
    void add(const std::string& path) { paths_.push_back(path); }
    ~TempFiles() {
        for (const std::string& path : paths_) {
            std::error_code ec;
            fs::remove(path, ec);
        }
    }
private:
    std::vector<std::string> paths_;
};

std::string tsv_escape(const std::string& value) {
    std::string result;
    result.reserve(value.size());
    for (const char ch : value) {
        if (ch == '\t') result += "\\t";
        else if (ch == '\n') result += "\\n";
        else if (ch == '\r') result += "\\r";
        else result.push_back(ch);
    }
    return result;
}

void write_metadata(const std::string& temporary,
                    const std::vector<std::pair<std::string, std::string>>& values) {
    std::ofstream out(temporary, std::ios::out | std::ios::trunc);
    if (!out) throw UserError("cannot open metadata temporary file: " + temporary);
    out << "key\tvalue\n";
    for (const auto& item : values) {
        if (has_forbidden_tsv_char(item.first)) {
            throw UserError("internal error: invalid metadata key");
        }
        out << item.first << '\t' << tsv_escape(item.second) << '\n';
    }
    out.flush();
    if (!out) throw UserError("failed while writing metadata: " + temporary);
    out.close();
    if (!out) throw UserError("failed while closing metadata: " + temporary);
}

void checked_add(uint64_t& destination, uint64_t value,
                 const std::string& context) {
    if (value > std::numeric_limits<uint64_t>::max() - destination) {
        throw UserError("64-bit count overflow while " + context);
    }
    destination += value;
}

class BgzfWriter {
public:
    explicit BgzfWriter(const std::string& path) : path_(path) {
        fp_ = bgzf_open(path.c_str(), "w");
        if (fp_ == nullptr) throw UserError("cannot open BGZF output: " + path);
    }

    BgzfWriter(const BgzfWriter&) = delete;
    BgzfWriter& operator=(const BgzfWriter&) = delete;

    ~BgzfWriter() {
        if (fp_ != nullptr) bgzf_close(fp_);
    }

    void attach_pool(hts_tpool* pool) {
        if (pool != nullptr && bgzf_thread_pool(fp_, pool, 256) != 0) {
            throw UserError("cannot attach BGZF thread pool for " + path_);
        }
    }

    void write(const std::string& data) {
        const ssize_t written = bgzf_write(fp_, data.data(), data.size());
        if (written < 0 || static_cast<size_t>(written) != data.size()) {
            throw UserError("failed while writing BGZF output: " + path_);
        }
    }

    void close() {
        if (fp_ == nullptr) return;
        BGZF* closing = fp_;
        fp_ = nullptr;
        if (bgzf_close(closing) != 0) {
            throw UserError("failed while closing BGZF output: " + path_);
        }
    }

private:
    std::string path_;
    BGZF* fp_ = nullptr;
};

void read_gzip_transparent_lines(
        const std::string& path,
        const std::function<void(const std::string&, uint64_t)>& consume) {
    gzFile input = gzopen(path.c_str(), "rb");
    if (input == nullptr) throw UserError("cannot open text input: " + path);
    if (gzbuffer(input, 1U << 20) != 0) {
        gzclose(input);
        throw UserError("cannot allocate zlib input buffer for " + path);
    }

    std::array<char, 65536> buffer{};
    std::string pending;
    uint64_t line_number = 0;
    while (true) {
        char* result = gzgets(input, buffer.data(), static_cast<int>(buffer.size()));
        if (result == nullptr) break;
        pending.append(result);
        if (!pending.empty() && pending.back() == '\n') {
            pending.pop_back();
            if (!pending.empty() && pending.back() == '\r') pending.pop_back();
            consume(pending, ++line_number);
            pending.clear();
        }
        if (pending.size() > (1U << 20)) {
            gzclose(input);
            throw UserError("text line exceeds 1 MiB in " + path);
        }
    }
    if (!pending.empty()) {
        if (pending.back() == '\r') pending.pop_back();
        consume(pending, ++line_number);
    }

    int zlib_error = Z_OK;
    const char* zlib_message = gzerror(input, &zlib_error);
    if (zlib_error != Z_OK && zlib_error != Z_STREAM_END) {
        const std::string message = zlib_message == nullptr ? "unknown zlib error" : zlib_message;
        gzclose(input);
        throw UserError("failed while reading " + path + ": " + message);
    }
    if (gzclose(input) != Z_OK) {
        throw UserError("failed while closing text input: " + path);
    }
}

std::string normalize_10x(const std::string& barcode) {
    const size_t dash = barcode.rfind('-');
    if (dash == std::string::npos || dash + 1 == barcode.size()) return barcode;
    for (size_t i = dash + 1; i < barcode.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(barcode[i]))) return barcode;
    }
    return barcode.substr(0, dash);
}

struct BarcodeRoster {
    std::unordered_set<std::string> exact;
    std::unordered_set<std::string> normalized;
    uint64_t data_lines = 0;
    uint64_t duplicate_exact = 0;
    uint64_t normalization_collisions = 0;
};

BarcodeRoster load_roster(const std::string& path, bool normalize) {
    BarcodeRoster roster;
    bool saw_data = false;
    read_gzip_transparent_lines(path, [&](const std::string& raw, uint64_t line_number) {
        const std::string stripped = trim(raw);
        if (stripped.empty() || stripped[0] == '#') return;
        const size_t whitespace = stripped.find_first_of("\t ");
        const std::string barcode = stripped.substr(0, whitespace);
        if (barcode.empty()) return;
        if (!saw_data && (barcode == "barcode" || barcode == "cell_barcode")) {
            saw_data = true;
            return;
        }
        saw_data = true;
        if (barcode.find_first_of("\t\r\n ") != std::string::npos) {
            throw UserError("invalid barcode at " + path + ":" +
                            std::to_string(line_number));
        }
        ++roster.data_lines;
        if (!roster.exact.insert(barcode).second) ++roster.duplicate_exact;
        const std::string key = normalize ? normalize_10x(barcode) : barcode;
        if (!roster.normalized.insert(key).second && roster.exact.count(barcode) == 1) {
            // This includes exact duplicates as well; the final collision count is
            // recomputed below from unique exact and normalized cardinalities.
        }
    });
    if (roster.exact.empty()) {
        throw UserError("barcode roster has no data rows: " + path);
    }
    roster.normalization_collisions = roster.exact.size() - roster.normalized.size();
    return roster;
}

void ensure_rosters_disjoint(const BarcodeRoster& cells,
                             const BarcodeRoster& empties,
                             bool normalize) {
    const auto& left = normalize ? cells.normalized : cells.exact;
    const auto& right = normalize ? empties.normalized : empties.exact;
    const auto* smaller = &left;
    const auto* larger = &right;
    if (smaller->size() > larger->size()) std::swap(smaller, larger);
    size_t shown = 0;
    std::ostringstream examples;
    uint64_t overlap = 0;
    for (const std::string& barcode : *smaller) {
        if (larger->count(barcode) != 0) {
            ++overlap;
            if (shown++ < 5) {
                if (examples.tellp() > 0) examples << ", ";
                examples << barcode;
            }
        }
    }
    if (overlap != 0) {
        throw UserError("cell and empty barcode rosters overlap after " +
                        std::string(normalize ? "10x normalization" : "exact matching") +
                        " (" + std::to_string(overlap) + " keys; examples: " +
                        examples.str() + ")");
    }
}

enum class BarcodeClass {
    Missing,
    CellExact,
    CellNormalized,
    EmptyExact,
    EmptyNormalized,
    Noncell
};

BarcodeClass classify_barcode(const std::string& barcode,
                              const BarcodeRoster& cells,
                              const BarcodeRoster* empties,
                              bool normalize) {
    if (cells.exact.count(barcode) != 0) return BarcodeClass::CellExact;
    if (empties != nullptr && empties->exact.count(barcode) != 0) {
        return BarcodeClass::EmptyExact;
    }
    if (normalize) {
        const std::string key = normalize_10x(barcode);
        if (cells.normalized.count(key) != 0) return BarcodeClass::CellNormalized;
        if (empties != nullptr && empties->normalized.count(key) != 0) {
            return BarcodeClass::EmptyNormalized;
        }
    }
    return BarcodeClass::Noncell;
}

struct Counts {
    std::array<uint64_t, 5> value{};

    bool any() const {
        return std::any_of(value.begin(), value.end(), [](uint64_t x) { return x != 0; });
    }

    void add(uint64_t bases, BarcodeClass barcode_class) {
        checked_add(value[0], bases, "counting bam_all aligned bases");
        if (barcode_class == BarcodeClass::Missing) return;
        checked_add(value[1], bases, "counting all_barcoded aligned bases");
        if (barcode_class == BarcodeClass::CellExact ||
            barcode_class == BarcodeClass::CellNormalized) {
            checked_add(value[2], bases, "counting cell aligned bases");
        } else if (barcode_class == BarcodeClass::EmptyExact ||
                   barcode_class == BarcodeClass::EmptyNormalized) {
            checked_add(value[4], bases, "counting empty aligned bases");
        } else {
            checked_add(value[3], bases, "counting noncell aligned bases");
        }
    }
};

struct Contig {
    std::string name;
    uint64_t length = 0;
};

struct BuildOptions {
    std::string bam;
    std::string cell_barcodes;
    std::string empty_barcodes;
    std::string output;
    std::string metadata;
    uint64_t window_size = kDefaultWindowSize;
    int min_mapq = 0;
    uint32_t exclude_flags = kDefaultExcludeFlags;
    std::string tag = "CB";
    int threads = kDefaultThreads;
    bool normalize = true;
    bool force = false;
};

struct MergeOptions {
    std::string input_list;
    std::string output_prefix;
    std::string metadata;
    std::vector<std::string> populations;
    int threads = kDefaultThreads;
    bool force = false;
};

std::string take_option_value(int& index, int argc, char** argv,
                              const std::string& option) {
    if (index + 1 >= argc) throw UserError(option + " requires a value");
    ++index;
    return argv[index];
}

void mark_once(std::set<std::string>& seen, const std::string& option) {
    if (!seen.insert(option).second) throw UserError("option specified more than once: " + option);
}

bool valid_tag_name(const std::string& tag) {
    if (tag.size() != 2) return false;
    const unsigned char first = static_cast<unsigned char>(tag[0]);
    const unsigned char second = static_cast<unsigned char>(tag[1]);
    return std::isalpha(first) && std::isalnum(second);
}

BuildOptions parse_build_options(int argc, char** argv) {
    BuildOptions options;
    std::set<std::string> seen;
    bool normalization_choice = false;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout << build_help();
            std::exit(0);
        } else if (arg == "--bam") {
            mark_once(seen, arg); options.bam = take_option_value(i, argc, argv, arg);
        } else if (arg == "--cell-barcodes") {
            mark_once(seen, arg); options.cell_barcodes = take_option_value(i, argc, argv, arg);
        } else if (arg == "--empty-barcodes") {
            mark_once(seen, arg); options.empty_barcodes = take_option_value(i, argc, argv, arg);
        } else if (arg == "--output") {
            mark_once(seen, arg); options.output = take_option_value(i, argc, argv, arg);
        } else if (arg == "--metadata") {
            mark_once(seen, arg); options.metadata = take_option_value(i, argc, argv, arg);
        } else if (arg == "--window-size") {
            mark_once(seen, arg);
            options.window_size = parse_u64(take_option_value(i, argc, argv, arg), arg);
        } else if (arg == "--min-mapq") {
            mark_once(seen, arg);
            const uint64_t value = parse_u64(take_option_value(i, argc, argv, arg), arg);
            if (value > 255) throw UserError("--min-mapq must be between 0 and 255");
            options.min_mapq = static_cast<int>(value);
        } else if (arg == "--exclude-flags") {
            mark_once(seen, arg);
            const std::string text = take_option_value(i, argc, argv, arg);
            const bool hexadecimal = text.size() > 2 && text[0] == '0' &&
                (text[1] == 'x' || text[1] == 'X');
            const uint64_t value = parse_u64(text, arg, hexadecimal ? 16 : 10);
            if (value > 0xFFFFU) throw UserError("--exclude-flags must fit in 16 bits");
            options.exclude_flags = static_cast<uint32_t>(value);
        } else if (arg == "--tag") {
            mark_once(seen, arg); options.tag = take_option_value(i, argc, argv, arg);
        } else if (arg == "--threads") {
            mark_once(seen, arg); options.threads = parse_threads(take_option_value(i, argc, argv, arg));
        } else if (arg == "--no-normalize-10x" || arg == "--normalize-10x") {
            if (normalization_choice) throw UserError("choose only one 10x normalization option");
            normalization_choice = true;
            options.normalize = arg == "--normalize-10x";
        } else if (arg == "--force") {
            mark_once(seen, arg); options.force = true;
        } else {
            throw UserError("unknown build option: " + arg);
        }
    }

    if (options.bam.empty()) throw UserError("build requires --bam");
    if (options.cell_barcodes.empty()) throw UserError("build requires --cell-barcodes");
    if (options.output.empty()) throw UserError("build requires --output");
    if (options.metadata.empty()) throw UserError("build requires --metadata");
    if (options.window_size == 0 || options.window_size >
        static_cast<uint64_t>(std::numeric_limits<int32_t>::max())) {
        throw UserError("--window-size must be between 1 and 2147483647");
    }
    if (!valid_tag_name(options.tag)) {
        throw UserError("--tag must be a valid two-character SAM tag");
    }

    require_absolute("--bam", options.bam);
    require_absolute("--cell-barcodes", options.cell_barcodes);
    require_absolute("--output", options.output);
    require_absolute("--metadata", options.metadata);
    if (!options.empty_barcodes.empty()) {
        require_absolute("--empty-barcodes", options.empty_barcodes);
    }
    require_input_file("--bam", options.bam);
    require_input_file("--cell-barcodes", options.cell_barcodes);
    if (!options.empty_barcodes.empty()) {
        require_input_file("--empty-barcodes", options.empty_barcodes);
    }
    require_output_parent("--output", options.output);
    require_output_parent("--metadata", options.metadata);
    std::vector<std::pair<std::string, std::string>> distinct = {
        {"--bam", options.bam}, {"--cell-barcodes", options.cell_barcodes},
        {"--output", options.output}, {"--metadata", options.metadata}};
    if (!options.empty_barcodes.empty()) {
        distinct.push_back({"--empty-barcodes", options.empty_barcodes});
    }
    require_distinct_paths(distinct);
    require_replaceable(options.output, options.force);
    require_replaceable(options.metadata, options.force);
    return options;
}

MergeOptions parse_merge_options(int argc, char** argv) {
    MergeOptions options;
    std::set<std::string> seen;
    std::string population_text;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout << merge_help();
            std::exit(0);
        } else if (arg == "--input-list") {
            mark_once(seen, arg); options.input_list = take_option_value(i, argc, argv, arg);
        } else if (arg == "--output-prefix") {
            mark_once(seen, arg); options.output_prefix = take_option_value(i, argc, argv, arg);
        } else if (arg == "--metadata") {
            mark_once(seen, arg); options.metadata = take_option_value(i, argc, argv, arg);
        } else if (arg == "--populations") {
            mark_once(seen, arg); population_text = take_option_value(i, argc, argv, arg);
        } else if (arg == "--threads") {
            mark_once(seen, arg); options.threads = parse_threads(take_option_value(i, argc, argv, arg));
        } else if (arg == "--force") {
            mark_once(seen, arg); options.force = true;
        } else {
            throw UserError("unknown merge option: " + arg);
        }
    }
    if (options.input_list.empty()) throw UserError("merge requires --input-list");
    if (options.output_prefix.empty()) throw UserError("merge requires --output-prefix");
    if (options.metadata.empty()) throw UserError("merge requires --metadata");
    require_absolute("--input-list", options.input_list);
    require_absolute("--output-prefix", options.output_prefix);
    require_absolute("--metadata", options.metadata);
    require_input_file("--input-list", options.input_list);
    require_output_parent("--output-prefix", options.output_prefix);
    require_output_parent("--metadata", options.metadata);

    if (population_text.empty()) {
        for (const char* population : kPopulations) options.populations.emplace_back(population);
    } else {
        const std::vector<std::string> requested = split_commas(population_text);
        std::set<std::string> requested_set;
        for (const std::string& population : requested) {
            if (population.empty()) throw UserError("--populations contains an empty name");
            if (!requested_set.insert(population).second) {
                throw UserError("duplicate population in --populations: " + population);
            }
            if (std::find_if(kPopulations.begin(), kPopulations.end(),
                    [&](const char* known) { return population == known; }) == kPopulations.end()) {
                throw UserError("unknown population in --populations: " + population);
            }
        }
        for (const char* known : kPopulations) {
            if (requested_set.count(known) != 0) options.populations.emplace_back(known);
        }
    }
    return options;
}

struct BuildStats {
    uint64_t records_total = 0;
    uint64_t records_filtered_flags = 0;
    uint64_t records_filtered_mapq = 0;
    uint64_t records_passing_filters = 0;
    uint64_t records_no_reference = 0;
    uint64_t records_with_aligned_bases = 0;
    uint64_t records_without_aligned_bases = 0;
    uint64_t records_without_barcode = 0;
    uint64_t records_cell_exact = 0;
    uint64_t records_cell_normalized = 0;
    uint64_t records_empty_exact = 0;
    uint64_t records_empty_normalized = 0;
    uint64_t records_noncell = 0;
    uint64_t windows_written = 0;
    Counts aligned_bases;
};

class WindowAccumulator {
public:
    WindowAccumulator(BgzfWriter& writer, const std::vector<Contig>& contigs,
                      uint64_t window_size, BuildStats& stats)
        : writer_(writer), contigs_(contigs), window_size_(window_size), stats_(stats) {}

    void switch_contig(int32_t tid) {
        if (current_tid_ == tid) return;
        flush_all();
        current_tid_ = tid;
    }

    void flush_before(uint64_t reference_position) {
        while (!windows_.empty()) {
            const uint64_t index = windows_.begin()->first;
            const uint64_t start = index * window_size_;
            const uint64_t end = std::min(start + window_size_,
                                          contigs_.at(current_tid_).length);
            if (end > reference_position) break;
            write_first();
        }
    }

    void add_segment(uint64_t start, uint64_t end, BarcodeClass barcode_class) {
        while (start < end) {
            const uint64_t index = start / window_size_;
            const uint64_t boundary = std::min(end, (index + 1) * window_size_);
            const uint64_t bases = boundary - start;
            windows_[index].add(bases, barcode_class);
            stats_.aligned_bases.add(bases, barcode_class);
            start = boundary;
        }
    }

    void finish() { flush_all(); }

private:
    void write_first() {
        const auto iterator = windows_.begin();
        const uint64_t index = iterator->first;
        const Counts& counts = iterator->second;
        if (counts.any()) {
            const uint64_t start = index * window_size_;
            const uint64_t end = std::min(start + window_size_,
                                          contigs_.at(current_tid_).length);
            std::string line = contigs_.at(current_tid_).name + "\t" +
                std::to_string(start) + "\t" + std::to_string(end);
            for (const uint64_t value : counts.value) line += "\t" + std::to_string(value);
            line.push_back('\n');
            writer_.write(line);
            ++stats_.windows_written;
        }
        windows_.erase(iterator);
    }

    void flush_all() {
        while (!windows_.empty()) write_first();
    }

    BgzfWriter& writer_;
    const std::vector<Contig>& contigs_;
    uint64_t window_size_;
    BuildStats& stats_;
    int32_t current_tid_ = -1;
    std::map<uint64_t, Counts> windows_;
};

std::string bam_record_label(const bam1_t* record, uint64_t ordinal) {
    const char* qname = bam_get_qname(record);
    return "record " + std::to_string(ordinal) +
           (qname == nullptr ? std::string() : " ('" + std::string(qname) + "')");
}

std::vector<Contig> extract_contigs(const sam_hdr_t* header) {
    const int count = sam_hdr_nref(header);
    if (count <= 0) throw UserError("BAM header contains no reference contigs");
    std::vector<Contig> contigs;
    contigs.reserve(static_cast<size_t>(count));
    std::unordered_set<std::string> names;
    for (int tid = 0; tid < count; ++tid) {
        const char* name_pointer = sam_hdr_tid2name(header, tid);
        const hts_pos_t length_value = sam_hdr_tid2len(header, tid);
        if (name_pointer == nullptr || *name_pointer == '\0') {
            throw UserError("BAM header has an unnamed contig at index " + std::to_string(tid));
        }
        const std::string name(name_pointer);
        if (has_forbidden_tsv_char(name)) {
            throw UserError("BAM contig name cannot be represented in TSV: " + name);
        }
        if (!names.insert(name).second) throw UserError("duplicate BAM contig name: " + name);
        if (length_value <= 0) {
            throw UserError("BAM contig has non-positive length: " + name);
        }
        contigs.push_back({name, static_cast<uint64_t>(length_value)});
    }
    return contigs;
}

std::string coverage_header(const std::vector<Contig>& contigs,
                            const BuildOptions& options) {
    std::string header;
    header += "##schema\t" + std::string(kSchema) + "\n";
    header += "##tool_version\t" + version_string() + "\n";
    header += "##window_size\t" + std::to_string(options.window_size) + "\n";
    header += "##coordinate_system\t0-based-half-open\n";
    header += "##count_unit\taligned_reference_bases\n";
    header += "##population_columns\tbam_all,all_barcoded,cells,noncell,empty\n";
    header += "##population_partition\tall_barcoded=cells+noncell+empty\n";
    header += "##min_mapq\t" + std::to_string(options.min_mapq) + "\n";
    header += "##exclude_flags\t" + std::to_string(options.exclude_flags) + "\n";
    header += "##barcode_tag\t" + options.tag + "\n";
    header += "##normalize_10x\t" +
              std::string(options.normalize ? "true" : "false") + "\n";
    header += "##empty_roster\t" +
              std::string(options.empty_barcodes.empty() ? "absent" : "present") + "\n";
    for (size_t i = 0; i < contigs.size(); ++i) {
        header += "##contig\t" + std::to_string(i) + "\t" + contigs[i].name +
                  "\t" + std::to_string(contigs[i].length) + "\n";
    }
    header += "#contig\tstart\tend\tbam_all\tall_barcoded\tcells\tnoncell\tempty\n";
    return header;
}

void increment_class_stat(BuildStats& stats, BarcodeClass barcode_class) {
    switch (barcode_class) {
        case BarcodeClass::Missing: ++stats.records_without_barcode; break;
        case BarcodeClass::CellExact: ++stats.records_cell_exact; break;
        case BarcodeClass::CellNormalized: ++stats.records_cell_normalized; break;
        case BarcodeClass::EmptyExact: ++stats.records_empty_exact; break;
        case BarcodeClass::EmptyNormalized: ++stats.records_empty_normalized; break;
        case BarcodeClass::Noncell: ++stats.records_noncell; break;
    }
}

std::vector<std::pair<std::string, std::string>> build_metadata_values(
        const BuildOptions& options, const BarcodeRoster& cells,
        const BarcodeRoster* empties, const std::vector<Contig>& contigs,
        const BuildStats& stats) {
    const uint64_t records_barcoded = stats.records_cell_exact +
        stats.records_cell_normalized + stats.records_empty_exact +
        stats.records_empty_normalized + stats.records_noncell;
    std::vector<std::pair<std::string, std::string>> values = {
        {"schema", kSchema}, {"tool_version", version_string()}, {"command", "build"},
        {"status", "complete"}, {"bam", options.bam},
        {"cell_barcodes", options.cell_barcodes},
        {"empty_barcodes", options.empty_barcodes.empty() ? "NONE" : options.empty_barcodes},
        {"coverage_output", options.output}, {"window_size", std::to_string(options.window_size)},
        {"min_mapq", std::to_string(options.min_mapq)},
        {"exclude_flags", "0x" + [&]() {
            std::ostringstream stream; stream << std::uppercase << std::hex << options.exclude_flags;
            return stream.str(); }()},
        {"barcode_tag", options.tag}, {"normalize_10x", options.normalize ? "true" : "false"},
        {"bam_reader_threads", std::to_string(options.threads)},
        {"contigs", std::to_string(contigs.size())},
        {"cell_roster_rows", std::to_string(cells.data_lines)},
        {"cell_roster_unique_exact", std::to_string(cells.exact.size())},
        {"cell_roster_unique_normalized", std::to_string(cells.normalized.size())},
        {"cell_roster_duplicate_exact", std::to_string(cells.duplicate_exact)},
        {"cell_roster_normalization_collisions", std::to_string(cells.normalization_collisions)},
        {"empty_roster_rows", std::to_string(empties == nullptr ? 0 : empties->data_lines)},
        {"empty_roster_unique_exact", std::to_string(empties == nullptr ? 0 : empties->exact.size())},
        {"empty_roster_unique_normalized", std::to_string(empties == nullptr ? 0 : empties->normalized.size())},
        {"empty_roster_duplicate_exact", std::to_string(empties == nullptr ? 0 : empties->duplicate_exact)},
        {"empty_roster_normalization_collisions", std::to_string(empties == nullptr ? 0 : empties->normalization_collisions)},
        {"records_total", std::to_string(stats.records_total)},
        {"records_filtered_flags", std::to_string(stats.records_filtered_flags)},
        {"records_filtered_mapq", std::to_string(stats.records_filtered_mapq)},
        {"records_passing_filters", std::to_string(stats.records_passing_filters)},
        {"records_no_reference", std::to_string(stats.records_no_reference)},
        {"records_with_aligned_bases", std::to_string(stats.records_with_aligned_bases)},
        {"records_without_aligned_bases", std::to_string(stats.records_without_aligned_bases)},
        {"records_without_barcode", std::to_string(stats.records_without_barcode)},
        {"records_barcoded", std::to_string(records_barcoded)},
        {"records_cell_exact", std::to_string(stats.records_cell_exact)},
        {"records_cell_normalized", std::to_string(stats.records_cell_normalized)},
        {"records_empty_exact", std::to_string(stats.records_empty_exact)},
        {"records_empty_normalized", std::to_string(stats.records_empty_normalized)},
        {"records_noncell", std::to_string(stats.records_noncell)},
        {"windows_written", std::to_string(stats.windows_written)}
    };
    for (size_t i = 0; i < kPopulations.size(); ++i) {
        values.push_back({"aligned_bases_" + std::string(kPopulations[i]),
                          std::to_string(stats.aligned_bases.value[i])});
    }
    return values;
}

int run_build(const BuildOptions& options) {
    BarcodeRoster cells = load_roster(options.cell_barcodes, options.normalize);
    std::unique_ptr<BarcodeRoster> empties;
    if (!options.empty_barcodes.empty()) {
        empties = std::make_unique<BarcodeRoster>(
            load_roster(options.empty_barcodes, options.normalize));
        ensure_rosters_disjoint(cells, *empties, options.normalize);
    }

    samFile* raw_input = sam_open(options.bam.c_str(), "r");
    if (raw_input == nullptr) throw UserError("cannot open BAM input: " + options.bam);
    struct SamCloser {
        void operator()(samFile* input) const { if (input != nullptr) sam_close(input); }
    };
    std::unique_ptr<samFile, SamCloser> input(raw_input);
    if (options.threads > 1 && hts_set_threads(input.get(), options.threads) != 0) {
        throw UserError("cannot enable HTSlib reader threads for " + options.bam);
    }
    sam_hdr_t* raw_header = sam_hdr_read(input.get());
    if (raw_header == nullptr) throw UserError("cannot read BAM header: " + options.bam);
    struct HeaderCloser {
        void operator()(sam_hdr_t* header) const { if (header != nullptr) sam_hdr_destroy(header); }
    };
    std::unique_ptr<sam_hdr_t, HeaderCloser> header(raw_header);
    const std::vector<Contig> contigs = extract_contigs(header.get());

    const std::string output_temporary = temporary_path_for(options.output);
    const std::string metadata_temporary = temporary_path_for(options.metadata);
    TempFiles temporary_files;
    temporary_files.add(output_temporary);
    temporary_files.add(metadata_temporary);

    BgzfWriter output(output_temporary);
    output.write(coverage_header(contigs, options));

    BuildStats stats;
    WindowAccumulator accumulator(output, contigs, options.window_size, stats);
    bam1_t* raw_record = bam_init1();
    if (raw_record == nullptr) throw UserError("cannot allocate BAM record");
    struct RecordCloser {
        void operator()(bam1_t* record) const { if (record != nullptr) bam_destroy1(record); }
    };
    std::unique_ptr<bam1_t, RecordCloser> record(raw_record);

    int32_t last_tid = -1;
    hts_pos_t last_position = -1;
    bool saw_unplaced = false;
    int read_status = 0;
    while ((read_status = sam_read1(input.get(), header.get(), record.get())) >= 0) {
        ++stats.records_total;
        const int32_t tid = record->core.tid;
        const hts_pos_t position = record->core.pos;
        if (tid < 0) {
            saw_unplaced = true;
        } else {
            if (tid >= static_cast<int32_t>(contigs.size()) || position < 0) {
                throw UserError("invalid reference coordinate at " +
                                bam_record_label(record.get(), stats.records_total));
            }
            if (saw_unplaced || tid < last_tid ||
                (tid == last_tid && position < last_position)) {
                throw UserError("BAM is not coordinate sorted at " +
                                bam_record_label(record.get(), stats.records_total));
            }
            last_tid = tid;
            last_position = position;
        }

        if ((static_cast<uint32_t>(record->core.flag) & options.exclude_flags) != 0) {
            ++stats.records_filtered_flags;
            continue;
        }
        if (record->core.qual < options.min_mapq) {
            ++stats.records_filtered_mapq;
            continue;
        }
        ++stats.records_passing_filters;
        if (tid < 0) {
            ++stats.records_no_reference;
            continue;
        }
        if (static_cast<uint64_t>(position) >= contigs[tid].length) {
            throw UserError("alignment starts beyond contig " + contigs[tid].name +
                            " at " + bam_record_label(record.get(), stats.records_total));
        }

        BarcodeClass barcode_class = BarcodeClass::Missing;
        uint8_t* auxiliary = bam_aux_get(record.get(), options.tag.c_str());
        if (auxiliary != nullptr) {
            if (*auxiliary != 'Z') {
                throw UserError("barcode tag " + options.tag + " is type '" +
                                std::string(1, static_cast<char>(*auxiliary)) +
                                "', expected 'Z', at " +
                                bam_record_label(record.get(), stats.records_total));
            }
            const char* barcode_pointer = bam_aux2Z(auxiliary);
            if (barcode_pointer == nullptr || *barcode_pointer == '\0') {
                throw UserError("barcode tag " + options.tag + " is empty or malformed at " +
                                bam_record_label(record.get(), stats.records_total));
            }
            barcode_class = classify_barcode(barcode_pointer, cells, empties.get(),
                                             options.normalize);
        }
        increment_class_stat(stats, barcode_class);

        accumulator.switch_contig(tid);
        accumulator.flush_before(static_cast<uint64_t>(position));
        uint64_t reference_position = static_cast<uint64_t>(position);
        uint64_t aligned_bases_this_record = 0;
        const uint32_t* cigar = bam_get_cigar(record.get());
        for (uint32_t cigar_index = 0; cigar_index < record->core.n_cigar; ++cigar_index) {
            const uint32_t operation = bam_cigar_op(cigar[cigar_index]);
            const uint64_t length = bam_cigar_oplen(cigar[cigar_index]);
            switch (operation) {
                case BAM_CMATCH:
                case BAM_CEQUAL:
                case BAM_CDIFF: {
                    if (length > contigs[tid].length - reference_position) {
                        throw UserError("CIGAR extends beyond contig " + contigs[tid].name +
                                        " at " + bam_record_label(record.get(), stats.records_total));
                    }
                    accumulator.add_segment(reference_position,
                                            reference_position + length, barcode_class);
                    checked_add(aligned_bases_this_record, length,
                                "counting bases in one BAM record");
                    reference_position += length;
                    break;
                }
                case BAM_CDEL:
                case BAM_CREF_SKIP:
                    if (length > contigs[tid].length - reference_position) {
                        throw UserError("CIGAR extends beyond contig " + contigs[tid].name +
                                        " at " + bam_record_label(record.get(), stats.records_total));
                    }
                    reference_position += length;
                    break;
                case BAM_CINS:
                case BAM_CSOFT_CLIP:
                case BAM_CHARD_CLIP:
                case BAM_CPAD:
                    break;
                default:
                    throw UserError("unsupported or malformed CIGAR operation at " +
                                    bam_record_label(record.get(), stats.records_total));
            }
        }
        if (aligned_bases_this_record == 0) ++stats.records_without_aligned_bases;
        else ++stats.records_with_aligned_bases;
    }
    if (read_status < -1) {
        throw UserError("HTSlib reported a truncated or malformed BAM while reading " +
                        options.bam);
    }
    if (stats.aligned_bases.value[2] == 0) {
        throw UserError("zero accepted aligned bases matched --cell-barcodes; verify the "
                        "BAM barcode tag, barcode namespace, normalization, and roster");
    }
    if (empties != nullptr && stats.aligned_bases.value[4] == 0) {
        throw UserError("zero accepted aligned bases matched --empty-barcodes; verify the "
                        "BAM barcode tag, barcode namespace, normalization, and roster, "
                        "or omit --empty-barcodes for a cell-only BAM");
    }
    accumulator.finish();
    output.close();

    write_metadata(metadata_temporary,
                   build_metadata_values(options, cells, empties.get(), contigs, stats));
    invalidate_old_metadata(options.metadata, options.force);
    atomic_publish(output_temporary, options.output);
    atomic_publish(metadata_temporary, options.metadata);
    return 0;
}

struct InputSpec {
    std::string library_id;
    std::string path;
};

std::vector<InputSpec> read_input_list(const std::string& path) {
    std::vector<InputSpec> inputs;
    std::unordered_set<std::string> library_ids;
    std::unordered_set<std::string> paths;
    bool first_data_line = true;
    read_gzip_transparent_lines(path, [&](const std::string& raw, uint64_t line_number) {
        const std::string stripped = trim(raw);
        if (stripped.empty() || stripped[0] == '#') return;
        std::vector<std::string> fields = split_tabs(raw);
        for (std::string& field : fields) field = trim(field);
        if (first_data_line && fields.size() == 2 && fields[0] == "library_id" &&
            fields[1] == "coverage_path") {
            first_data_line = false;
            return;
        }
        first_data_line = false;
        if (fields.size() != 2 || fields[0].empty() || fields[1].empty()) {
            throw UserError("input list requires exactly two non-empty TSV fields at " +
                            path + ":" + std::to_string(line_number));
        }
        if (has_forbidden_tsv_char(fields[0])) {
            throw UserError("invalid library_id at " + path + ":" +
                            std::to_string(line_number));
        }
        if (!library_ids.insert(fields[0]).second) {
            throw UserError("duplicate library_id in input list: " + fields[0]);
        }
        std::error_code path_error;
        require_absolute("coverage_path", fields[1]);
        require_input_file("coverage_path", fields[1]);
        const std::string normalized_path =
            fs::absolute(fs::path(fields[1]), path_error).lexically_normal().string();
        if (path_error) {
            throw UserError("cannot resolve coverage path at " + path + ":" +
                            std::to_string(line_number) + ": " + path_error.message());
        }
        if (!paths.insert(normalized_path).second) {
            throw UserError("duplicate coverage path in input list: " + fields[1]);
        }
        for (const InputSpec& previous : inputs) {
            std::error_code equivalent_error;
            const bool equivalent = fs::equivalent(previous.path, fields[1], equivalent_error);
            if (equivalent_error) {
                throw UserError("cannot compare coverage input paths: " +
                                equivalent_error.message());
            }
            if (equivalent) {
                throw UserError("coverage paths refer to the same file: " + previous.path +
                                " and " + fields[1]);
            }
        }
        inputs.push_back({fields[0], fields[1]});
    });
    if (inputs.empty()) throw UserError("input list contains no coverage rows: " + path);
    return inputs;
}

struct CoverageRow {
    int32_t tid = -1;
    uint64_t start = 0;
    uint64_t end = 0;
    Counts counts;
};

class CoverageReader {
public:
    CoverageReader(std::string library_id, std::string path, hts_tpool* pool)
        : library_id_(std::move(library_id)), path_(std::move(path)) {
        fp_ = bgzf_open(path_.c_str(), "r");
        if (fp_ == nullptr) throw UserError("cannot open BGZF coverage input: " + path_);
        try {
            attach_pool(pool);
            parse_header();
            advance();
        } catch (...) {
            cleanup();
            throw;
        }
    }

    CoverageReader(const CoverageReader&) = delete;
    CoverageReader& operator=(const CoverageReader&) = delete;

    ~CoverageReader() {
        cleanup();
    }

    const std::string& library_id() const { return library_id_; }
    const std::string& path() const { return path_; }
    const std::string& schema() const { return schema_; }
    uint64_t window_size() const { return window_size_; }
    int min_mapq() const { return min_mapq_; }
    uint32_t exclude_flags() const { return exclude_flags_; }
    const std::string& barcode_tag() const { return barcode_tag_; }
    bool normalize_10x_enabled() const { return normalize_10x_; }
    bool has_empty_roster() const { return has_empty_roster_; }
    const std::vector<Contig>& contigs() const { return contigs_; }
    bool has_row() const { return has_row_; }
    const CoverageRow& row() const { return row_; }
    uint64_t rows_read() const { return rows_read_; }

    void attach_pool(hts_tpool* pool) {
        if (pool != nullptr && bgzf_thread_pool(fp_, pool, 256) != 0) {
            throw UserError("cannot attach BGZF thread pool for " + path_);
        }
    }

    void advance() {
        const ssize_t length = bgzf_getline(fp_, '\n', &line_);
        if (length == -1) {
            has_row_ = false;
            return;
        }
        if (length < -1) throw UserError("BGZF read failure in " + path_);
        size_t text_length = static_cast<size_t>(length);
        if (text_length > 0 && line_.s[text_length - 1] == '\r') --text_length;
        const std::string_view text(line_.s, text_length);
        if (text.empty() || text.front() == '#') {
            throw UserError("unexpected blank/comment line after data header in " + path_);
        }
        std::array<std::string_view, 8> fields;
        if (!split_eight_tabs(text, fields)) {
            throw UserError("coverage row does not have 8 columns in " + path_);
        }
        int32_t tid = -1;
        if (fields[0] == previous_contig_name_) {
            tid = previous_contig_name_tid_;
        } else {
            const std::string contig_name(fields[0]);
            const auto contig_iterator = contig_to_tid_.find(contig_name);
            if (contig_iterator == contig_to_tid_.end()) {
                throw UserError("coverage row uses unknown contig '" + contig_name +
                                "' in " + path_);
            }
            previous_contig_name_ = contig_name;
            previous_contig_name_tid_ = contig_iterator->second;
            tid = contig_iterator->second;
        }
        CoverageRow parsed;
        parsed.tid = tid;
        parsed.start = parse_u64_view(fields[1], "coverage start", path_);
        parsed.end = parse_u64_view(fields[2], "coverage end", path_);
        for (size_t i = 0; i < parsed.counts.value.size(); ++i) {
            parsed.counts.value[i] = parse_u64_view(
                fields[i + 3], "coverage count", path_);
        }
        const Contig& contig = contigs_.at(parsed.tid);
        if (parsed.start >= contig.length || parsed.start % window_size_ != 0) {
            throw UserError("invalid window start in " + path_ + ": " +
                            std::string(text));
        }
        const uint64_t expected_end = std::min(parsed.start + window_size_, contig.length);
        if (parsed.end != expected_end || parsed.end <= parsed.start) {
            throw UserError("invalid window end in " + path_ + ": " +
                            std::string(text));
        }
        if (!parsed.counts.any()) {
            throw UserError("non-sparse all-zero coverage row in " + path_ + ": " +
                            std::string(text));
        }
        uint64_t classified_barcoded = 0;
        checked_add(classified_barcoded, parsed.counts.value[2],
                    "validating cell coverage counts");
        checked_add(classified_barcoded, parsed.counts.value[3],
                    "validating noncell coverage counts");
        checked_add(classified_barcoded, parsed.counts.value[4],
                    "validating empty coverage counts");
        if (parsed.counts.value[1] != classified_barcoded) {
            throw UserError("population counts are inconsistent in " + path_ + ": " +
                            std::string(text));
        }
        if (parsed.counts.value[1] > parsed.counts.value[0]) {
            throw UserError("all_barcoded exceeds bam_all in " + path_ + ": " +
                            std::string(text));
        }
        if (!has_empty_roster_ && parsed.counts.value[4] != 0) {
            throw UserError("empty coverage is nonzero despite empty_roster=absent in " +
                            path_ + ": " + std::string(text));
        }
        if (have_previous_ && (parsed.tid < previous_tid_ ||
            (parsed.tid == previous_tid_ && parsed.start <= previous_start_))) {
            throw UserError("coverage rows are not strictly sorted in BAM contig order: " + path_);
        }
        previous_tid_ = parsed.tid;
        previous_start_ = parsed.start;
        have_previous_ = true;
        row_ = parsed;
        has_row_ = true;
        ++rows_read_;
    }

private:
    void cleanup() noexcept {
        if (line_.s != nullptr) {
            std::free(line_.s);
            line_ = {0, 0, nullptr};
        }
        if (fp_ != nullptr) {
            bgzf_close(fp_);
            fp_ = nullptr;
        }
    }

    std::string read_header_line() {
        const ssize_t length = bgzf_getline(fp_, '\n', &line_);
        if (length == -1) throw UserError("truncated coverage header in " + path_);
        if (length < -1) throw UserError("BGZF read failure in " + path_);
        std::string text(line_.s, static_cast<size_t>(length));
        if (!text.empty() && text.back() == '\r') text.pop_back();
        return text;
    }

    void parse_header() {
        bool saw_schema = false;
        bool saw_window = false;
        bool saw_coordinate = false;
        bool saw_unit = false;
        bool saw_populations = false;
        bool saw_partition = false;
        bool saw_min_mapq = false;
        bool saw_exclude_flags = false;
        bool saw_barcode_tag = false;
        bool saw_normalize = false;
        bool saw_empty_roster = false;
        while (true) {
            const std::string line = read_header_line();
            if (line == "#contig\tstart\tend\tbam_all\tall_barcoded\tcells\tnoncell\tempty") {
                break;
            }
            const std::vector<std::string> fields = split_tabs(line);
            if (fields.empty() || fields[0].rfind("##", 0) != 0) {
                throw UserError("malformed coverage header in " + path_ + ": " + line);
            }
            if (fields[0] == "##schema") {
                if (saw_schema || fields.size() != 2) throw UserError("malformed schema header in " + path_);
                schema_ = fields[1]; saw_schema = true;
            } else if (fields[0] == "##window_size") {
                if (saw_window || fields.size() != 2) throw UserError("malformed window header in " + path_);
                window_size_ = parse_u64(fields[1], "coverage window_size");
                if (window_size_ == 0) throw UserError("zero coverage window_size in " + path_);
                saw_window = true;
            } else if (fields[0] == "##coordinate_system") {
                if (saw_coordinate || fields.size() != 2 || fields[1] != "0-based-half-open") {
                    throw UserError("unsupported coordinate header in " + path_);
                }
                saw_coordinate = true;
            } else if (fields[0] == "##count_unit") {
                if (saw_unit || fields.size() != 2 || fields[1] != "aligned_reference_bases") {
                    throw UserError("unsupported count unit in " + path_);
                }
                saw_unit = true;
            } else if (fields[0] == "##population_columns") {
                if (saw_populations || fields.size() != 2 ||
                    fields[1] != "bam_all,all_barcoded,cells,noncell,empty") {
                    throw UserError("unsupported population columns in " + path_);
                }
                saw_populations = true;
            } else if (fields[0] == "##population_partition") {
                if (saw_partition || fields.size() != 2 ||
                    fields[1] != "all_barcoded=cells+noncell+empty") {
                    throw UserError("unsupported population partition in " + path_);
                }
                saw_partition = true;
            } else if (fields[0] == "##min_mapq") {
                if (saw_min_mapq || fields.size() != 2) {
                    throw UserError("malformed min_mapq header in " + path_);
                }
                const uint64_t value = parse_u64(fields[1], "coverage min_mapq");
                if (value > 255) throw UserError("invalid min_mapq header in " + path_);
                min_mapq_ = static_cast<int>(value);
                saw_min_mapq = true;
            } else if (fields[0] == "##exclude_flags") {
                if (saw_exclude_flags || fields.size() != 2) {
                    throw UserError("malformed exclude_flags header in " + path_);
                }
                const uint64_t value = parse_u64(fields[1], "coverage exclude_flags");
                if (value > 0xFFFFU) throw UserError("invalid exclude_flags header in " + path_);
                exclude_flags_ = static_cast<uint32_t>(value);
                saw_exclude_flags = true;
            } else if (fields[0] == "##barcode_tag") {
                if (saw_barcode_tag || fields.size() != 2 || !valid_tag_name(fields[1])) {
                    throw UserError("invalid barcode_tag header in " + path_);
                }
                barcode_tag_ = fields[1];
                saw_barcode_tag = true;
            } else if (fields[0] == "##normalize_10x") {
                if (saw_normalize || fields.size() != 2 ||
                    (fields[1] != "true" && fields[1] != "false")) {
                    throw UserError("invalid normalize_10x header in " + path_);
                }
                normalize_10x_ = fields[1] == "true";
                saw_normalize = true;
            } else if (fields[0] == "##empty_roster") {
                if (saw_empty_roster || fields.size() != 2 ||
                    (fields[1] != "present" && fields[1] != "absent")) {
                    throw UserError("invalid empty_roster header in " + path_);
                }
                has_empty_roster_ = fields[1] == "present";
                saw_empty_roster = true;
            } else if (fields[0] == "##contig") {
                if (fields.size() != 4) throw UserError("malformed contig header in " + path_);
                const uint64_t index = parse_u64(fields[1], "coverage contig index");
                const uint64_t length = parse_u64(fields[3], "coverage contig length");
                if (index != contigs_.size() || fields[2].empty() || length == 0 ||
                    contig_to_tid_.count(fields[2]) != 0) {
                    throw UserError("invalid contig dictionary in " + path_);
                }
                contig_to_tid_[fields[2]] = static_cast<int32_t>(index);
                contigs_.push_back({fields[2], length});
            } else if (fields[0] == "##tool_version") {
                if (fields.size() != 2) throw UserError("malformed tool version header in " + path_);
            } else {
                throw UserError("unknown coverage header field in " + path_ + ": " + fields[0]);
            }
        }
        if (!saw_schema || schema_ != kSchema || !saw_window || !saw_coordinate ||
            !saw_unit || !saw_populations || !saw_partition || !saw_min_mapq ||
            !saw_exclude_flags || !saw_barcode_tag || !saw_normalize ||
            !saw_empty_roster || contigs_.empty()) {
            throw UserError("incomplete or incompatible coverage header in " + path_);
        }
    }

    std::string library_id_;
    std::string path_;
    BGZF* fp_ = nullptr;
    kstring_t line_{0, 0, nullptr};
    std::string schema_;
    uint64_t window_size_ = 0;
    int min_mapq_ = 0;
    uint32_t exclude_flags_ = 0;
    std::string barcode_tag_;
    bool normalize_10x_ = false;
    bool has_empty_roster_ = false;
    std::vector<Contig> contigs_;
    std::unordered_map<std::string, int32_t> contig_to_tid_;
    std::string previous_contig_name_;
    int32_t previous_contig_name_tid_ = -1;
    CoverageRow row_;
    bool has_row_ = false;
    bool have_previous_ = false;
    int32_t previous_tid_ = -1;
    uint64_t previous_start_ = 0;
    uint64_t rows_read_ = 0;
};

bool same_dictionary(const std::vector<Contig>& left,
                     const std::vector<Contig>& right) {
    if (left.size() != right.size()) return false;
    for (size_t i = 0; i < left.size(); ++i) {
        if (left[i].name != right[i].name || left[i].length != right[i].length) return false;
    }
    return true;
}

size_t population_index(const std::string& population) {
    for (size_t i = 0; i < kPopulations.size(); ++i) {
        if (population == kPopulations[i]) return i;
    }
    throw UserError("internal error: unknown population " + population);
}

std::string join(const std::vector<std::string>& values, const std::string& separator) {
    std::string result;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) result += separator;
        result += values[i];
    }
    return result;
}

std::string format_coverage(uint64_t aligned_bases, uint64_t window_length) {
    if (aligned_bases % window_length == 0) {
        return std::to_string(aligned_bases / window_length);
    }
    std::ostringstream stream;
    stream << std::setprecision(12) << std::defaultfloat
           << static_cast<long double>(aligned_bases) /
              static_cast<long double>(window_length);
    return stream.str();
}

struct HeapEntry {
    int32_t tid;
    uint64_t start;
    size_t reader_index;
};

struct HeapLater {
    bool operator()(const HeapEntry& left, const HeapEntry& right) const {
        if (left.tid != right.tid) return left.tid > right.tid;
        if (left.start != right.start) return left.start > right.start;
        return left.reader_index > right.reader_index;
    }
};

int run_merge(const MergeOptions& options) {
    const std::vector<InputSpec> specs = read_input_list(options.input_list);
    hts_tpool* pool = nullptr;
    if (options.threads > 1) {
        pool = hts_tpool_init(options.threads);
        if (pool == nullptr) throw UserError("cannot create shared BGZF thread pool");
    }
    struct PoolCloser {
        void operator()(hts_tpool* value) const {
            if (value != nullptr) hts_tpool_destroy(value);
        }
    };
    std::unique_ptr<hts_tpool, PoolCloser> pool_guard(pool);

    std::vector<std::unique_ptr<CoverageReader>> readers;
    readers.reserve(specs.size());
    for (const InputSpec& spec : specs) {
        readers.push_back(std::make_unique<CoverageReader>(
            spec.library_id, spec.path, pool));
    }
    const uint64_t window_size = readers.front()->window_size();
    const std::vector<Contig>& contigs = readers.front()->contigs();
    for (size_t i = 1; i < readers.size(); ++i) {
        if (readers[i]->schema() != readers.front()->schema()) {
            throw UserError("coverage schema mismatch for library " + readers[i]->library_id());
        }
        if (readers[i]->window_size() != window_size) {
            throw UserError("window-size mismatch for library " + readers[i]->library_id());
        }
        if (!same_dictionary(readers[i]->contigs(), contigs)) {
            throw UserError("BAM contig dictionary mismatch for library " + readers[i]->library_id());
        }
        if (readers[i]->min_mapq() != readers.front()->min_mapq() ||
            readers[i]->exclude_flags() != readers.front()->exclude_flags() ||
            readers[i]->barcode_tag() != readers.front()->barcode_tag() ||
            readers[i]->normalize_10x_enabled() != readers.front()->normalize_10x_enabled() ||
            readers[i]->has_empty_roster() != readers.front()->has_empty_roster()) {
            throw UserError("coverage filter/barcode contract mismatch for library " +
                            readers[i]->library_id());
        }
    }

    std::vector<std::string> output_paths;
    std::vector<std::string> temporary_paths;
    for (const std::string& population : options.populations) {
        output_paths.push_back(options.output_prefix + "." + population + ".bedGraph.gz");
    }
    std::vector<std::pair<std::string, std::string>> distinct = {
        {"--input-list", options.input_list}, {"--metadata", options.metadata}};
    for (const InputSpec& spec : specs) {
        distinct.push_back({"coverage input " + spec.library_id, spec.path});
    }
    for (size_t i = 0; i < output_paths.size(); ++i) {
        require_output_parent("aggregate output", output_paths[i]);
        require_replaceable(output_paths[i], options.force);
        distinct.push_back({"aggregate output " + options.populations[i], output_paths[i]});
    }
    require_replaceable(options.metadata, options.force);
    require_distinct_paths(distinct);

    TempFiles temporary_files;
    for (const std::string& path : output_paths) {
        temporary_paths.push_back(temporary_path_for(path));
        temporary_files.add(temporary_paths.back());
    }
    const std::string metadata_temporary = temporary_path_for(options.metadata);
    temporary_files.add(metadata_temporary);

    std::vector<std::unique_ptr<BgzfWriter>> outputs;
    outputs.reserve(temporary_paths.size());
    for (const std::string& path : temporary_paths) {
        outputs.push_back(std::make_unique<BgzfWriter>(path));
        outputs.back()->attach_pool(pool);
    }

    std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapLater> heap;
    for (size_t i = 0; i < readers.size(); ++i) {
        if (readers[i]->has_row()) {
            heap.push({readers[i]->row().tid, readers[i]->row().start, i});
        }
    }

    uint64_t union_windows = 0;
    std::array<uint64_t, 5> aligned_totals{};
    std::array<uint64_t, 5> output_rows{};
    while (!heap.empty()) {
        const int32_t tid = heap.top().tid;
        const uint64_t start = heap.top().start;
        const uint64_t end = std::min(start + window_size, contigs.at(tid).length);
        Counts sums;
        while (!heap.empty() && heap.top().tid == tid && heap.top().start == start) {
            const size_t reader_index = heap.top().reader_index;
            heap.pop();
            const CoverageRow& row = readers[reader_index]->row();
            if (row.end != end) {
                throw UserError("window boundary mismatch for library " +
                                readers[reader_index]->library_id());
            }
            for (size_t population = 0; population < sums.value.size(); ++population) {
                checked_add(sums.value[population], row.counts.value[population],
                            "merging coverage counts");
            }
            readers[reader_index]->advance();
            if (readers[reader_index]->has_row()) {
                heap.push({readers[reader_index]->row().tid,
                           readers[reader_index]->row().start, reader_index});
            }
        }
        ++union_windows;
        const uint64_t actual_length = end - start;
        for (size_t population = 0; population < sums.value.size(); ++population) {
            checked_add(aligned_totals[population], sums.value[population],
                        "summing aggregate metadata counts");
        }
        for (size_t output_index = 0; output_index < options.populations.size(); ++output_index) {
            const size_t count_index = population_index(options.populations[output_index]);
            const uint64_t count = sums.value[count_index];
            if (count == 0) continue;
            outputs[output_index]->write(contigs.at(tid).name + "\t" +
                std::to_string(start) + "\t" + std::to_string(end) + "\t" +
                format_coverage(count, actual_length) + "\n");
            ++output_rows[count_index];
        }
    }
    for (auto& output : outputs) output->close();
    outputs.clear();

    uint64_t input_rows = 0;
    for (const auto& reader : readers) {
        checked_add(input_rows, reader->rows_read(), "summing input row counts");
    }
    std::vector<std::pair<std::string, std::string>> metadata = {
        {"schema", kSchema}, {"tool_version", version_string()}, {"command", "merge"},
        {"status", "complete"}, {"input_list", options.input_list},
        {"output_prefix", options.output_prefix}, {"window_size", std::to_string(window_size)},
        {"contigs", std::to_string(contigs.size())},
        {"libraries", std::to_string(readers.size())},
        {"populations", join(options.populations, ",")},
        {"min_mapq", std::to_string(readers.front()->min_mapq())},
        {"exclude_flags", std::to_string(readers.front()->exclude_flags())},
        {"barcode_tag", readers.front()->barcode_tag()},
        {"normalize_10x", readers.front()->normalize_10x_enabled() ? "true" : "false"},
        {"empty_roster", readers.front()->has_empty_roster() ? "present" : "absent"},
        {"bgzf_threads", std::to_string(options.threads)},
        {"input_sparse_rows", std::to_string(input_rows)},
        {"union_windows", std::to_string(union_windows)}
    };
    for (size_t i = 0; i < readers.size(); ++i) {
        metadata.push_back({"library_" + std::to_string(i + 1),
                            readers[i]->library_id() + "\t" + readers[i]->path()});
    }
    for (size_t i = 0; i < kPopulations.size(); ++i) {
        metadata.push_back({"aligned_bases_" + std::string(kPopulations[i]),
                            std::to_string(aligned_totals[i])});
        if (std::find(options.populations.begin(), options.populations.end(),
                      kPopulations[i]) != options.populations.end()) {
            metadata.push_back({"output_rows_" + std::string(kPopulations[i]),
                                std::to_string(output_rows[i])});
        }
    }
    write_metadata(metadata_temporary, metadata);
    invalidate_old_metadata(options.metadata, options.force);
    for (size_t i = 0; i < output_paths.size(); ++i) {
        atomic_publish(temporary_paths[i], output_paths[i]);
    }
    atomic_publish(metadata_temporary, options.metadata);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 1) {
            std::cout << main_help();
            return 0;
        }
        const std::string command = argv[1];
        if (command == "--version" || command == "-V") {
            if (argc != 2) throw UserError("--version does not accept arguments");
            std::cout << version_string() << '\n';
            return 0;
        }
        if (command == "--help" || command == "-h" || command == "help") {
            if (argc != 2) throw UserError("use '<subcommand> --help' for subcommand help");
            std::cout << main_help();
            return 0;
        }
        if (command == "build") return run_build(parse_build_options(argc, argv));
        if (command == "merge") return run_merge(parse_merge_options(argc, argv));
        throw UserError("unknown subcommand: " + command);
    } catch (const UserError& error) {
        std::cerr << "ERROR: " << error.what() << "\n\n";
        std::cerr << main_help();
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: unexpected failure: " << error.what() << '\n';
        return 1;
    }
}
