SHELL = bash
COMP = g++
CCOMP = gcc
PREFIX ?= /usr/local
SOURCE_REVISION ?= $(shell git rev-parse --short HEAD 2>/dev/null || echo archive)

# Cluster HTSlib installation. Override these variables for another system.
HTSLIB_PREFIX ?= /nvme/software/packages/htslib/1.20
HTSLIB_INCLUDE ?= $(HTSLIB_PREFIX)/include
HTSLIB_LIB ?= $(HTSLIB_PREFIX)/lib
HTSLIB_RPATH_FLAG = -Wl,-rpath,$(HTSLIB_LIB)

# One build is deployed across Zen3 and Zen4 cluster nodes.
CPU_ARCH ?= znver3
CPU_TUNE ?= znver3
ARCHFLAGS ?= -march=$(CPU_ARCH) -mtune=$(CPU_TUNE)

BC_LENX2 = 32
KX2 = 16
NBITS ?= 2048
MAX_SITES ?= 2000
MAKE = make
PROJROOT = $(shell pwd)

CXXFLAGS_STD = -std=c++11 -fPIC -D_REENTRANT -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2)
CXXFLAGS_PARALLEL = -std=c++11 -fPIC -D_REENTRANT -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -DNBITS=$(NBITS) \
                    -DCELLBOUNCER_SOURCE_REVISION=\"$(SOURCE_REVISION)\" \
                    -O3 $(ARCHFLAGS) -fopenmp
CXXFLAGS_TET = -std=c++11 -fPIC -D_REENTRANT -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -DNBITS=$(NBITS) \
               -O3 $(ARCHFLAGS) -fopenmp
CXXFLAGS_CACHE = -std=c++17 -fPIC -D_REENTRANT -O3 $(ARCHFLAGS) -Wall -Wextra -pedantic
CXXFLAGS_MT = -std=c++17 -fPIC -D_REENTRANT -O3 $(ARCHFLAGS) -Wall -Wextra -pedantic
CXXFLAGS_SCRUB = -std=c++11 -fPIC -D_REENTRANT -O3 $(ARCHFLAGS) -Wall -Wextra
CFLAGS = -fPIC -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -O3 $(ARCHFLAGS)

# Refactored sources include the unified headers directly.
CXXIFLAGS = -I$(HTSLIB_INCLUDE) -I$(PREFIX)/include -Iinclude -Isrc
CIFLAGS = -I$(HTSLIB_INCLUDE) -I$(PREFIX)/include -Iinclude
LFLAGS = -L$(HTSLIB_LIB) $(HTSLIB_RPATH_FLAG) -L$(PREFIX)/lib -Llib
LFLAGS_PARALLEL = -L$(HTSLIB_LIB) $(HTSLIB_RPATH_FLAG) -L$(PREFIX)/lib -Llib -fopenmp -flto=auto
LFLAGS_TET = -L$(HTSLIB_LIB) $(HTSLIB_RPATH_FLAG) -L$(PREFIX)/lib -Llib -fopenmp

ifeq ($(findstring cellbouncer,${CONDA_PREFIX}),cellbouncer)
    CXXIFLAGS += -I${CONDA_PREFIX}/include
    CIFLAGS += -I${CONDA_PREFIX}/include
    LFLAGS += -L${CONDA_PREFIX}/lib
    LFLAGS_PARALLEL += -L${CONDA_PREFIX}/lib
    LFLAGS_TET += -L${CONDA_PREFIX}/lib
endif

DEPS = lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a
DEPS2 = -lz -lhts -lpthread
DEPS2_PARALLEL = -lz -lhts -lpthread -lrt
DEPS_CACHE = -lhts -lz -lpthread
DEPS_MT = -lhts -lz -lpthread
DEPS_SCRUB = -lhts -lz -lpthread
BUILD_DIR_STAMP = build/.directory_ready
HTSWRAPPER_HEADER_STAMP = include/htswrapper/.headers_installed

# optimML is patched only in a private build copy. The vendored source remains
# unchanged, and the patch is rebuilt when this Makefile changes.
OPTIMML_PATCH_DIR = build/optimML_cellbouncer
OPTIMML_PATCH_STAMP = $(OPTIMML_PATCH_DIR)/.cellbouncer_patch_v3
OPTIMML_SOURCE_FILES = $(shell find dependencies/optimML -type f \
    -not -path '*/build/*' -not -path '*/lib/*' 2>/dev/null)

# Executable names consumed by orchestrate_pipeline.py are unchanged.
ORCHESTRATOR_BINS = demux_parallel vcf_loader_daemon tet_ambient_profile \
                    tet_contam_estimate legacy2c_contam_estimate \
                    tetra_score_calls tetra_refine
AUX_ROOT_BINS = demux_mt demux_species demux_tags doublet_dragon bulkprops \
                bam_ram_host_daemon genotype_scrub_bam snps_per_read \
                mt_fusion_ratio
DEPRECATED_BINS = demux_vcf quant_contam quant3_contam quant3_contam_ap \
                  quant3_contam_empty_drops downsample_vcf
UTIL_BINS = utils/refine_vcf utils/bam_indiv_rg utils/bam_split_bcs \
            utils/bam_cb_cache_extract utils/split_read_files \
            utils/atac_fq_preprocess utils/combine_species_counts \
            utils/composite_bam2counts utils/downsample_vcf_parallel

# Preserve the unrelated FASTK utility when the full upstream source subtree is
# present in the real checkout. The lean source bundle does not carry FASTK, so
# this target is selected only where it was already available.
ifneq ($(and $(wildcard src/FASTK/libfastk.c),$(wildcard src/FASTK/libfastk.h)),)
UTIL_BINS += utils/get_unique_kmers
endif

all: dependencies orchestrator_tools auxiliary_tools utilities

orchestrator_tools: $(ORCHESTRATOR_BINS)
auxiliary_tools: $(AUX_ROOT_BINS)
utilities: $(UTIL_BINS)
dependencies: $(DEPS)

source_inventory:
	@echo "Unified implementation modules:"; \
	for f in src/io.cpp src/vcf_hts.cpp src/genotype_llr.cpp src/ambient_rna_three_ap.cpp; do \
	    if [ -s "$$f" ]; then echo "  PRESENT $$f"; else echo "  MISSING $$f"; fi; \
	done; \
	echo; echo "Orchestrator executables:"; \
	for f in $(ORCHESTRATOR_BINS); do echo "  $$f"; done

print_flags:
	@echo "HTSLIB_PREFIX=$(HTSLIB_PREFIX)"
	@echo "CPU_ARCH=$(CPU_ARCH)"
	@echo "CPU_TUNE=$(CPU_TUNE)"
	@echo "ARCHFLAGS=$(ARCHFLAGS)"
	@echo "SOURCE_REVISION=$(SOURCE_REVISION)"
	@echo "CXXFLAGS_PARALLEL=$(CXXFLAGS_PARALLEL)"
	@echo "CXXFLAGS_TET=$(CXXFLAGS_TET)"

check_required_ambient_source:
	@if [ ! -s src/ambient_rna_three_ap.cpp ]; then \
	    echo "ERROR: src/ambient_rna_three_ap.cpp is missing." >&2; \
	    exit 2; \
	fi

# -----------------------------------------------------------------------------
# Orchestrator-facing executables
# -----------------------------------------------------------------------------

demux_parallel: src/demux_parallel.cpp src/io.h src/vcf_hts.h src/genotype_llr.h \
    build/common_parallel.o build/io_parallel.o build/vcf_hts.o build/genotype_llr.o \
    $(DEPS) $(HTSWRAPPER_HEADER_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g \
	    build/common_parallel.o build/io_parallel.o build/vcf_hts.o build/genotype_llr.o \
	    src/demux_parallel.cpp -o $@ $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

vcf_loader_daemon: src/vcf_loader_daemon.cpp src/vcf_hts.h \
    build/common_parallel.o build/vcf_hts.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g \
	    build/common_parallel.o build/vcf_hts.o src/vcf_loader_daemon.cpp \
	    -o $@ $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

build/ambient_rna_three_ap_tet.o: check_required_ambient_source \
    src/ambient_rna_three_ap.cpp src/ambient_rna_three_ap.h src/genotype_llr.h \
    src/common.h $(DEPS) | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -g \
	    src/ambient_rna_three_ap.cpp -c -o $@

tet_ambient_profile: src/tet_ambient_profile.cpp src/ambient_rna_three_ap.h src/io.h \
    build/common_parallel.o build/io_parallel.o build/vcf_hts.o build/genotype_llr.o \
    build/ambient_rna_three_ap_tet.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g \
	    build/common_parallel.o build/io_parallel.o build/vcf_hts.o build/genotype_llr.o \
	    build/ambient_rna_three_ap_tet.o src/tet_ambient_profile.cpp \
	    -o $@ $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

tet_contam_estimate: src/tet_contam_estimate.cpp src/ambient_rna_three_ap.h \
    src/ambient_rna_gex.h src/io.h build/common_parallel.o build/io_parallel.o \
    build/vcf_hts.o build/genotype_llr.o build/ambient_rna_three_ap_tet.o \
    build/ambient_rna_gex.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g \
	    build/common_parallel.o build/io_parallel.o build/vcf_hts.o build/genotype_llr.o \
	    build/ambient_rna_three_ap_tet.o build/ambient_rna_gex.o \
	    src/tet_contam_estimate.cpp -o $@ $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

build/legacy2c_model.o: src/legacy2c_model.cpp src/legacy2c_model.h \
    src/genotype_llr.h src/common.h $(DEPS) | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g src/legacy2c_model.cpp -c -o $@

legacy2c_contam_estimate: src/legacy2c_contam_estimate.cpp src/legacy2c_model.h src/io.h \
    build/common_parallel.o build/io_parallel.o build/vcf_hts.o build/genotype_llr.o \
    build/legacy2c_model.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g \
	    build/common_parallel.o build/io_parallel.o build/vcf_hts.o build/genotype_llr.o \
	    build/legacy2c_model.o src/legacy2c_contam_estimate.cpp \
	    -o $@ $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

tetra_refine: src/tetra_refine.cpp lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 src/tetra_refine.cpp \
	    -o $@ $(LFLAGS) lib/libhtswrapper.a -lz

tetra_score_calls: src/tetra_score_calls.cpp lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -g src/tetra_score_calls.cpp \
	    -o $@ $(LFLAGS_TET) lib/libhtswrapper.a -lz

# -----------------------------------------------------------------------------
# Build-target compatibility for the existing compcb shell function
# -----------------------------------------------------------------------------
# The retired quant3 source files and executables are intentionally gone.  These
# names remain as phony make targets only because the established compcb command
# invokes them explicitly after `make all`.  They build the production
# replacements without recreating or installing deprecated quant3 binaries.
quant3_contam: tet_contam_estimate
	@echo "quant3_contam compatibility target satisfied by tet_contam_estimate"

quant3_contam_ap: tet_contam_estimate
	@echo "quant3_contam_ap compatibility target satisfied by tet_contam_estimate"

quant3_contam_empty_drops: tet_ambient_profile
	@echo "quant3_contam_empty_drops compatibility target satisfied by tet_ambient_profile"

# -----------------------------------------------------------------------------
# Retained auxiliary tools, all linked through the unified modules where relevant
# -----------------------------------------------------------------------------

demux_mt: src/demux_mt.cpp src/genotype_llr.h build/common_parallel.o \
    build/vcf_hts.o build/genotype_llr.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -D MAX_SITES=$(MAX_SITES) \
	    build/common_parallel.o build/vcf_hts.o build/genotype_llr.o src/demux_mt.cpp \
	    -o $@ $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

demux_species: src/demux_species.cpp build/common.o build/demux_species_io.o \
    build/species_kmers.o build/reads_demux.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 -g build/common.o \
	    build/demux_species_io.o build/species_kmers.o build/reads_demux.o \
	    src/demux_species.cpp $(LFLAGS) $(DEPS) -pthread -o $@ $(DEPS2)

demux_tags: src/demux_tags.cpp build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/demux_tags.cpp \
	    $(LFLAGS) $(DEPS) -D PROJ_ROOT=$(PROJROOT) -o $@ $(DEPS2)

doublet_dragon: src/doublet_dragon.cpp build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/doublet_dragon.cpp \
	    $(LFLAGS) $(DEPS) -o $@ $(DEPS2)

bulkprops: src/bulkprops.cpp src/vcf_hts.h src/io.h build/common_parallel.o \
    build/io_parallel.o build/vcf_hts.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) build/common_parallel.o \
	    build/io_parallel.o build/vcf_hts.o src/bulkprops.cpp \
	    -o $@ $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

bam_ram_host_daemon: src/bam_ram_host_daemon.cpp
	$(COMP) -std=c++11 -O3 $(ARCHFLAGS) src/bam_ram_host_daemon.cpp -o $@ -lpthread

build/genotype_scrub_bam.o: src/genotype_scrub_bam.cpp Makefile | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_SCRUB) -c src/genotype_scrub_bam.cpp -o $@

genotype_scrub_bam: build/genotype_scrub_bam.o
	$(COMP) $< -o $@ $(LFLAGS) $(DEPS_SCRUB)

snps_per_read: src/snps_per_read.cpp | lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) -std=c++17 -fPIC -D_REENTRANT -O3 $(ARCHFLAGS) \
	    -fopenmp -Wall -Wextra src/snps_per_read.cpp -o $@ \
	    $(LFLAGS) -lhts -lz -lpthread

mt_fusion_ratio: src/mt_fusion_ratio.cpp
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_MT) src/mt_fusion_ratio.cpp \
	    -o $@ $(LFLAGS) $(DEPS_MT)

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

utils/refine_vcf: src/refine_vcf.cpp src/refine_vcf.h src/vcf_hts.h \
    build/common_parallel.o build/vcf_hts.o $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g build/common_parallel.o \
	    build/vcf_hts.o src/refine_vcf.cpp -o $@ \
	    $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

utils/downsample_vcf_parallel: src/downsample_vcf_parallel.cpp \
    src/downsample_vcf_parallel.h build/common_parallel.o $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g build/common_parallel.o \
	    src/downsample_vcf_parallel.cpp -o $@ \
	    $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

utils/get_unique_kmers: src/get_unique_kmers.c src/FASTK/libfastk.c \
    src/FASTK/libfastk.h build/libfastk.o
	mkdir -p utils
	$(CCOMP) $(CIFLAGS) $(CFLAGS) build/libfastk.o src/get_unique_kmers.c \
	    -o $@ $(LFLAGS) -lz

utils/bam_indiv_rg: src/bam_indiv_rg.cpp build/common.o $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/bam_indiv_rg.cpp \
	    $(LFLAGS) $(DEPS) -o $@ $(DEPS2)

utils/bam_split_bcs: src/bam_split_bcs.cpp build/common.o $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/bam_split_bcs.cpp \
	    $(LFLAGS) $(DEPS) -o $@ $(DEPS2)

utils/bam_cb_cache_extract: src/bam_cb_cache_extract.cpp
	mkdir -p utils
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_CACHE) src/bam_cb_cache_extract.cpp \
	    -o $@ $(LFLAGS) $(DEPS_CACHE)

utils/atac_fq_preprocess: src/atac_fq_preprocess.cpp build/common.o $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) build/common.o \
	    src/atac_fq_preprocess.cpp $(LFLAGS) $(DEPS) -o $@ $(DEPS2)

utils/split_read_files: src/split_read_files.cpp build/common.o $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) build/common.o \
	    src/split_read_files.cpp $(LFLAGS) $(DEPS) -o $@ $(DEPS2)

utils/combine_species_counts: src/combine_species_counts.cpp build/common.o $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) build/common.o \
	    src/combine_species_counts.cpp $(LFLAGS) $(DEPS) -o $@ $(DEPS2)

utils/composite_bam2counts: src/composite_bam2counts.cpp $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) src/composite_bam2counts.cpp \
	    $(LFLAGS) $(DEPS) -o $@ $(DEPS2)

# -----------------------------------------------------------------------------
# Shared object modules
# -----------------------------------------------------------------------------

$(BUILD_DIR_STAMP):
	mkdir -p build
	touch $@

build/common.o: src/common.cpp src/common.h $(DEPS) | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/common.cpp -c -o $@

build/common_parallel.o: src/common.cpp src/common.h $(DEPS) \
    $(HTSWRAPPER_HEADER_STAMP) | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) src/common.cpp -c -o $@

build/io_parallel.o: src/io.cpp src/io.h src/common.h \
    $(HTSWRAPPER_HEADER_STAMP) | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) src/io.cpp -c -o $@

build/vcf_hts.o: src/vcf_hts.cpp src/vcf_hts.h src/common.h \
    lib/libhtswrapper.a $(HTSWRAPPER_HEADER_STAMP) | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) src/vcf_hts.cpp -c -o $@

build/genotype_llr.o: src/genotype_llr.cpp src/genotype_llr.h src/vcf_hts.h \
    src/common.h $(DEPS) $(HTSWRAPPER_HEADER_STAMP) | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) src/genotype_llr.cpp -c -o $@

build/ambient_rna_gex.o: src/ambient_rna_gex.cpp src/ambient_rna_gex.h \
    src/common.h $(DEPS) | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/ambient_rna_gex.cpp -c -o $@

build/species_kmers.o: src/species_kmers.cpp src/species_kmers.h src/common.h | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 -g src/species_kmers.cpp -c -o $@

build/reads_demux.o: src/reads_demux.cpp src/reads_demux.h src/common.h \
    lib/libhtswrapper.a | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 -g src/reads_demux.cpp -c -o $@

build/demux_species_io.o: src/demux_species_io.cpp src/demux_species_io.h \
    src/common.h lib/libhtswrapper.a | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g src/demux_species_io.cpp -c -o $@

build/libfastk.o: src/FASTK/libfastk.c src/FASTK/libfastk.h | $(BUILD_DIR_STAMP)
	$(CCOMP) $(CIFLAGS) $(CFLAGS) src/FASTK/libfastk.c -c -o $@

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

lib/libhtswrapper.a:
	mkdir -p dependencies/htswrapper/build dependencies/htswrapper/lib
	cd dependencies/htswrapper && $(MAKE) PREFIX=../.. BC_LENX2=$(BC_LENX2) KX2=$(KX2)
	cd dependencies/htswrapper && $(MAKE) install PREFIX=../..

$(HTSWRAPPER_HEADER_STAMP): lib/libhtswrapper.a
	mkdir -p include/htswrapper
	cd dependencies/htswrapper && $(MAKE) install PREFIX=../..
	test -s include/htswrapper/bc.h
	test -s include/htswrapper/bam.h
	test -s include/htswrapper/gzreader.h
	test -s include/htswrapper/robin_hood/robin_hood.h
	touch $@

lib/libmixturedist.a:
	mkdir -p dependencies/mixtureDist/build dependencies/mixtureDist/lib
	cd dependencies/mixtureDist && $(MAKE) PREFIX=../..
	cd dependencies/mixtureDist && $(MAKE) install PREFIX=../..

$(OPTIMML_PATCH_STAMP): $(OPTIMML_SOURCE_FILES) Makefile
	rm -rf $(OPTIMML_PATCH_DIR)
	mkdir -p build
	cp -a dependencies/optimML $(OPTIMML_PATCH_DIR)
	rm -rf $(OPTIMML_PATCH_DIR)/build $(OPTIMML_PATCH_DIR)/lib
	mkdir -p $(OPTIMML_PATCH_DIR)/build $(OPTIMML_PATCH_DIR)/lib
	sed -i \
	    -e 's/assert(-dot(g, p)<0);/if (-dot(g, p) >= 0) throw 1;/' \
	    -e 's/assert(-dot(g, p) < 0);/if (-dot(g, p) >= 0) throw 1;/' \
	    -e 's/assert(std::isfinite(alpha));/if (!std::isfinite(alpha)) throw 2;/' \
	    $(OPTIMML_PATCH_DIR)/src/stlbfgs/stlbfgs.cpp
	sed -i \
	    -e 's|mixcompsum_f += mixcompfracs_sparse\[i\]\[k\] / (exp(-x\[k\]) + 1);|mixcompsum_f += mixcompfracs_sparse[i][k] / (exp(-x[n_param - nmixcomp + k]) + 1);|' \
	    -e 's|mixcompsum_f += mixcompfracs\[i\]\[k\] / (exp(-x\[k\]) + 1);|mixcompsum_f += mixcompfracs[i][k] / (exp(-x[n_param - nmixcomp + k]) + 1);|' \
	    $(OPTIMML_PATCH_DIR)/src/multivar_ml.cpp
	sed -i \
	    -e 's|mixcompsum_f += mixcompfracs_sparse\[jid\]\[k\] / (exp(-x\[k\]) + 1);|mixcompsum_f += mixcompfracs_sparse[jid][k] / (exp(-x[n_param - nmixcomp + k]) + 1);|' \
	    $(OPTIMML_PATCH_DIR)/src/multivar.cpp
	grep -Fq 'if (-dot(g, p) >= 0) throw 1;' $(OPTIMML_PATCH_DIR)/src/stlbfgs/stlbfgs.cpp
	grep -Fq 'if (!std::isfinite(alpha)) throw 2;' $(OPTIMML_PATCH_DIR)/src/stlbfgs/stlbfgs.cpp
	grep -Fq 'mixcompsum_f += mixcompfracs_sparse[i][k] / (exp(-x[n_param - nmixcomp + k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar_ml.cpp
	grep -Fq 'mixcompsum_f += mixcompfracs[i][k] / (exp(-x[n_param - nmixcomp + k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar_ml.cpp
	grep -Fq 'mixcompsum_f += mixcompfracs_sparse[jid][k] / (exp(-x[n_param - nmixcomp + k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar.cpp
	touch $@

lib/liboptimml.a: $(OPTIMML_PATCH_STAMP)
	mkdir -p lib include/optimML
	$(MAKE) -C $(OPTIMML_PATCH_DIR) PREFIX=$(abspath .)
	$(MAKE) -C $(OPTIMML_PATCH_DIR) install PREFIX=$(abspath .)
	test -s $@

# -----------------------------------------------------------------------------
# Ambient mathematical tests
# -----------------------------------------------------------------------------

include/htswrapper/robin_hood/robin_hood.h: dependencies/htswrapper/src/robin_hood/robin_hood.h
	mkdir -p include/htswrapper/robin_hood
	cp $< $@

include/htswrapper/bc.h: dependencies/htswrapper/src/bc.h include/htswrapper/robin_hood/robin_hood.h
	mkdir -p include/htswrapper
	cp $< $@

build/ambient_rna_three_ap_test.o: check_required_ambient_source \
    src/ambient_rna_three_ap.h src/common.h lib/liboptimml.a lib/libmixturedist.a \
    include/htswrapper/bc.h include/htswrapper/robin_hood/robin_hood.h | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -ffunction-sections -fdata-sections \
	    src/ambient_rna_three_ap.cpp -c -o $@

build/test_ambient_support.o: src/test_ambient_support.cpp | $(BUILD_DIR_STAMP)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -ffunction-sections -fdata-sections \
	    src/test_ambient_support.cpp -c -o $@

test_ambient_math: src/test_ambient_gradients.cpp build/ambient_rna_three_ap_test.o \
    build/test_ambient_support.o lib/liboptimml.a lib/libmixturedist.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -ffunction-sections -fdata-sections \
	    src/test_ambient_gradients.cpp build/ambient_rna_three_ap_test.o \
	    build/test_ambient_support.o -o $@ $(LFLAGS_TET) -Wl,--gc-sections \
	    lib/liboptimml.a lib/libmixturedist.a -lz -lpthread
	./test_ambient_math

# -----------------------------------------------------------------------------
# Clean and install
# -----------------------------------------------------------------------------

clean: clean_build clean_binaries

clean_build:
	rm -f build/*.o

clean_binaries:
	rm -f $(ORCHESTRATOR_BINS) $(AUX_ROOT_BINS) $(DEPRECATED_BINS) test_ambient_math
	rm -f $(UTIL_BINS) utils/get_unique_kmers

clean_deps:
	cd dependencies/htswrapper && $(MAKE) clean || true
	cd dependencies/mixtureDist && $(MAKE) clean || true
	rm -rf $(OPTIMML_PATCH_DIR)

clean_all: clean clean_deps
	rm -f lib/libmixturedist.a lib/liboptimml.a lib/libhtswrapper.a

install: all
	mkdir -p $(PREFIX)/bin
	cp $(ORCHESTRATOR_BINS) $(AUX_ROOT_BINS) $(UTIL_BINS) $(PREFIX)/bin/
	@if [ -d scripts ]; then \
	    cp scripts/*.py $(PREFIX)/bin/ 2>/dev/null || true; \
	    chmod +x $(PREFIX)/bin/*.py 2>/dev/null || true; \
	fi

.PHONY: all orchestrator_tools auxiliary_tools utilities dependencies \
        quant3_contam quant3_contam_ap quant3_contam_empty_drops \
        source_inventory print_flags check_required_ambient_source \
        test_ambient_math clean clean_build clean_binaries clean_deps clean_all install
