SHELL = bash
COMP = g++
CCOMP = gcc
PREFIX ?= /usr/local

# -----------------------------------------------------------------------------
# Cluster HTSlib location
# -----------------------------------------------------------------------------
# Compile and link against the shared HTSlib installation used by the cluster,
# and embed its library directory as a runtime search path.  LIBRARY_PATH helps
# only at link time; without an rpath deployed executables fail at startup with:
#   libhts.so.3: cannot open shared object file
HTSLIB_PREFIX ?= /nvme/software/packages/htslib/1.20
HTSLIB_INCLUDE ?= $(HTSLIB_PREFIX)/include
HTSLIB_LIB ?= $(HTSLIB_PREFIX)/lib
HTSLIB_RPATH_FLAG = -Wl,-rpath,$(HTSLIB_LIB)

# -----------------------------------------------------------------------------
# Cluster-portable CPU target
# -----------------------------------------------------------------------------
# This tree is compiled once on ash and synced to pika/char/squirtle.
# ash/pika/char are Zen4-class EPYC systems, but squirtle is Zen3-class EPYC.
# Therefore the shared binaries must NOT use -march=native when compiled on ash.
# Default to Zen3 so the same binaries run safely on every node.
#
# Override only if intentionally building a node-specific binary, e.g.:
#   make CPU_ARCH=znver4 CPU_TUNE=znver4 ...
CPU_ARCH ?= znver3
CPU_TUNE ?= znver3
ARCHFLAGS ?= -march=$(CPU_ARCH) -mtune=$(CPU_TUNE)

# Standard flags for original CellBouncer tools
CXXFLAGS_STD = -std=c++11 -fPIC -D_REENTRANT -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2)

# Optimized flags for parallel tools:
# -O3 for high optimization without enabling fast-math semantics
# $(ARCHFLAGS) for a cluster-portable CPU target; default is znver3.
# -fopenmp for parallelism
CXXFLAGS_PARALLEL = -std=c++11 -fPIC -D_REENTRANT -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -DNBITS=$(NBITS) \
                    -O3 $(ARCHFLAGS) -fopenmp

CFLAGS = -fPIC -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -O3 $(ARCHFLAGS)
CXXIFLAGS = -I$(HTSLIB_INCLUDE) -I$(PREFIX)/include -Iinclude
CIFLAGS = -I$(HTSLIB_INCLUDE) -I$(PREFIX)/include -Iinclude
LFLAGS = -L$(HTSLIB_LIB) $(HTSLIB_RPATH_FLAG) -L$(PREFIX)/lib -Llib
LFLAGS_PARALLEL = -L$(HTSLIB_LIB) $(HTSLIB_RPATH_FLAG) -L$(PREFIX)/lib -Llib -fopenmp -flto=auto

# Optimized flags for tet contamination tools.  These deliberately use a
# separate ambient_rna_three_ap_tet.o object so legacy quant3_contam_ap remains
# buildable with the original standard flags, while the production tet tools get
# optimization and OpenMP linkage.
CXXFLAGS_TET = -std=c++11 -fPIC -D_REENTRANT -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) \
               -O3 $(ARCHFLAGS) -fopenmp
LFLAGS_TET = -L$(HTSLIB_LIB) $(HTSLIB_RPATH_FLAG) -L$(PREFIX)/lib -Llib -fopenmp

# Standalone barcode-cache extractor.  It uses std::filesystem, so it must
# be compiled as C++17 without changing the legacy CellBouncer C++11 tools.
CXXFLAGS_CACHE = -std=c++17 -fPIC -D_REENTRANT -O3 $(ARCHFLAGS) -Wall -Wextra -pedantic
DEPS_CACHE = -lhts -lz -lpthread

# Dedicated mitochondrial fusion-ratio estimator.  Mito panel selection itself
# is compiled directly into utils/downsample_vcf_parallel and is enabled there
# with --mt_output; there is no separate panel-selection executable.
CXXFLAGS_MT = -std=c++17 -fPIC -D_REENTRANT -O3 $(ARCHFLAGS) -Wall -Wextra -pedantic
DEPS_MT = -lhts -lz -lpthread

# Standalone genotype scrubber used when materializing synthetic benchmark BAMs.
# It uses htslib's codec thread pool but does not link CellBouncer objects.
CXXFLAGS_SCRUB = -std=c++11 -fPIC -D_REENTRANT -O3 $(ARCHFLAGS) -Wall -Wextra
DEPS_SCRUB = -lhts -lz -lpthread

NBITS ?= 2048

ifeq ($(findstring cellbouncer, ${CONDA_PREFIX}), cellbouncer)
    CXXIFLAGS += -I${CONDA_PREFIX}/include
    CIFLAGS += -I${CONDA_PREFIX}/include
    LFLAGS += -L${CONDA_PREFIX}/lib
    LFLAGS_PARALLEL += -L${CONDA_PREFIX}/lib
    LFLAGS_TET += -L${CONDA_PREFIX}/lib
endif

MAX_SITES ?= 2000
MAKE = make
PROJROOT = $(shell pwd)
BC_LENX2 = 32
KX2 = 16
DEPS = lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a
DEPS2 = -lz -lhts -lpthread
DEPS2_PARALLEL = -lz -lhts -lpthread -lrt

# optimML is an external dependency, but CellBouncer needs two STLBFGS
# assertions converted into catchable exceptions.  Never patch the vendored
# dependency tree in place: build an isolated patched copy under build/ so a
# clean checkout remains byte-for-byte unchanged.  The versioned stamp makes an
# already-existing lib/liboptimml.a rebuild whenever this patch recipe changes.
OPTIMML_PATCH_DIR = build/optimML_cellbouncer
OPTIMML_PATCH_STAMP = $(OPTIMML_PATCH_DIR)/.cellbouncer_patch_v3
OPTIMML_SOURCE_FILES = $(shell find dependencies/optimML -type f \
    -not -path '*/build/*' -not -path '*/lib/*')

# ============================================================================
# MAIN TARGETS
# ============================================================================

all: dependencies original_tools parallel_tools tet_tools qc_tools mitochondrial_tools

original_tools: demux_vcf demux_mt demux_tags demux_species quant_contam doublet_dragon bulkprops utils

parallel_tools: demux_parallel vcf_loader_daemon bam_ram_host_daemon genotype_scrub_bam utils/downsample_vcf_parallel tetra_refine

tet_tools: tet_ambient_profile tet_contam_estimate legacy2c_contam_estimate

qc_tools: tetra_score_calls snps_per_read

mitochondrial_tools: mt_fusion_ratio

utils: utils/refine_vcf utils/bam_indiv_rg utils/bam_split_bcs utils/bam_cb_cache_extract utils/get_unique_kmers utils/split_read_files utils/atac_fq_preprocess utils/combine_species_counts utils/composite_bam2counts utils/downsample_vcf

dependencies: lib/libhtswrapper.a lib/libmixturedist.a lib/liboptimml.a

print_flags:
	@echo "HTSLIB_PREFIX=$(HTSLIB_PREFIX)"
	@echo "HTSLIB_INCLUDE=$(HTSLIB_INCLUDE)"
	@echo "HTSLIB_LIB=$(HTSLIB_LIB)"
	@echo "HTSLIB_RPATH_FLAG=$(HTSLIB_RPATH_FLAG)"
	@echo "CPU_ARCH=$(CPU_ARCH)"
	@echo "CPU_TUNE=$(CPU_TUNE)"
	@echo "ARCHFLAGS=$(ARCHFLAGS)"
	@echo "CXXFLAGS_PARALLEL=$(CXXFLAGS_PARALLEL)"
	@echo "CXXFLAGS_TET=$(CXXFLAGS_TET)"
	@echo "CXXFLAGS_CACHE=$(CXXFLAGS_CACHE)"
	@echo "CXXFLAGS_SCRUB=$(CXXFLAGS_SCRUB)"
	@echo "DEPS_SCRUB=$(DEPS_SCRUB)"
	@echo "CFLAGS=$(CFLAGS)"

# ============================================================================
# ORIGINAL CELLBOUNCER TOOLS
# ============================================================================

demux_vcf: src/demux_vcf.cpp build/common.o build/demux_vcf_io.o build/demux_vcf_hts.o build/demux_vcf_llr.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g build/common.o build/demux_vcf_io.o build/demux_vcf_hts.o build/demux_vcf_llr.o src/demux_vcf.cpp -o demux_vcf $(LFLAGS) $(DEPS) $(DEPS2)

demux_mt: src/demux_mt.cpp src/common.h build/common.o build/demux_vcf_llr.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -D MAX_SITES=$(MAX_SITES) build/common.o build/demux_vcf_llr.o src/demux_mt.cpp -o demux_mt $(LFLAGS) $(DEPS) $(DEPS2)

demux_species: src/demux_species.cpp src/common.h build/common.o build/demux_species_io.o build/species_kmers.o build/reads_demux.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 -g build/common.o build/demux_species_io.o build/species_kmers.o build/reads_demux.o src/demux_species.cpp $(LFLAGS) $(DEPS) -pthread -o demux_species $(DEPS2)

demux_tags: src/demux_tags.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/demux_tags.cpp $(LFLAGS) $(DEPS) -D PROJ_ROOT=$(PROJROOT) -o demux_tags $(DEPS2)

quant_contam: src/common.h src/quant_contam.cpp src/ambient_rna.h build/common.o build/demux_vcf_io.o build/demux_vcf_llr.o build/ambient_rna.o build/ambient_rna_gex.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g build/common.o build/demux_vcf_io.o build/demux_vcf_llr.o build/ambient_rna.o build/ambient_rna_gex.o src/quant_contam.cpp $(LFLAGS) $(DEPS) -o quant_contam $(DEPS2)

# Compiled current-code two-component compatibility estimator.  This is a
# standalone C++ executable linked directly to the hardened ambient_rna model;
# it does not invoke or wrap quant_contam at runtime.
legacy2c_contam_estimate: src/legacy2c_contam_estimate.cpp src/ambient_rna.h src/common.h src/demux_vcf_io.h build/common.o build/demux_vcf_io.o build/demux_vcf_llr.o build/ambient_rna.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 $(ARCHFLAGS) -g build/common.o build/demux_vcf_io.o build/demux_vcf_llr.o build/ambient_rna.o src/legacy2c_contam_estimate.cpp $(LFLAGS) $(DEPS) -o legacy2c_contam_estimate $(DEPS2)

doublet_dragon: src/doublet_dragon.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/doublet_dragon.cpp $(LFLAGS) $(DEPS) -o doublet_dragon $(DEPS2)

bulkprops: src/bulkprops.cpp src/common.h build/common.o src/demux_vcf_hts.h build/demux_vcf_hts.o src/demux_vcf_io.h build/demux_vcf_io.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o build/demux_vcf_io.o build/demux_vcf_hts.o src/bulkprops.cpp $(LFLAGS) $(DEPS) -o bulkprops $(DEPS2)

# ============================================================================
# PARALLEL TOOLS
# ============================================================================

demux_parallel: src/demux_parallel.cpp build/common_parallel.o build/demux_vcf_io_parallel.o build/demux_parallel_hts.o build/demux_parallel_llr.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g build/common_parallel.o build/demux_vcf_io_parallel.o build/demux_parallel_hts.o build/demux_parallel_llr.o src/demux_parallel.cpp -o demux_parallel $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

vcf_loader_daemon: src/vcf_loader_daemon.cpp build/common_parallel.o build/demux_parallel_hts.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g build/common_parallel.o build/demux_parallel_hts.o src/vcf_loader_daemon.cpp -o vcf_loader_daemon $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

# bam_ram_host_daemon stages cache BAM files into node-local tmpfs (/dev/shm) as
# real files and publishes hosted_manifest.tsv. It is a standalone file-staging
# daemon (POSIX I/O + pthreads only, no htslib/CellBouncer deps), built C++11 with
# the cluster-portable ARCHFLAGS so the one binary runs on pika/char (Zen4) and
# squirtle (Zen3).
bam_ram_host_daemon: src/bam_ram_host_daemon.cpp
	$(COMP) -std=c++11 -O3 $(ARCHFLAGS) src/bam_ram_host_daemon.cpp -o bam_ram_host_daemon -lpthread

# Two-pass, per-cell dosage-projection scrubber for synthetic benchmark BAMs.
# Standalone executable: htslib + zlib + pthread only. Build through an object
# so the existing `make clean_build` workflow always forces a fresh rebuild.
build/genotype_scrub_bam.o: src/genotype_scrub_bam.cpp Makefile
	mkdir -p build
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_SCRUB) -c src/genotype_scrub_bam.cpp \
	    -o build/genotype_scrub_bam.o

genotype_scrub_bam: build/genotype_scrub_bam.o
	$(COMP) build/genotype_scrub_bam.o -o genotype_scrub_bam \
	    $(LFLAGS) $(DEPS_SCRUB)

utils/downsample_vcf_parallel: src/downsample_vcf_parallel.cpp src/downsample_vcf_parallel.h build/common_parallel.o $(DEPS)
	mkdir -p utils
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g build/common_parallel.o src/downsample_vcf_parallel.cpp -o utils/downsample_vcf_parallel $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

mt_fusion_ratio: src/mt_fusion_ratio.cpp
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_MT) src/mt_fusion_ratio.cpp -o mt_fusion_ratio $(LFLAGS) $(DEPS_MT)

tetra_refine: src/tetra_refine.cpp lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 src/tetra_refine.cpp -o tetra_refine $(LFLAGS) lib/libhtswrapper.a -lz

# ============================================================================
# THREE-COMPONENT MODEL (quant3_contam) - legacy
# ============================================================================

build/ambient_rna_three.o: src/ambient_rna_three.cpp src/ambient_rna_three.h src/common.h $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/ambient_rna_three.cpp -c -o build/ambient_rna_three.o

quant3_contam: src/quant3_contam.cpp src/ambient_rna_three.h src/common.h build/common.o build/demux_vcf_io.o build/demux_vcf_llr.o build/ambient_rna_three.o build/ambient_rna_gex.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g build/common.o build/demux_vcf_io.o build/demux_vcf_llr.o build/ambient_rna_three.o build/ambient_rna_gex.o src/quant3_contam.cpp -o quant3_contam $(LFLAGS) $(DEPS) $(DEPS2)

# ============================================================================
# THREE-COMPONENT MODEL + AP REFINEMENTS (quant3_contam_ap) - legacy
# These remain buildable for comparison but are superseded by tet_tools.
# ============================================================================

build/ambient_rna_three_ap.o: src/ambient_rna_three_ap.cpp src/ambient_rna_three_ap.h src/common.h $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g src/ambient_rna_three_ap.cpp -c -o build/ambient_rna_three_ap.o

build/demux_vcf_io_species.o: src/demux_vcf_io_species.cpp src/demux_vcf_io.h src/common.h
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g src/demux_vcf_io_species.cpp -c -o build/demux_vcf_io_species.o

quant3_contam_ap: src/quant3_contam_ap.cpp src/ambient_rna_three_ap.h src/common.h build/common.o build/demux_vcf_io.o build/demux_vcf_io_species.o build/demux_vcf_llr.o build/ambient_rna_three_ap.o build/ambient_rna_gex.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g build/common.o build/demux_vcf_io.o build/demux_vcf_io_species.o build/demux_vcf_llr.o build/ambient_rna_three_ap.o build/ambient_rna_gex.o src/quant3_contam_ap.cpp $(LFLAGS) $(DEPS) -o quant3_contam_ap $(DEPS2)

quant3_contam_empty_drops: src/quant3_contam_empty_drops.cpp src/ambient_rna_three_ap.h src/common.h build/common.o build/demux_vcf_io.o build/demux_vcf_io_species.o build/demux_vcf_llr.o build/ambient_rna_three_ap.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g build/common.o build/demux_vcf_io.o build/demux_vcf_io_species.o build/demux_vcf_llr.o build/ambient_rna_three_ap.o src/quant3_contam_empty_drops.cpp $(LFLAGS) $(DEPS) -o quant3_contam_empty_drops $(DEPS2)

# ============================================================================
# TET TOOLS (new pipeline, replaces quant3_contam_ap / quant3_contam_empty_drops)
# Same shared library (ambient_rna_three_ap) with WP0 fixes applied.
# ============================================================================

build/ambient_rna_three_ap_tet.o: src/ambient_rna_three_ap.cpp src/ambient_rna_three_ap.h src/common.h $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -g src/ambient_rna_three_ap.cpp -c -o build/ambient_rna_three_ap_tet.o

tet_ambient_profile: src/tet_ambient_profile.cpp src/ambient_rna_three_ap.h src/common.h \
    build/common.o build/demux_vcf_io.o build/demux_vcf_io_species.o build/demux_vcf_llr.o \
    build/ambient_rna_three_ap_tet.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -g build/common.o build/demux_vcf_io.o \
	    build/demux_vcf_io_species.o build/demux_vcf_llr.o build/ambient_rna_three_ap_tet.o \
	    src/tet_ambient_profile.cpp $(LFLAGS_TET) $(DEPS) -o tet_ambient_profile $(DEPS2)

tet_contam_estimate: src/tet_contam_estimate.cpp src/ambient_rna_three_ap.h src/common.h \
    build/common.o build/demux_vcf_io.o build/demux_vcf_io_species.o build/demux_vcf_llr.o \
    build/ambient_rna_three_ap_tet.o build/ambient_rna_gex.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -g build/common.o build/demux_vcf_io.o \
	    build/demux_vcf_io_species.o build/demux_vcf_llr.o build/ambient_rna_three_ap_tet.o \
	    build/ambient_rna_gex.o src/tet_contam_estimate.cpp $(LFLAGS_TET) $(DEPS) \
	    -o tet_contam_estimate $(DEPS2)

# Per-cell post-hoc QC scorer used by the swap-audit layer.
tetra_score_calls: src/tetra_score_calls.cpp lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -g src/tetra_score_calls.cpp \
	    -o tetra_score_calls $(LFLAGS_TET) lib/libhtswrapper.a -lz

# Standalone per-read SNP-coverage counter: htslib + zlib + pthread + OpenMP.
# robin_hood is header-only, resolved via -Iinclude from the installed htswrapper
# headers. Built as C++17 like bam_cb_cache_extract. The order-only
# lib/libhtswrapper.a prerequisite installs include/htswrapper/robin_hood/robin_hood.h;
# the archive itself is not linked.
snps_per_read: src/snps_per_read.cpp | lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) -std=c++17 -fPIC -D_REENTRANT -O3 $(ARCHFLAGS) -fopenmp -Wall -Wextra \
	    src/snps_per_read.cpp -o snps_per_read $(LFLAGS) -lhts -lz -lpthread

# ============================================================================
# UTILITY TOOLS
# ============================================================================

utils/refine_vcf: src/refine_vcf.cpp src/refine_vcf.h src/common.h build/common.o src/demux_vcf_hts.h build/demux_vcf_hts.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g build/common.o build/demux_vcf_hts.o src/refine_vcf.cpp $(LFLAGS) $(DEPS) -o utils/refine_vcf $(DEPS2)

utils/bam_indiv_rg: src/bam_indiv_rg.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/bam_indiv_rg.cpp $(LFLAGS) $(DEPS) -o utils/bam_indiv_rg $(DEPS2)

utils/bam_split_bcs: src/bam_split_bcs.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/bam_split_bcs.cpp $(LFLAGS) $(DEPS) -o utils/bam_split_bcs $(DEPS2)

utils/bam_cb_cache_extract: src/bam_cb_cache_extract.cpp
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_CACHE) src/bam_cb_cache_extract.cpp -o utils/bam_cb_cache_extract $(LFLAGS) $(DEPS_CACHE)

utils/get_unique_kmers: src/get_unique_kmers.c src/FASTK/libfastk.c build/libfastk.o
	$(CCOMP) $(CIFLAGS) $(CFLAGS) build/libfastk.o src/get_unique_kmers.c -o utils/get_unique_kmers $(LFLAGS) -lz

utils/atac_fq_preprocess: src/atac_fq_preprocess.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) build/common.o src/atac_fq_preprocess.cpp $(LFLAGS) $(DEPS) -o utils/atac_fq_preprocess $(DEPS2)

utils/split_read_files: src/split_read_files.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) build/common.o src/split_read_files.cpp $(LFLAGS) $(DEPS) -o utils/split_read_files $(DEPS2)

utils/combine_species_counts: src/combine_species_counts.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) build/common.o src/combine_species_counts.cpp $(LFLAGS) $(DEPS) -o utils/combine_species_counts $(DEPS2)

utils/composite_bam2counts: src/composite_bam2counts.cpp lib/libhtswrapper.a $(DEPS)
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) src/composite_bam2counts.cpp $(LFLAGS) $(DEPS) -o utils/composite_bam2counts $(DEPS2)

utils/downsample_vcf: src/downsample_vcf.cpp src/downsample_vcf.h build/common.o $(DEPS)
	$(COMP) $(CXXFLAGS_STD) $(CXXIFLAGS) -DNBITS=$(NBITS) src/downsample_vcf.cpp build/common.o $(LFLAGS) $(DEPS) -o utils/downsample_vcf $(DEPS2)

# ============================================================================
# OBJECT FILES - ORIGINAL (standard flags)
# ============================================================================

build/common.o: src/common.cpp src/common.h lib/libhtswrapper.a lib/libmixturedist.a lib/liboptimml.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/common.cpp -c -o build/common.o

build/demux_vcf_io.o: src/demux_vcf_io.cpp src/demux_vcf_io.h src/common.h
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/demux_vcf_io.cpp -c -o build/demux_vcf_io.o

build/demux_vcf_hts.o: src/demux_vcf_hts.cpp src/demux_vcf_hts.h src/common.h lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/demux_vcf_hts.cpp -c -o build/demux_vcf_hts.o

build/demux_vcf_llr.o: src/demux_vcf_llr.cpp src/demux_vcf_llr.h src/common.h
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/demux_vcf_llr.cpp -c -o build/demux_vcf_llr.o

build/ambient_rna.o: src/ambient_rna.cpp src/ambient_rna.h src/common.h $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/ambient_rna.cpp -c -o build/ambient_rna.o

build/ambient_rna_gex.o: src/ambient_rna_gex.cpp src/ambient_rna_gex.h src/common.h $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) src/ambient_rna_gex.cpp -c -o build/ambient_rna_gex.o 

build/species_kmers.o: src/species_kmers.cpp src/species_kmers.h src/common.h
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 -g src/species_kmers.cpp -c -o build/species_kmers.o

build/reads_demux.o: src/reads_demux.cpp src/reads_demux.h src/common.h lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -O3 -g src/reads_demux.cpp -c -o build/reads_demux.o

build/demux_species_io.o: src/demux_species_io.cpp src/demux_species_io.h src/common.h lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g src/demux_species_io.cpp -c -o build/demux_species_io.o

build/libfastk.o: src/FASTK/libfastk.c src/FASTK/libfastk.h
	$(CCOMP) $(CIFLAGS) $(CFLAGS) src/FASTK/libfastk.c -c -o build/libfastk.o

build/gene_core.o: src/FASTK/gene_core.c src/FASTK/gene_core.h
	$(CCOMP) $(CIFLAGS) $(CFLAGS) src/FASTK/gene_core.c -c -o build/gene_core.o

# ============================================================================
# OBJECT FILES - PARALLEL (optimized flags with OpenMP)
# ============================================================================

build/common_parallel.o: src/common.cpp src/common.h lib/libhtswrapper.a lib/libmixturedist.a lib/liboptimml.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) src/common.cpp -c -o build/common_parallel.o

build/demux_vcf_io_parallel.o: src/demux_vcf_io.cpp src/demux_vcf_io.h src/common.h
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) src/demux_vcf_io.cpp -c -o build/demux_vcf_io_parallel.o

build/demux_parallel_hts.o: src/demux_parallel_hts.cpp src/demux_parallel_hts.h src/common.h lib/libhtswrapper.a
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) src/demux_parallel_hts.cpp -c -o build/demux_parallel_hts.o

build/demux_parallel_llr.o: src/demux_parallel_llr.cpp src/demux_parallel_llr.h src/common.h
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) src/demux_parallel_llr.cpp -c -o build/demux_parallel_llr.o

# ============================================================================
# DEPENDENCIES
#
# optimML/STLBFGS handling is intentionally implemented only here, in the
# top-level build.  The vendored dependencies/optimML source is never edited.
# A private copy is created in build/, patched there, compiled, and installed
# into this project's lib/ and include/ directories.
# ============================================================================

lib/libhtswrapper.a:
	mkdir -p dependencies/htswrapper/build dependencies/htswrapper/lib
	cd dependencies/htswrapper && $(MAKE) PREFIX=../.. BC_LENX2=$(BC_LENX2) KX2=$(KX2)
	cd dependencies/htswrapper && $(MAKE) install PREFIX=../..

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
	# CellBouncer requires optimML's mixture logistic coordinates to use the
	# actual mixture block when one or more external parameters precede it.
	# Apply this only to the private build copy; never edit dependencies/optimML.
	sed -i \
	    -e 's|mixcompsum_f += mixcompfracs_sparse\[i\]\[k\] / (exp(-x\[k\]) + 1);|mixcompsum_f += mixcompfracs_sparse[i][k] / (exp(-x[n_param - nmixcomp + k]) + 1);|' \
	    -e 's|mixcompsum_f += mixcompfracs\[i\]\[k\] / (exp(-x\[k\]) + 1);|mixcompsum_f += mixcompfracs[i][k] / (exp(-x[n_param - nmixcomp + k]) + 1);|' \
	    $(OPTIMML_PATCH_DIR)/src/multivar_ml.cpp
	sed -i \
	    -e 's|mixcompsum_f += mixcompfracs_sparse\[jid\]\[k\] / (exp(-x\[k\]) + 1);|mixcompsum_f += mixcompfracs_sparse[jid][k] / (exp(-x[n_param - nmixcomp + k]) + 1);|' \
	    $(OPTIMML_PATCH_DIR)/src/multivar.cpp
	grep -Fq 'if (-dot(g, p) >= 0) throw 1;' $(OPTIMML_PATCH_DIR)/src/stlbfgs/stlbfgs.cpp
	grep -Fq 'if (!std::isfinite(alpha)) throw 2;' $(OPTIMML_PATCH_DIR)/src/stlbfgs/stlbfgs.cpp
	! grep -Fq 'assert(std::isfinite(alpha));' $(OPTIMML_PATCH_DIR)/src/stlbfgs/stlbfgs.cpp
	grep -Fq 'mixcompsum_f += mixcompfracs_sparse[i][k] / (exp(-x[n_param - nmixcomp + k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar_ml.cpp
	grep -Fq 'mixcompsum_f += mixcompfracs[i][k] / (exp(-x[n_param - nmixcomp + k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar_ml.cpp
	grep -Fq 'mixcompsum_f += mixcompfracs_sparse[jid][k] / (exp(-x[n_param - nmixcomp + k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar.cpp
	! grep -Fq 'mixcompsum_f += mixcompfracs_sparse[i][k] / (exp(-x[k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar_ml.cpp
	! grep -Fq 'mixcompsum_f += mixcompfracs[i][k] / (exp(-x[k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar_ml.cpp
	! grep -Fq 'mixcompsum_f += mixcompfracs_sparse[jid][k] / (exp(-x[k]) + 1);' $(OPTIMML_PATCH_DIR)/src/multivar.cpp
	touch $@

lib/liboptimml.a: $(OPTIMML_PATCH_STAMP)
	mkdir -p lib include/optimML
	$(MAKE) -C $(OPTIMML_PATCH_DIR) PREFIX=$(abspath .)
	$(MAKE) -C $(OPTIMML_PATCH_DIR) install PREFIX=$(abspath .)
	test -s $@


# ============================================================================
# AMBIENT MODEL MATHEMATICAL TESTS
# ============================================================================

# Stage only the headers required by the isolated math-test build without
# compiling or modifying the htswrapper dependency.
include/htswrapper/robin_hood/robin_hood.h: dependencies/htswrapper/src/robin_hood/robin_hood.h
	mkdir -p include/htswrapper/robin_hood
	cp $< $@

include/htswrapper/bc.h: dependencies/htswrapper/src/bc.h include/htswrapper/robin_hood/robin_hood.h
	mkdir -p include/htswrapper
	cp $< $@

build/ambient_rna_three_ap_test.o: src/ambient_rna_three_ap.cpp src/ambient_rna_three_ap.h src/common.h lib/liboptimml.a lib/libmixturedist.a include/htswrapper/bc.h include/htswrapper/robin_hood/robin_hood.h
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -ffunction-sections -fdata-sections src/ambient_rna_three_ap.cpp -c -o build/ambient_rna_three_ap_test.o

build/test_ambient_support.o: src/test_ambient_support.cpp
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_TET) -ffunction-sections -fdata-sections src/test_ambient_support.cpp -c -o build/test_ambient_support.o

test_ambient_math: src/test_ambient_gradients.cpp build/ambient_rna_three_ap_test.o build/test_ambient_support.o lib/liboptimml.a lib/libmixturedist.a
	$(COMP) $(CXXIFLAGS) -Isrc $(CXXFLAGS_TET) -ffunction-sections -fdata-sections src/test_ambient_gradients.cpp build/ambient_rna_three_ap_test.o build/test_ambient_support.o -o test_ambient_math $(LFLAGS_TET) -Wl,--gc-sections lib/liboptimml.a lib/libmixturedist.a -lz -lpthread
	./test_ambient_math

# ============================================================================
# CLEAN TARGETS
# ============================================================================

clean: clean_build clean_binaries

clean_build:
	rm -f build/*.o
	rm -rf $(OPTIMML_PATCH_DIR)
	rm -f lib/liboptimml.a
	rm -rf include/optimML

clean_binaries:
	rm -f demux_vcf demux_mt demux_species demux_tags quant_contam doublet_dragon bulkprops
	rm -f demux_parallel vcf_loader_daemon bam_ram_host_daemon genotype_scrub_bam tetra_refine
	rm -f quant3_contam quant3_contam_ap quant3_contam_empty_drops
	rm -f tet_ambient_profile tet_contam_estimate legacy2c_contam_estimate tetra_score_calls snps_per_read test_ambient_math mt_fusion_ratio
	rm -f utils/refine_vcf utils/bam_indiv_rg utils/bam_split_bcs utils/bam_cb_cache_extract utils/get_unique_kmers
	rm -f utils/split_read_files utils/atac_fq_preprocess utils/combine_species_counts
	rm -f utils/composite_bam2counts utils/downsample_vcf utils/downsample_vcf_parallel

clean_deps:
	cd dependencies/htswrapper && $(MAKE) clean || true
	cd dependencies/mixtureDist && $(MAKE) clean || true
	rm -rf $(OPTIMML_PATCH_DIR)

clean_all: clean clean_deps
	rm -f lib/libmixturedist.a lib/liboptimml.a lib/libhtswrapper.a

# ============================================================================
# INSTALL
# ============================================================================

install: all quant3_contam quant3_contam_ap quant3_contam_empty_drops install_scripts
	mkdir -p $(PREFIX)/bin
	cp demux_vcf demux_mt demux_species demux_tags quant_contam doublet_dragon bulkprops $(PREFIX)/bin/
	cp demux_parallel vcf_loader_daemon bam_ram_host_daemon genotype_scrub_bam tetra_refine $(PREFIX)/bin/
	cp quant3_contam quant3_contam_ap quant3_contam_empty_drops $(PREFIX)/bin/
	cp tet_ambient_profile tet_contam_estimate legacy2c_contam_estimate tetra_score_calls snps_per_read mt_fusion_ratio $(PREFIX)/bin/
	cp utils/refine_vcf utils/bam_indiv_rg utils/bam_split_bcs utils/bam_cb_cache_extract utils/get_unique_kmers $(PREFIX)/bin/
	cp utils/split_read_files utils/atac_fq_preprocess utils/combine_species_counts $(PREFIX)/bin/
	cp utils/composite_bam2counts utils/downsample_vcf utils/downsample_vcf_parallel $(PREFIX)/bin/


install_scripts:
	mkdir -p $(PREFIX)/bin
	cp scripts/*.py $(PREFIX)/bin/
	chmod +x $(PREFIX)/bin/*.py

.PHONY: all original_tools parallel_tools tet_tools qc_tools mitochondrial_tools utils dependencies print_flags test_ambient_math clean clean_build clean_binaries clean_deps clean_all install install_scripts
