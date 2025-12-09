SHELL = bash
COMP = g++
CCOMP = gcc
PREFIX ?= /usr/local

# Standard flags for original CellBouncer tools
CXXFLAGS_STD = -std=c++11 -fPIC -D_REENTRANT -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2)

# Optimized flags for parallel tools: -Ofast for aggressive optimization, 
# -march=native for CPU-specific instructions, -ffast-math for faster floating point, 
# -fopenmp for parallelism
CXXFLAGS_PARALLEL = -std=c++11 -fPIC -D_REENTRANT -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -DNBITS=$(NBITS) \
                    -Ofast -march=native -ffast-math -fopenmp

CFLAGS = -fPIC -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -O3 -march=native
CXXIFLAGS = -I$(PREFIX)/include -Iinclude
CIFLAGS = -I$(PREFIX)/include -Iinclude
LFLAGS = -L$(PREFIX)/lib -Llib
LFLAGS_PARALLEL = -L$(PREFIX)/lib -Llib -fopenmp -flto
NBITS ?= 2048

ifeq ($(findstring cellbouncer, ${CONDA_PREFIX}), cellbouncer)
    CXXIFLAGS += -I${CONDA_PREFIX}/include
    CIFLAGS += -I${CONDA_PREFIX}/include
    LFLAGS += -L${CONDA_PREFIX}/lib
    LFLAGS_PARALLEL += -L${CONDA_PREFIX}/lib
endif

MAX_SITES ?= 2000
MAKE = make
PROJROOT = $(shell pwd)
BC_LENX2 = 32
KX2 = 16
DEPS = lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a
DEPS2 = -lz -lhts -lpthread
DEPS2_PARALLEL = -lz -lhts -lpthread -lrt

# ============================================================================
# MAIN TARGETS
# ============================================================================

all: dependencies original_tools parallel_tools

original_tools: demux_vcf demux_mt demux_tags demux_species quant_contam doublet_dragon bulkprops utils

parallel_tools: demux_parallel vcf_loader_daemon

utils: utils/refine_vcf utils/bam_indiv_rg utils/bam_split_bcs utils/get_unique_kmers utils/split_read_files utils/atac_fq_preprocess utils/combine_species_counts utils/composite_bam2counts utils/downsample_vcf

dependencies: lib/libhtswrapper.a lib/libmixturedist.a lib/liboptimml.a

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

doublet_dragon: src/doublet_dragon.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/doublet_dragon.cpp $(LFLAGS) $(DEPS) -o doublet_dragon $(DEPS2)

bulkprops: src/bulkprops.cpp src/common.h build/common.o src/demux_vcf_hts.h build/demux_vcf_hts.o src/demux_vcf_io.h build/demux_vcf_io.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o build/demux_vcf_io.o build/demux_vcf_hts.o src/bulkprops.cpp $(LFLAGS) $(DEPS) -o bulkprops $(DEPS2)

# ============================================================================
# PARALLEL TOOLS (NEW)
# ============================================================================

demux_parallel: src/demux_parallel.cpp build/common_parallel.o build/demux_vcf_io_parallel.o build/demux_parallel_hts.o build/demux_parallel_llr.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g build/common_parallel.o build/demux_vcf_io_parallel.o build/demux_parallel_hts.o build/demux_parallel_llr.o src/demux_parallel.cpp -o demux_parallel $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

vcf_loader_daemon: src/vcf_loader_daemon.cpp build/common_parallel.o build/demux_parallel_hts.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_PARALLEL) -g build/common_parallel.o build/demux_parallel_hts.o src/vcf_loader_daemon.cpp -o vcf_loader_daemon $(LFLAGS_PARALLEL) $(DEPS) $(DEPS2_PARALLEL)

# ============================================================================
# UTILITY TOOLS
# ============================================================================

utils/refine_vcf: src/refine_vcf.cpp src/refine_vcf.h src/common.h build/common.o src/demux_vcf_hts.h build/demux_vcf_hts.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) -g build/common.o build/demux_vcf_hts.o src/refine_vcf.cpp $(LFLAGS) $(DEPS) -o utils/refine_vcf $(DEPS2)

utils/bam_indiv_rg: src/bam_indiv_rg.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/bam_indiv_rg.cpp $(LFLAGS) $(DEPS) -o utils/bam_indiv_rg $(DEPS2)

utils/bam_split_bcs: src/bam_split_bcs.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(CXXIFLAGS) $(CXXFLAGS_STD) build/common.o src/bam_split_bcs.cpp $(LFLAGS) $(DEPS) -o utils/bam_split_bcs $(DEPS2)

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
# ============================================================================

lib/libhtswrapper.a:
	cd dependencies/htswrapper && $(MAKE) PREFIX=../.. BC_LENX2=$(BC_LENX2) KX2=$(KX2)
	cd dependencies/htswrapper && $(MAKE) install PREFIX=../..

lib/libmixturedist.a:
	cd dependencies/mixtureDist && $(MAKE) PREFIX=../..
	cd dependencies/mixtureDist && $(MAKE) install PREFIX=../..

lib/liboptimml.a:
	cd dependencies/optimML && $(MAKE) PREFIX=../..
	cd dependencies/optimML && $(MAKE) install PREFIX=../..

# ============================================================================
# CLEAN TARGETS
# ============================================================================

clean: clean_build clean_binaries

clean_build:
	rm -f build/*.o

clean_binaries:
	rm -f demux_vcf demux_mt demux_species demux_tags quant_contam doublet_dragon bulkprops
	rm -f demux_parallel vcf_loader_daemon
	rm -f utils/refine_vcf utils/bam_indiv_rg utils/bam_split_bcs utils/get_unique_kmers
	rm -f utils/split_read_files utils/atac_fq_preprocess utils/combine_species_counts
	rm -f utils/composite_bam2counts utils/downsample_vcf

clean_deps:
	cd dependencies/htswrapper && $(MAKE) clean || true
	cd dependencies/mixtureDist && $(MAKE) clean || true
	cd dependencies/optimML && $(MAKE) clean || true

clean_all: clean clean_deps
	rm -f lib/libmixturedist.a lib/liboptimml.a lib/libhtswrapper.a

# ============================================================================
# INSTALL
# ============================================================================

install: all
	mkdir -p $(PREFIX)/bin
	cp demux_vcf demux_mt demux_species demux_tags quant_contam doublet_dragon bulkprops $(PREFIX)/bin/
	cp demux_parallel vcf_loader_daemon $(PREFIX)/bin/
	cp utils/refine_vcf utils/bam_indiv_rg utils/bam_split_bcs utils/get_unique_kmers $(PREFIX)/bin/
	cp utils/split_read_files utils/atac_fq_preprocess utils/combine_species_counts $(PREFIX)/bin/
	cp utils/composite_bam2counts utils/downsample_vcf $(PREFIX)/bin/

.PHONY: all original_tools parallel_tools utils dependencies clean clean_build clean_binaries clean_deps clean_all install
