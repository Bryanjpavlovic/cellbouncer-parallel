#!/bin/bash
#SBATCH --job-name=downsample_vcf
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=downsample_vcf_%j.out
#SBATCH --error=downsample_vcf_%j.err

##############################################################################
# Unified VCF Downsampling Pipeline (Parallel)
#
# Runs downsample_vcf_parallel with all scoring in one step:
#   - DEMUX_SCORE: Clade rarity score (parallel by chromosome)
#   - ANNOT_SCORE: Gene body annotation (1.0 genic, 0.1 intergenic)
#   - COV_SCORE: RNA-seq coverage at each position
#
# Uses RAM disk for optimal I/O performance with large VCF files.
##############################################################################

set -euo pipefail

echo "=============================================="
echo "Unified VCF Downsampling Pipeline (Parallel)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_JOB_NODELIST}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo "=============================================="

# === CONFIGURATION (EDIT THESE) ===
VCF_IN="/path/to/your/input.vcf.gz"
COV_FILE="/path/to/your/coverage.bedgraph.gz"
GTF_FILE="/path/to/your/annotations.gtf.gz"

# Target number of SNPs in output
TARGET_NUM=20000000

# Output directory
OUT_DIR="./output"

# Number of threads (should match cpus-per-task)
THREADS=96

# Random seed for reproducibility (comment out for random)
SEED=42

# RAM disk configuration
RAMDISK="/dev/shm/vcf_downsample_${SLURM_JOB_ID}"
# ======================

# Human-readable number for filename
if [[ ${TARGET_NUM} -ge 1000000 ]]; then
    NUM_STR="$((TARGET_NUM / 1000000))M"
elif [[ ${TARGET_NUM} -ge 1000 ]]; then
    NUM_STR="$((TARGET_NUM / 1000))K"
else
    NUM_STR="${TARGET_NUM}"
fi

OUTPUT_VCF="${OUT_DIR}/downsampled_${NUM_STR}.vcf.gz"

# Create directories
mkdir -p "${OUT_DIR}"

# Load modules (adjust for your cluster)
# module purge
# module load htslib/1.20
# module load bcftools/1.20

echo ""
echo "Configuration:"
echo "  Input VCF: ${VCF_IN}"
echo "  GTF file: ${GTF_FILE}"
echo "  Coverage file: ${COV_FILE}"
echo "  Target SNPs: ${TARGET_NUM} (${NUM_STR})"
echo "  Threads: ${THREADS}"
echo "  Output: ${OUTPUT_VCF}"
echo "  RAM disk: ${RAMDISK}"
echo ""

# Verify inputs exist
for f in "${VCF_IN}" "${COV_FILE}" "${GTF_FILE}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: File not found: ${f}"
        exit 1
    fi
done

# Check for VCF index (needed for bcftools view to convert efficiently)
if [[ ! -f "${VCF_IN}.tbi" ]] && [[ ! -f "${VCF_IN}.csi" ]]; then
    echo "WARNING: No index found for input VCF, conversion may be slower"
fi

# ============================================================================
# Setup RAM disk
# ============================================================================
echo "Setting up RAM disk..."

mkdir -p "${RAMDISK}"

# Cleanup function
cleanup() {
    echo ""
    echo "Cleaning up RAM disk..."
    rm -rf "${RAMDISK}"
    echo "Cleanup complete."
}
trap cleanup EXIT

echo "  RAM disk created at: ${RAMDISK}"
echo "  Available space: $(df -h "${RAMDISK}" | tail -1 | awk '{print $4}')"

# ============================================================================
# Copy files to RAM disk
# ============================================================================
echo ""
echo "Copying files to RAM disk..."

START_COPY=$(date +%s)

# Copy VCF and convert to BCF for parallel queries
# BCF is cached on persistent storage to avoid reconversion on subsequent runs

BCF_CACHE="${VCF_IN%.vcf.gz}.bcf"
BCF_CACHE_INDEX="${BCF_CACHE}.csi"

if [[ -f "${BCF_CACHE}" ]] && [[ -f "${BCF_CACHE_INDEX}" ]]; then
    echo "  Found cached BCF: ${BCF_CACHE}"
    echo "  Copying BCF to ramdisk..."
    cp "${BCF_CACHE}" "${RAMDISK}/input.bcf"
    cp "${BCF_CACHE_INDEX}" "${RAMDISK}/input.bcf.csi"
else
    echo "  No cached BCF found, converting VCF.gz to BCF..."
    echo "  This may take several minutes (will be cached for future runs)..."
    
    # Convert VCF.gz to BCF on ramdisk (BCF supports proper random access)
    bcftools view -Ob -o "${RAMDISK}/input.bcf" --threads "${THREADS}" "${VCF_IN}"
    
    echo "  Creating CSI index for BCF..."
    bcftools index -c "${RAMDISK}/input.bcf"
    
    echo "  Saving BCF to persistent storage for future runs..."
    cp "${RAMDISK}/input.bcf" "${BCF_CACHE}"
    cp "${RAMDISK}/input.bcf.csi" "${BCF_CACHE_INDEX}"
    echo "    -> Cached: ${BCF_CACHE}"
fi

echo "    -> input.bcf ready"

# Verify BCF exists
if [[ ! -f "${RAMDISK}/input.bcf" ]] || [[ ! -f "${RAMDISK}/input.bcf.csi" ]]; then
    echo "ERROR: BCF not found on ramdisk!"
    ls -la "${RAMDISK}/"
    exit 1
fi

# Copy coverage file and index
echo "  Copying coverage file..."
cp "${COV_FILE}" "${RAMDISK}/coverage.bg.gz"
if [[ -f "${COV_FILE}.tbi" ]]; then
    cp "${COV_FILE}.tbi" "${RAMDISK}/coverage.bg.gz.tbi"
fi

# Copy GTF (cache decompressed version for future runs)
echo "  Preparing GTF file..."
GTF_CACHE="${GTF_FILE%.gz}"

if [[ "${GTF_FILE}" == *.gz ]]; then
    if [[ -f "${GTF_CACHE}" ]]; then
        echo "    Found cached decompressed GTF: ${GTF_CACHE}"
        cp "${GTF_CACHE}" "${RAMDISK}/annotations.gtf"
    else
        echo "    Decompressing GTF (will be cached for future runs)..."
        zcat "${GTF_FILE}" > "${RAMDISK}/annotations.gtf"
        cp "${RAMDISK}/annotations.gtf" "${GTF_CACHE}"
        echo "    -> Cached: ${GTF_CACHE}"
    fi
else
    cp "${GTF_FILE}" "${RAMDISK}/annotations.gtf"
fi

END_COPY=$(date +%s)
echo "  Copy complete in $((END_COPY - START_COPY)) seconds"
echo "  RAM disk usage: $(du -sh "${RAMDISK}" | cut -f1)"
echo ""
echo "  RAM disk contents:"
ls -lh "${RAMDISK}/"

# ============================================================================
# Run parallel downsampling
# ============================================================================
echo ""
echo "=============================================="
echo "Running downsample_vcf_parallel..."
echo "=============================================="

# Count input variants
echo "Input variants: $(bcftools index -n "${RAMDISK}/input.bcf")"
echo ""

START_RUN=$(date +%s)

downsample_vcf_parallel \
    --vcf "${RAMDISK}/input.bcf" \
    --num "${TARGET_NUM}" \
    --output "${RAMDISK}/output.vcf.gz" \
    --threads "${THREADS}" \
    --gtf "${RAMDISK}/annotations.gtf" \
    --cov "${RAMDISK}/coverage.bg.gz" \
    ${SEED:+--seed "${SEED}"}

END_RUN=$(date +%s)
ELAPSED=$((END_RUN - START_RUN))

echo ""
echo "Processing complete in $((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m $((ELAPSED % 60))s"

# ============================================================================
# Index output and verify
# ============================================================================
echo ""
echo "Indexing output VCF (CSI)..."
bcftools index -c "${RAMDISK}/output.vcf.gz"

OUT_COUNT=$(bcftools index -n "${RAMDISK}/output.vcf.gz")
echo "Output variants: ${OUT_COUNT}"

# ============================================================================
# Copy output back to storage
# ============================================================================
echo ""
echo "Copying output to final location..."
cp "${RAMDISK}/output.vcf.gz" "${OUTPUT_VCF}"
cp "${RAMDISK}/output.vcf.gz.csi" "${OUTPUT_VCF}.csi"

# ============================================================================
# Summary statistics
# ============================================================================
echo ""
echo "=============================================="
echo "Summary Statistics"
echo "=============================================="

# Score distributions
echo ""
echo "Score distributions (sample):"
echo "CHROM	POS	DEMUX_SCORE	ANNOT_SCORE	COV_SCORE"
bcftools query -f '%CHROM\t%POS\t%DEMUX_SCORE\t%ANNOT_SCORE\t%COV_SCORE\n' "${OUTPUT_VCF}" 2>/dev/null | head -20 || \
bcftools query -f '%CHROM\t%POS\t%DEMUX_SCORE\n' "${OUTPUT_VCF}" | head -20

# Chromosome distribution
echo ""
echo "Top 10 chromosomes by SNP count:"
bcftools query -f '%CHROM\n' "${OUTPUT_VCF}" | \
    sort | uniq -c | sort -rn | head -10 | \
    awk '{printf "  %s: %d\n", $2, $1}'

# Coverage score stats (if present)
echo ""
echo "Coverage score statistics:"
bcftools query -f '%COV_SCORE\n' "${OUTPUT_VCF}" 2>/dev/null | \
    awk 'BEGIN{min=999999; max=0; sum=0; n=0; zero=0}
         {v=$1+0; if(v<min)min=v; if(v>max)max=v; sum+=v; n++; if(v==0)zero++}
         END{if(n>0) printf "  Min: %.2f\n  Max: %.2f\n  Mean: %.2f\n  Zero coverage: %d (%.1f%%)\n", 
             min, max, sum/n, zero, 100*zero/n}' || echo "  (COV_SCORE not present)"

# Annotation score stats (if present)
echo ""
echo "Annotation score statistics:"
bcftools query -f '%ANNOT_SCORE\n' "${OUTPUT_VCF}" 2>/dev/null | \
    awk 'BEGIN{genic=0; intergenic=0}
         {if($1>=1.0)genic++; else intergenic++}
         END{total=genic+intergenic; 
             if(total>0) printf "  Genic (1.0): %d (%.1f%%)\n  Intergenic (0.1): %d (%.1f%%)\n", 
             genic, 100*genic/total, intergenic, 100*intergenic/total}' || echo "  (ANNOT_SCORE not present)"

echo ""
echo "=============================================="
echo "Pipeline Complete!"
echo "=============================================="
echo "Output: ${OUTPUT_VCF}"
echo "SNPs: ${OUT_COUNT}"
echo "Total runtime: $((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m $((ELAPSED % 60))s"
echo "Finished: $(date)"
echo "=============================================="
