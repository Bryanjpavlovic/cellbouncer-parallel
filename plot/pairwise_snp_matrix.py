#!/usr/bin/env python3
"""
Fast pairwise distinguishing SNP analysis using numpy.
Outputs text report and heatmap visualization.
"""

import sys
import subprocess
import os

try:
    import numpy as np
except ImportError:
    print("ERROR: numpy not installed. Run: pip install numpy", file=sys.stderr)
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
except ImportError:
    print("WARNING: matplotlib not installed. Heatmap will be skipped.", file=sys.stderr)
    plt = None

def get_samples(vcf_file):
    result = subprocess.run(['bcftools', 'query', '-l', vcf_file], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running bcftools: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip().split('\n')

def make_heatmap(pair_matrix, samples, output_file, total_snps):
    """Create a heatmap visualization of pairwise distinguishing SNPs"""
    if plt is None:
        print("Skipping heatmap (matplotlib not available)", flush=True)
        return
    
    n = len(samples)
    
    # Convert to millions for readability
    matrix_millions = pair_matrix / 1_000_000
    
    # Create figure with appropriate size
    fig_size = max(12, n * 0.5)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    
    # Create heatmap with viridis
    im = ax.imshow(matrix_millions, cmap='viridis', aspect='equal')
    
    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.8)
    cbar.ax.set_ylabel('Distinguishing SNPs (millions)', rotation=-90, va="bottom", fontsize=12)
    
    # Set ticks and labels
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(samples, fontsize=9)
    ax.set_yticklabels(samples, fontsize=9)
    
    # Rotate x labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # Add text annotations for each cell
    for i in range(n):
        for j in range(n):
            if i == j:
                text = "-"
                color = "white"
            else:
                val = matrix_millions[i, j]
                text = f"{val:.1f}M"
                # Choose text color based on background
                color = "white" if val < (matrix_millions.max() * 0.5) else "black"
            ax.text(j, i, text, ha="center", va="center", color=color, fontsize=7)
    
    ax.set_title(f'Pairwise Distinguishing SNPs\n({total_snps:,} total SNPs analyzed)', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Heatmap saved to: {output_file}", flush=True)

def main():
    print("Starting pairwise analysis...", flush=True)
    
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <vcf_file> [max_snps] [output_prefix]")
        print(f"  output_prefix: base name for output files (default: pairwise_analysis)")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    max_snps = None
    output_prefix = "pairwise_analysis"
    
    # Parse optional arguments
    for arg in sys.argv[2:]:
        if arg.isdigit():
            max_snps = int(arg)
        else:
            output_prefix = arg
    
    report_file = f"{output_prefix}.txt"
    heatmap_file = f"{output_prefix}.png"
    
    if not os.path.exists(vcf_file):
        print(f"ERROR: File not found: {vcf_file}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Analyzing: {vcf_file}", flush=True)
    print(f"Report will be saved to: {report_file}", flush=True)
    print(f"Heatmap will be saved to: {heatmap_file}", flush=True)
    
    samples = get_samples(vcf_file)
    n = len(samples)
    print(f"Found {n} samples: {', '.join(samples)}", flush=True)
    
    # Initialize pairwise counts matrix
    pair_matrix = np.zeros((n, n), dtype=np.int64)
    
    # Use bcftools to output genotypes
    fmt = '[%GT\\t]\\n'
    cmd = ['bcftools', 'query', '-f', fmt, vcf_file]
    
    print(f"Running: {' '.join(cmd)}", flush=True)
    
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1024*1024)
    
    # Process in batches for speed
    batch_size = 50000
    batch = []
    total_snps = 0
    
    gt_map = {
        '.': -1, './.': -1, '.|.': -1,
        '0/0': 0, '0|0': 0, '0': 0,
        '0/1': 1, '1/0': 1, '0|1': 1, '1|0': 1,
        '1/1': 2, '1|1': 2, '1': 2,
    }
    
    def gt_to_num(gt):
        return gt_map.get(gt, -1)
    
    def process_batch(batch, pair_matrix):
        if not batch:
            return
        
        arr = np.array(batch, dtype=np.int8)
        called = arr >= 0
        
        arr_3d = arr[:, :, np.newaxis]
        arr_3d_t = arr[:, np.newaxis, :]
        
        called_3d = called[:, :, np.newaxis]
        called_3d_t = called[:, np.newaxis, :]
        
        both_called = called_3d & called_3d_t
        different = arr_3d != arr_3d_t
        distinguishing = both_called & different
        
        pair_matrix += np.sum(distinguishing, axis=0)
    
    print("Reading VCF...", flush=True)
    
    for line in proc.stdout:
        fields = line.strip().rstrip('\t').split('\t')
        if len(fields) != n:
            continue
        
        gts = [gt_to_num(g) for g in fields]
        batch.append(gts)
        total_snps += 1
        
        if len(batch) >= batch_size:
            process_batch(batch, pair_matrix)
            batch = []
            print(f"  Processed {total_snps:,} SNPs...", flush=True)
        
        if max_snps and total_snps >= max_snps:
            break
    
    process_batch(batch, pair_matrix)
    proc.wait()
    
    if total_snps == 0:
        print("ERROR: No SNPs processed!", file=sys.stderr)
        sys.exit(1)
    
    # Build output text
    output_lines = []
    output_lines.append(f"=== Pairwise Distinguishing SNP Analysis ===")
    output_lines.append(f"Input: {vcf_file}")
    output_lines.append(f"Total SNPs analyzed: {total_snps:,}")
    output_lines.append(f"Samples: {n}")
    output_lines.append("")
    
    # Pairwise - find lowest and highest
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((samples[i], samples[j], int(pair_matrix[i,j])))
    
    pairs.sort(key=lambda x: x[2])
    
    output_lines.append("=== LOWEST 20 (hardest to distinguish) ===")
    for s1, s2, count in pairs[:20]:
        output_lines.append(f"  {s1} vs {s2}: {count:,} ({100*count/total_snps:.1f}%)")
    
    output_lines.append("")
    output_lines.append("=== HIGHEST 20 (easiest to distinguish) ===")
    for s1, s2, count in pairs[-20:]:
        output_lines.append(f"  {s1} vs {s2}: {count:,} ({100*count/total_snps:.1f}%)")
    
    output_lines.append("")
    output_lines.append("=== ALL PAIRS (sorted by distinguishing SNPs) ===")
    for s1, s2, count in pairs:
        output_lines.append(f"  {s1} vs {s2}: {count:,} ({100*count/total_snps:.1f}%)")
    
    # Matrix - compact version
    output_lines.append("")
    output_lines.append("=== Pairwise matrix (thousands of distinguishing SNPs) ===")
    max_name_len = min(12, max(len(s) for s in samples))
    header = "".ljust(max_name_len+2) + "".join(s[:6].ljust(7) for s in samples)
    output_lines.append(header)
    for i, s1 in enumerate(samples):
        row = s1[:max_name_len].ljust(max_name_len+2)
        for j in range(n):
            if i == j:
                row += "-".ljust(7)
            else:
                k = pair_matrix[i,j] // 1000
                row += f"{k}k".ljust(7)
        output_lines.append(row)
    
    # Print to stdout
    report_text = "\n".join(output_lines)
    print(f"\n{report_text}", flush=True)
    
    # Write to file
    with open(report_file, 'w') as f:
        f.write(report_text)
    print(f"\nReport saved to: {report_file}", flush=True)
    
    # Create heatmap
    make_heatmap(pair_matrix, samples, heatmap_file, total_snps)
    
    print("\nDone!", flush=True)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
