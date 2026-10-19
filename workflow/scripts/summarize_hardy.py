#!/usr/bin/env python3
"""
Summarize Hardy-Weinberg equilibrium test results from vcftools --hardy output.
"""
import argparse
import sys
import os
import numpy as np
import gzip
import re
from collections import defaultdict, Counter

# Import common functions
sys.path.insert(0, os.path.dirname(__file__))
from qc_common import create_ascii_histogram, extract_rad_locus

def build_vcf_position_map(vcf_file):
    """Build a mapping from (CHR, POS) to variant ID and RAD locus."""
    position_map = {}
    opener = gzip.open if vcf_file.endswith('.gz') else open

    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = fields[1]
            variant_id = fields[2]

            # Extract RAD locus from variant ID
            rad_locus = extract_rad_locus(variant_id)

            position_map[(chrom, pos)] = {
                'id': variant_id,
                'rad_locus': rad_locus
            }

    return position_map

def parse_hardy_file(hardy_file, position_map=None):
    """Parse vcftools .hwe file."""
    loci = []

    with open(hardy_file, 'r') as f:
        header = f.readline().strip()

        for line in f:
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 6:
                    chrom = parts[0]
                    pos = parts[1]
                    obs_counts = parts[2]
                    exp_counts = parts[3]
                    chisq = float(parts[4]) if parts[4] != 'nan' else np.nan
                    p_hwe = float(parts[5]) if parts[5] != 'nan' else np.nan
                    p_deficit = float(parts[6]) if len(parts) > 6 and parts[6] != 'nan' else np.nan
                    p_excess = float(parts[7]) if len(parts) > 7 and parts[7] != 'nan' else np.nan

                    # Get variant ID and RAD locus from position map if available
                    variant_id = None
                    rad_locus = None
                    if position_map and (chrom, pos) in position_map:
                        variant_id = position_map[(chrom, pos)]['id']
                        rad_locus = position_map[(chrom, pos)]['rad_locus']

                    loci.append({
                        'chrom': chrom,
                        'pos': pos,
                        'variant_id': variant_id,
                        'rad_locus': rad_locus,
                        'obs_counts': obs_counts,
                        'exp_counts': exp_counts,
                        'chisq': chisq,
                        'p_hwe': p_hwe,
                        'p_deficit': p_deficit,
                        'p_excess': p_excess
                    })

    return loci

def summarize_hwe(loci, alpha=0.05):
    """Generate summary statistics for HWE test results."""
    valid_loci = [l for l in loci if not np.isnan(l['p_hwe'])]

    total_loci = len(loci)
    valid_loci_count = len(valid_loci)

    if valid_loci_count == 0:
        return {
            'total_loci': total_loci,
            'valid_loci': 0,
            'sig_deviations': 0,
            'sig_deficit': 0,
            'sig_excess': 0,
            'mean_p': np.nan,
            'median_p': np.nan,
            'p_values': [],
            'top_deviations': [],
            'rad_locus_stats': {}
        }

    sig_deviations = sum(1 for l in valid_loci if l['p_hwe'] < alpha)
    sig_deficit = sum(1 for l in valid_loci if not np.isnan(l['p_deficit']) and l['p_deficit'] < alpha)
    sig_excess = sum(1 for l in valid_loci if not np.isnan(l['p_excess']) and l['p_excess'] < alpha)

    p_values = [l['p_hwe'] for l in valid_loci]
    mean_p = np.mean(p_values)
    median_p = np.median(p_values)

    sorted_loci = sorted(valid_loci, key=lambda x: x['p_hwe'])
    top_deviations = sorted_loci[:100]

    # Analyze RAD locus statistics
    rad_locus_stats = {}
    loci_with_rad = [l for l in valid_loci if l.get('rad_locus') is not None]

    if loci_with_rad:
        rad_locus_total = Counter(l['rad_locus'] for l in loci_with_rad)
        rad_locus_sig = Counter(l['rad_locus'] for l in loci_with_rad if l['p_hwe'] < alpha)

        for rad_locus, total_snps in rad_locus_total.items():
            sig_snps = rad_locus_sig.get(rad_locus, 0)
            rad_locus_stats[rad_locus] = {
                'total_snps': total_snps,
                'sig_snps': sig_snps,
                'prop_sig': sig_snps / total_snps if total_snps > 0 else 0
            }

    return {
        'total_loci': total_loci,
        'valid_loci': valid_loci_count,
        'sig_deviations': sig_deviations,
        'sig_deficit': sig_deficit,
        'sig_excess': sig_excess,
        'mean_p': mean_p,
        'median_p': median_p,
        'p_values': p_values,
        'top_deviations': top_deviations,
        'rad_locus_stats': rad_locus_stats
    }

def write_report(summary, output_file, alpha=0.05, vcf_file=None):
    """Write summary report."""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("HARDY-WEINBERG EQUILIBRIUM SUMMARY\n")
        f.write("=" * 80 + "\n\n")

        if vcf_file:
            f.write(f"Input VCF: {vcf_file}\n")
        f.write(f"Significance threshold (alpha): {alpha}\n\n")

        f.write("-" * 80 + "\n")
        f.write("STATISTICS\n")
        f.write("-" * 80 + "\n\n")

        f.write(f"Total loci tested: {summary['total_loci']}\n")
        f.write(f"Valid loci (with p-values): {summary['valid_loci']}\n")
        if summary['valid_loci'] > 0:
            f.write(f"HWE deviations (p < {alpha}): {summary['sig_deviations']} ({100.0 * summary['sig_deviations'] / summary['valid_loci']:.2f}%)\n")
            f.write(f"Heterozygote deficit (p < {alpha}): {summary['sig_deficit']} ({100.0 * summary['sig_deficit'] / summary['valid_loci']:.2f}%)\n")
            f.write(f"Heterozygote excess (p < {alpha}): {summary['sig_excess']} ({100.0 * summary['sig_excess'] / summary['valid_loci']:.2f}%)\n\n")
        f.write(f"Mean p-value: {summary['mean_p']:.6f}\n")
        f.write(f"Median p-value: {summary['median_p']:.6f}\n\n")

        # P-value distribution histogram
        if summary['p_values']:
            f.write("-" * 80 + "\n")
            f.write("P-VALUE DISTRIBUTION\n")
            f.write("-" * 80 + "\n\n")
            f.write(create_ascii_histogram(summary['p_values'], bins=20, label="P-values"))

        f.write("-" * 80 + "\n")
        f.write("TOP 100 MOST SIGNIFICANT DEVIATIONS\n")
        f.write("-" * 80 + "\n\n")

        if summary['top_deviations']:
            has_variant_ids = summary['top_deviations'][0].get('variant_id') is not None

            if has_variant_ids:
                f.write(f"{'Variant ID':<20} {'RAD Locus':<15} {'CHR':<10} {'POS':<10} {'P_HWE':<12} {'Observed':<20} {'Expected':<20}\n")
                f.write("-" * 120 + "\n")
                for locus in summary['top_deviations']:
                    variant_id = locus.get('variant_id', 'N/A')
                    rad_locus = locus.get('rad_locus', 'N/A')
                    f.write(f"{variant_id:<20} {rad_locus:<15} {locus['chrom']:<10} {locus['pos']:<10} {locus['p_hwe']:<12.6e} ")
                    f.write(f"{locus['obs_counts']:<20} {locus['exp_counts']:<20}\n")
            else:
                f.write(f"{'CHR':<15} {'POS':<10} {'P_HWE':<12} {'Observed':<20} {'Expected':<20}\n")
                f.write("-" * 80 + "\n")
                for locus in summary['top_deviations']:
                    f.write(f"{locus['chrom']:<15} {locus['pos']:<10} {locus['p_hwe']:<12.6e} ")
                    f.write(f"{locus['obs_counts']:<20} {locus['exp_counts']:<20}\n")
        else:
            f.write("No significant deviations found.\n")

        f.write("\n")

        # RAD locus analysis
        if summary['rad_locus_stats']:
            f.write("-" * 80 + "\n")
            f.write("RAD LOCUS ANALYSIS\n")
            f.write("-" * 80 + "\n\n")

            total_rad_loci = len(summary['rad_locus_stats'])
            rad_loci_with_sig = sum(1 for stats in summary['rad_locus_stats'].values() if stats['sig_snps'] > 0)
            rad_loci_all_sig = sum(1 for stats in summary['rad_locus_stats'].values()
                                   if stats['sig_snps'] == stats['total_snps'] and stats['total_snps'] > 0)

            f.write(f"Total RAD loci: {total_rad_loci}\n")
            f.write(f"Loci with ≥1 significant SNP: {rad_loci_with_sig} ({100.0 * rad_loci_with_sig / total_rad_loci:.1f}%)\n")
            f.write(f"Loci with all SNPs significant: {rad_loci_all_sig} ({100.0 * rad_loci_all_sig / total_rad_loci:.1f}%)\n\n")

            # Loci where 100% of SNPs deviate
            loci_100_percent = sorted(
                [(locus, stats) for locus, stats in summary['rad_locus_stats'].items()
                 if stats['sig_snps'] == stats['total_snps'] and stats['total_snps'] > 0],
                key=lambda x: x[0]
            )

            if loci_100_percent:
                f.write("-" * 80 + "\n")
                f.write(f"RAD loci with all SNPs deviating ({len(loci_100_percent)} loci):\n")
                f.write("-" * 80 + "\n\n")
                f.write(f"{'RAD Locus':<15} {'Total SNPs':<12} {'Sig SNPs':<12} {'Proportion':<12}\n")
                f.write("-" * 80 + "\n")

                for rad_locus, stats in loci_100_percent:
                    f.write(f"{rad_locus:<15} {stats['total_snps']:<12} {stats['sig_snps']:<12} {stats['prop_sig']:<12.3f}\n")

                f.write("\n")

            # Distribution of deviation proportions per RAD locus
            prop_sig_values = [stats['prop_sig'] for stats in summary['rad_locus_stats'].values()]
            if prop_sig_values:
                f.write("-" * 80 + "\n")
                f.write("DEVIATION RATE DISTRIBUTION (per RAD locus)\n")
                f.write("-" * 80 + "\n\n")
                f.write(create_ascii_histogram(prop_sig_values, bins=20, label="Proportion of significant SNPs per RAD locus"))

        f.write("=" * 80 + "\n")

def main():
    parser = argparse.ArgumentParser(
        description='Summarize Hardy-Weinberg equilibrium test results'
    )
    parser.add_argument('--input', required=True,
                        help='Input .hwe file from vcftools --hardy')
    parser.add_argument('--vcf', required=False,
                        help='Optional VCF file to extract variant IDs and RAD locus information')
    parser.add_argument('--output', required=True,
                        help='Output summary report file')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='Significance threshold (default: 0.05)')

    args = parser.parse_args()

    position_map = None
    if args.vcf:
        print(f"Reading VCF file: {args.vcf}", file=sys.stderr)
        position_map = build_vcf_position_map(args.vcf)

    print(f"Reading Hardy-Weinberg results: {args.input}", file=sys.stderr)
    loci = parse_hardy_file(args.input, position_map=position_map)

    if not loci:
        print("ERROR: No loci found in input file", file=sys.stderr)
        sys.exit(1)

    print(f"Summarizing {len(loci)} loci...", file=sys.stderr)
    summary = summarize_hwe(loci, alpha=args.alpha)

    print(f"Writing report: {args.output}", file=sys.stderr)
    write_report(summary, args.output, alpha=args.alpha, vcf_file=args.vcf)

    print(f"Total loci: {summary['total_loci']}", file=sys.stderr)
    if summary['valid_loci'] > 0:
        print(f"Significant deviations: {summary['sig_deviations']} ({100.0 * summary['sig_deviations'] / summary['valid_loci']:.1f}%)", file=sys.stderr)

if __name__ == '__main__':
    main()
