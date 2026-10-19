#!/usr/bin/env python3
"""
Simple QC report for ipyrad assemblies.
Shows basic assembly statistics.
"""
import sys
import os
import re
import gzip
from datetime import datetime
from collections import defaultdict
import numpy as np

# Import common QC functions
sys.path.insert(0, os.path.dirname(__file__))
from qc_common import (
    count_loci_from_file,
    analyze_vcf,
    create_ascii_histogram,
    compute_replicate_discordance_details,
    parse_whoa_summary,
    calculate_pca_variance,
    calculate_genetic_distance_missingness_from_matrix,
)
from extract_sample_stats import extract_vcf_stats
from summarize_hardy import parse_hardy_file, summarize_hwe, build_vcf_position_map
from shared_loci_matrix import parse_loci_presence, shared_loci_matrix
from plot_from_csv import plot_full_assembly_qc
import pandas as pd

# Import ipyrad for JSON loading
try:
    import ipyrad as ip
    HAS_IPYRAD = True
except ImportError:
    HAS_IPYRAD = False

# Redirect output to log
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

# Parse inputs (handle both dict and list formats)
if isinstance(snakemake.input, dict):
    vcf_file = snakemake.input.get('vcf', None)
    loci_file = snakemake.input.get('loci', None)
    hardy_hwe = snakemake.input.get('hardy_hwe', None)
    whoa_report = snakemake.input.get('whoa_report', None)
else:
    # If inputs are a list, try to infer from paths
    vcf_file = None
    loci_file = None
    hardy_hwe = None
    whoa_report = None
    for inp in snakemake.input:
        if inp.endswith('.vcf'):
            vcf_file = inp
        elif inp.endswith('.loci'):
            loci_file = inp
        elif inp.endswith('.hardy.hwe'):
            hardy_hwe = inp
        elif inp.endswith('.whoa_report.txt'):
            whoa_report = inp

hardy_hwe = getattr(snakemake.input, "hardy_hwe", None) or hardy_hwe
whoa_report = getattr(snakemake.input, "whoa_report", None) or whoa_report
vcf_file = getattr(snakemake.input, "vcf", None) or vcf_file
loci_file = getattr(snakemake.input, "loci", None) or loci_file

project_config = snakemake.params.get('project_config', {})
replicate_suffix = (
    snakemake.params.get('replicate_suffix')
    or project_config.get('qc', {}).get('replicate_suffix')
    or project_config.get('optimization', {}).get('replicate_suffix', '_repl')
)
id_pattern = project_config.get('qc', {}).get('id_pattern', r'loc(\d+)_')

assembly_name = snakemake.wildcards.get('assembly_name', 'unknown')
project_dir = snakemake.wildcards.get('project_dir', '')
assembly_dir = snakemake.wildcards.get('assembly_dir', '')

# Construct JSON file path (similar to s7_filter_export.py)
json_path = f"{project_dir}/{assembly_dir}/{assembly_name}/{assembly_name}.json"
# Construct s7 stats file path
s7_stats_path = f"{project_dir}/{assembly_dir}/{assembly_name}/{assembly_name}_outfiles/{assembly_name}_stats.txt"

def count_vcf_loci(vcf_file):
    """
    Count unique loci (RAD fragments) from VCF file.
    Loci are extracted from the ID column (3rd column) by taking unique values
    of the 'loc' part before '_pos'.
    """
    if not vcf_file or not os.path.exists(vcf_file):
        return 0
    
    unique_loci = set()
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    try:
        with opener(vcf_file, 'rt') as f:
            for line in f:
                # Skip comment lines
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                
                # Extract locus ID from the 3rd column (ID field)
                variant_id = fields[2]
                # Extract the 'loc' part before '_pos' (e.g., 'loc0' from 'loc0_pos0')
                match = re.match(r'(loc\d+)_pos', variant_id)
                if match:
                    locus_id = match.group(1)
                    unique_loci.add(locus_id)
                else:
                    # Fallback: try to extract just 'loc' part
                    match = re.match(r'(loc\d+)', variant_id)
                    if match:
                        unique_loci.add(match.group(1))
        
        return len(unique_loci)
    except Exception as e:
        print(f"Warning: Error counting loci from VCF {vcf_file}: {e}", file=sys.stderr)
        return 0

def count_vcf_variants(vcf_file):
    """
    Count variants from VCF file (non-comment lines).
    """
    if not vcf_file or not os.path.exists(vcf_file):
        return 0
    
    count = 0
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    try:
        with opener(vcf_file, 'rt') as f:
            for line in f:
                # Skip comment lines (lines starting with #)
                if not line.startswith('#'):
                    count += 1
        return count
    except Exception as e:
        print(f"Warning: Error counting variants from VCF {vcf_file}: {e}", file=sys.stderr)
        return 0

def extract_positions_within_loci(vcf_file):
    """
    Extract positions within loci from VCF file.
    Returns a list of position values (e.g., [0, 29, 84, ...] from 'loc0_pos0', 'loc0_pos29', etc.)
    """
    if not vcf_file or not os.path.exists(vcf_file):
        return []


def extract_locus_lengths_from_loci(loci_file):
    """
    Estimate locus lengths from ipyrad .loci file.
    Returns:
      - locus_length_map: dict like {"loc0": 145, ...}
      - locus_lengths: list of lengths
    """
    if not loci_file or not os.path.exists(loci_file):
        return {}, []

    locus_length_map = {}
    locus_lengths = []
    current_lengths = []
    locus_idx = 0

    try:
        with open(loci_file, 'r') as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith('//'):
                    if current_lengths:
                        locus_len = int(np.median(current_lengths))
                        locus_id = f"loc{locus_idx}"
                        locus_length_map[locus_id] = locus_len
                        locus_lengths.append(locus_len)
                        locus_idx += 1
                        current_lengths = []
                    continue

                # Typical line: "<sample_name><whitespace><sequence>"
                parts = line.split()
                if len(parts) < 2:
                    continue
                seq = parts[-1].strip()
                if seq:
                    current_lengths.append(len(seq))

        # Handle file without trailing "//"
        if current_lengths:
            locus_len = int(np.median(current_lengths))
            locus_id = f"loc{locus_idx}"
            locus_length_map[locus_id] = locus_len
            locus_lengths.append(locus_len)

        return locus_length_map, locus_lengths

    except Exception as e:
        print(f"Warning: Error extracting locus lengths from {loci_file}: {e}", file=sys.stderr)
        return {}, []


def parse_locus_id_from_variant(variant_id):
    """Extract locus ID token like 'loc123' from variant ID."""
    if not variant_id:
        return None
    match = re.search(r'(loc\d+)', variant_id)
    if match:
        return match.group(1)
    return None


def analyze_variant_structure(vcf_file, locus_length_map=None):
    """
    Analyze SNP composition and position patterns.
    Returns a dict with:
      - allele class counts
      - variants per locus and density
      - normalized SNP positions and enrichment summaries
    """
    summary = {
        'n_variants': 0,
        'n_biallelic': 0,
        'n_triallelic': 0,
        'n_quadriallelic': 0,
        'n_gt4_allelic': 0,
        'variants_per_locus': defaultdict(int),
        'normalized_positions': [],
        'distance_to_edge': [],
        'n_with_norm_pos': 0,
        'start_bin_frac': np.nan,
        'center_bin_frac': np.nan,
        'end_bin_frac': np.nan,
        'center_vs_edges_enrichment': np.nan,
        'edge_vs_center_enrichment': np.nan,
        'snp_density_per_100bp': [],
    }

    if not vcf_file or not os.path.exists(vcf_file):
        return summary

    opener = gzip.open if vcf_file.endswith('.gz') else open

    try:
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue

                summary['n_variants'] += 1
                variant_id = fields[2]
                alt = fields[4]

                # Allelicity
                if alt == '.':
                    n_alleles = 1
                else:
                    n_alleles = 1 + alt.count(',') + 1
                if n_alleles == 2:
                    summary['n_biallelic'] += 1
                elif n_alleles == 3:
                    summary['n_triallelic'] += 1
                elif n_alleles == 4:
                    summary['n_quadriallelic'] += 1
                elif n_alleles > 4:
                    summary['n_gt4_allelic'] += 1

                # Per-locus variant count
                locus_id = parse_locus_id_from_variant(variant_id)
                if locus_id:
                    summary['variants_per_locus'][locus_id] += 1

                # Position within locus from ID token "_posNNN"
                pos_match = re.search(r'_pos(\d+)', variant_id)
                if pos_match and locus_id and locus_length_map:
                    try:
                        pos = int(pos_match.group(1))
                        loc_len = locus_length_map.get(locus_id)
                        if loc_len and loc_len > 1:
                            # Normalize to [0, 1]
                            r = pos / (loc_len - 1)
                            r = max(0.0, min(1.0, r))
                            summary['normalized_positions'].append(r)
                            summary['distance_to_edge'].append(min(r, 1.0 - r))
                            summary['n_with_norm_pos'] += 1
                    except ValueError:
                        pass

        # Density stats when lengths are available
        if summary['variants_per_locus'] and locus_length_map:
            for locus_id, n_snps in summary['variants_per_locus'].items():
                loc_len = locus_length_map.get(locus_id)
                if loc_len and loc_len > 0:
                    density = (n_snps / loc_len) * 100.0
                    summary['snp_density_per_100bp'].append(density)

        # Position enrichment summaries (for concatenation-artifact checks)
        if summary['normalized_positions']:
            pos = np.array(summary['normalized_positions'], dtype=float)
            start_frac = np.mean(pos <= 0.2)
            center_frac = np.mean((pos >= 0.4) & (pos <= 0.6))
            end_frac = np.mean(pos >= 0.8)
            edge_mean = np.mean([start_frac, end_frac])

            summary['start_bin_frac'] = float(start_frac)
            summary['center_bin_frac'] = float(center_frac)
            summary['end_bin_frac'] = float(end_frac)

            if edge_mean > 0:
                summary['center_vs_edges_enrichment'] = float(center_frac / edge_mean)
            if center_frac > 0:
                summary['edge_vs_center_enrichment'] = float(edge_mean / center_frac)

        return summary

    except Exception as e:
        print(f"Warning: Error analyzing variant structure from {vcf_file}: {e}", file=sys.stderr)
        return summary
    
    positions = []
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    try:
        with opener(vcf_file, 'rt') as f:
            for line in f:
                # Skip comment lines
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                
                # Extract position from the ID column (3rd column)
                variant_id = fields[2]
                # Match pattern like 'loc0_pos29' -> extract 29
                match = re.search(r'_pos(\d+)', variant_id)
                if match:
                    pos = int(match.group(1))
                    positions.append(pos)
        
        return positions
    except Exception as e:
        print(f"Warning: Error extracting positions from VCF {vcf_file}: {e}", file=sys.stderr)
        return []

def get_reads_consens_from_json(json_file, replicate_suffix="_repl"):
    """
    Extract reads_consens per sample from assembly JSON.
    Returns a list of tuples: [(sample_name, reads_consens), ...] sorted low to high.
    """
    if not HAS_IPYRAD or not json_file or not os.path.exists(json_file):
        return []

    try:
        data = ip.load_json(json_file)
        reads_consens = []

        if hasattr(data, 'stats') and data.stats is not None and 'reads_consens' in data.stats.columns:
            for sample_name in data.stats.index:
                if replicate_suffix and sample_name.endswith(replicate_suffix):
                    continue
                try:
                    val = data.stats.loc[sample_name, 'reads_consens']
                    if val is None or np.isnan(val):
                        continue
                    reads_consens.append((sample_name, int(float(val))))
                except (ValueError, KeyError, IndexError, TypeError):
                    continue

        reads_consens.sort(key=lambda x: x[1])
        return reads_consens

    except Exception as e:
        print(f"Warning: Error extracting reads_consens from JSON {json_file}: {e}", file=sys.stderr)
        return []


def get_rad_loci_per_sample_from_s7_stats(stats_file, replicate_suffix="_repl"):
    """
    Extract sample_coverage from step7 *_stats.txt.
    Returns a list of tuples: [(sample_name, loci_count), ...] sorted low to high.
    """
    if not stats_file or not os.path.exists(stats_file):
        return []

    try:
        sample_loci = []
        in_samples_section = False
        seen_header = False

        with open(stats_file, 'r') as f:
            for raw_line in f:
                line = raw_line.rstrip('\n')
                stripped = line.strip()

                if "## The number of loci recovered for each Sample." in stripped:
                    in_samples_section = True
                    continue

                if not in_samples_section:
                    continue

                # Keep comment/API lines inside the section until we reach the table header.
                if stripped.startswith("##"):
                    continue

                if not stripped:
                    continue

                if not seen_header:
                    if stripped == "sample_coverage":
                        seen_header = True
                    continue

                # End section at the next titled block.
                if stripped.startswith("The distribution of SNPs"):
                    break

                match = re.match(r'^(\S+)\s+(\d+)$', stripped)
                if not match:
                    continue

                sample_name = match.group(1)
                if replicate_suffix and sample_name.endswith(replicate_suffix):
                    continue

                loci_count = int(match.group(2))
                sample_loci.append((sample_name, loci_count))

        sample_loci.sort(key=lambda x: x[1])
        return sample_loci

    except Exception as e:
        print(f"Warning: Error extracting sample_coverage from stats {stats_file}: {e}", file=sys.stderr)
        return []


def create_sorted_sample_plot(sample_values, title, max_bar_length=60):
    """
    Create an ASCII bar chart showing sample values, sorted from low to high.
    Format: SampleName  ReadCount  ####################
    """
    if not sample_values:
        return []

    lines = []
    lines.append(title)
    lines.append("")

    # Find max value for scaling
    max_value = max(count for _, count in sample_values) if sample_values else 1

    # Find longest sample name for alignment
    max_name_len = max(len(name) for name, _ in sample_values) if sample_values else 20

    for sample_name, value in sample_values:
        # Calculate bar length
        if max_value > 0:
            bar_length = int((value / max_value) * max_bar_length)
        else:
            bar_length = 0

        bar = '#' * bar_length
        # Format: SampleName (padded)  Value (right-aligned)  Bar
        lines.append(f"  {sample_name:<{max_name_len}}  {value:>10,}  {bar}")

    return lines

# Generate report
report_lines = []
report_lines.append("=" * 80)
report_lines.append("QC REPORT")
report_lines.append("=" * 80)
report_lines.append(f"Assembly: {assembly_name}")
report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
report_lines.append("")

# Basic statistics from .loci file (if available)
n_loci_from_file = 0
locus_length_map = {}
locus_lengths = []
if loci_file and os.path.exists(loci_file):
    n_loci_from_file = count_loci_from_file(loci_file)
    report_lines.append(f"Loci (from .loci file): {n_loci_from_file:,}")
    locus_length_map, locus_lengths = extract_locus_lengths_from_loci(loci_file)

# VCF Analysis
vcf_stats = {}
if vcf_file and os.path.exists(vcf_file):
    # Count variants and loci directly from VCF
    n_variants = count_vcf_variants(vcf_file)
    n_loci_vcf = count_vcf_loci(vcf_file)
    
    report_lines.append("")
    report_lines.append("VCF Statistics:")
    report_lines.append(f"  Variants: {n_variants:,}")
    report_lines.append(f"  Loci (RAD fragments): {n_loci_vcf:,}")
    
    # Get detailed VCF statistics
    vcf_stats = analyze_vcf(vcf_file, id_pattern, replicate_suffix)
    
    # Verify variant count matches
    if vcf_stats['total_variants'] != n_variants:
        print(f"Warning: Variant count mismatch: direct count={n_variants}, analyze_vcf={vcf_stats['total_variants']}", file=sys.stderr)
    
    # Verify loci count matches
    n_loci_from_stats = len(vcf_stats['variants_per_locus'])
    if n_loci_from_stats != n_loci_vcf:
        print(f"Warning: Loci count mismatch: direct count={n_loci_vcf}, analyze_vcf={n_loci_from_stats}", file=sys.stderr)
    
    report_lines.append(f"  Samples: {len(vcf_stats['sample_names'])}")
    
    if vcf_stats['non_replicate_samples']:
        report_lines.append(f"  Non-replicate samples: {len(vcf_stats['non_replicate_samples'])}")
    
    if vcf_stats['replicate_pairs']:
        report_lines.append(f"  Replicate pairs: {len(vcf_stats['replicate_pairs'])}")

    # Detailed replicate concordance metrics
    rep = compute_replicate_discordance_details(vcf_file, replicate_suffix)
    rep_summary = rep.get("summary", {})
    rep_pairs = rep.get("pair_rows", [])
    if rep_summary.get("n_pairs", 0) > 0:
        report_lines.append("")
        report_lines.append("Replicate Concordance (technical reproducibility):")
        report_lines.append(f"  Replicate pairs analyzed: {rep_summary['n_pairs']}")
        report_lines.append(f"  Mean discordance: {rep_summary['mean_discordance']:.4f}")
        report_lines.append(f"  SD discordance: {rep_summary['sd_discordance']:.4f}")
        report_lines.append(f"  Mean loci dropout: {rep_summary['mean_loci_dropout']:.4f}")
        report_lines.append(
            f"  Mean allele dropout (hom-het): {rep_summary['mean_allele_dropout_hom_het']:.4f}"
        )
        report_lines.append(
            f"  Mean allele dropout (hom-hom): {rep_summary['mean_allele_dropout_hom_hom']:.4f}"
        )
        report_lines.append(f"  Mean SNP error: {rep_summary['mean_snp_error']:.4f}")

        # Show worst pairs by discordance for targeted QC follow-up
        rep_pairs_sorted = sorted(
            [p for p in rep_pairs if not np.isnan(p.get("replicate_discordance", np.nan))],
            key=lambda x: x["replicate_discordance"],
            reverse=True
        )
        if rep_pairs_sorted:
            report_lines.append("")
            report_lines.append("  Highest-discordance replicate pairs:")
            report_lines.append("    Sample\tReplicate\tDiscordance\tCompared_loci")
            for pair in rep_pairs_sorted[:10]:
                report_lines.append(
                    f"    {pair['sample']}\t{pair['replicate']}\t"
                    f"{pair['replicate_discordance']:.4f}\t{pair['compared_loci']}"
                )
    
    # Missing data statistics
    if vcf_stats['missing_data_per_variant']:
        missing_data = vcf_stats['missing_data_per_variant']
        mean_missing = np.mean(missing_data)
        median_missing = np.median(missing_data)
        report_lines.append("")
        report_lines.append("Missing Data Statistics:")
        report_lines.append(f"  Mean missing data per variant: {mean_missing:.2%}")
        report_lines.append(f"  Median missing data per variant: {median_missing:.2%}")
        
        if vcf_stats['missing_data_per_variant_no_repl']:
            missing_no_repl = vcf_stats['missing_data_per_variant_no_repl']
            mean_missing_no_repl = np.mean(missing_no_repl)
            report_lines.append(f"  Mean missing data per variant (excluding replicates): {mean_missing_no_repl:.2%}")
    
    # Heterozygosity statistics
    if vcf_stats['sample_het_sites'] and vcf_stats['sample_hom_sites']:
        het_rates = []
        for sample in vcf_stats['non_replicate_samples'] if vcf_stats['non_replicate_samples'] else vcf_stats['sample_names']:
            het = vcf_stats['sample_het_sites'].get(sample, 0)
            hom = vcf_stats['sample_hom_sites'].get(sample, 0)
            total = het + hom
            if total > 0:
                het_rates.append(het / total)
        
        if het_rates:
            mean_het = np.mean(het_rates)
            median_het = np.median(het_rates)
            report_lines.append("")
            report_lines.append("Heterozygosity Statistics:")
            report_lines.append(f"  Mean heterozygosity per sample: {mean_het:.4f}")
            report_lines.append(f"  Median heterozygosity per sample: {median_het:.4f}")
    
    # Variants per locus statistics
    if vcf_stats['variants_per_locus']:
        variants_per_locus = [len(variants) for variants in vcf_stats['variants_per_locus'].values()]
        if variants_per_locus:
            mean_variants_per_locus = np.mean(variants_per_locus)
            median_variants_per_locus = np.median(variants_per_locus)
            max_variants_per_locus = np.max(variants_per_locus)
            report_lines.append("")
            report_lines.append("Variants per Locus:")
            report_lines.append(f"  Mean: {mean_variants_per_locus:.2f}")
            report_lines.append(f"  Median: {median_variants_per_locus:.1f}")
            report_lines.append(f"  Max: {max_variants_per_locus}")

    # Locus length distribution from .loci (more interpretable than raw SNP position)
    if locus_lengths:
        report_lines.append("")
        report_lines.append("Locus Length Distribution (.loci):")
        report_lines.append(f"  Loci with inferred lengths: {len(locus_lengths):,}")
        report_lines.append(f"  Mean length: {np.mean(locus_lengths):.1f} bp")
        report_lines.append(f"  Median length: {np.median(locus_lengths):.1f} bp")
        report_lines.append(f"  Min/Max: {np.min(locus_lengths)} / {np.max(locus_lengths)} bp")
        length_hist = create_ascii_histogram(locus_lengths, bins=20, width=60, label="Locus length (bp)")
        for line in length_hist.split('\n'):
            if line.strip():
                report_lines.append(f"  {line}")
            else:
                report_lines.append("")

    # Variant structure and SNP-position diagnostics
    variant_struct = analyze_variant_structure(vcf_file, locus_length_map=locus_length_map)
    if variant_struct['n_variants'] > 0:
        report_lines.append("")
        report_lines.append("Allelicity of Variant Sites:")
        n_total = variant_struct['n_variants']
        n_multi = (
            variant_struct['n_triallelic'] +
            variant_struct['n_quadriallelic'] +
            variant_struct['n_gt4_allelic']
        )
        report_lines.append(
            f"  Biallelic: {variant_struct['n_biallelic']:,} ({variant_struct['n_biallelic']/n_total:.2%})"
        )
        report_lines.append(
            f"  Multiallelic: {n_multi:,} ({n_multi/n_total:.2%}) "
            f"[tri={variant_struct['n_triallelic']:,}, quad={variant_struct['n_quadriallelic']:,}, >4={variant_struct['n_gt4_allelic']:,}]"
        )

        if variant_struct['variants_per_locus']:
            vpl = list(variant_struct['variants_per_locus'].values())
            report_lines.append("")
            report_lines.append("SNPs per Locus:")
            report_lines.append(f"  Mean: {np.mean(vpl):.2f}")
            report_lines.append(f"  Median: {np.median(vpl):.1f}")
            report_lines.append(f"  95th percentile: {np.percentile(vpl, 95):.1f}")
            report_lines.append(f"  Max: {np.max(vpl)}")
            vpl_hist = create_ascii_histogram(vpl, bins=20, width=60, label="SNPs per locus")
            for line in vpl_hist.split('\n'):
                if line.strip():
                    report_lines.append(f"  {line}")
                else:
                    report_lines.append("")

        if variant_struct['snp_density_per_100bp']:
            dens = variant_struct['snp_density_per_100bp']
            report_lines.append("")
            report_lines.append("SNP Density per Locus Length:")
            report_lines.append(f"  Mean: {np.mean(dens):.2f} SNPs/100bp")
            report_lines.append(f"  Median: {np.median(dens):.2f} SNPs/100bp")
            report_lines.append(f"  95th percentile: {np.percentile(dens, 95):.2f} SNPs/100bp")
            dens_hist = create_ascii_histogram(dens, bins=20, width=60, label="SNP density (SNPs/100bp)")
            for line in dens_hist.split('\n'):
                if line.strip():
                    report_lines.append(f"  {line}")
                else:
                    report_lines.append("")

        # Normalized SNP position checks for edge/center enrichment artifacts
        if variant_struct['normalized_positions']:
            report_lines.append("")
            report_lines.append("Normalized SNP Position Pattern (artifact check):")
            report_lines.append(
                f"  SNPs with valid normalized positions: {variant_struct['n_with_norm_pos']:,}"
            )
            report_lines.append(
                f"  Start bin (0.00-0.20): {variant_struct['start_bin_frac']:.2%}"
            )
            report_lines.append(
                f"  Center bin (0.40-0.60): {variant_struct['center_bin_frac']:.2%}"
            )
            report_lines.append(
                f"  End bin (0.80-1.00): {variant_struct['end_bin_frac']:.2%}"
            )

            cve = variant_struct['center_vs_edges_enrichment']
            evc = variant_struct['edge_vs_center_enrichment']
            if not np.isnan(cve):
                report_lines.append(f"  Center vs mean(edges) enrichment: {cve:.2f}x")
            if not np.isnan(evc):
                report_lines.append(f"  Mean(edges) vs center enrichment: {evc:.2f}x")

            norm_hist = create_ascii_histogram(
                variant_struct['normalized_positions'],
                bins=20,
                width=60,
                label="Normalized position in locus (0=start, 1=end)"
            )
            for line in norm_hist.split('\n'):
                if line.strip():
                    report_lines.append(f"  {line}")
                else:
                    report_lines.append("")

            edge_hist = create_ascii_histogram(
                variant_struct['distance_to_edge'],
                bins=20,
                width=60,
                label="Distance to nearest edge (normalized)"
            )
            for line in edge_hist.split('\n'):
                if line.strip():
                    report_lines.append(f"  {line}")
                else:
                    report_lines.append("")

            # Simple interpretation flags
            report_lines.append("")
            report_lines.append("Position Pattern Flags:")
            if not np.isnan(cve) and cve >= 1.5:
                report_lines.append("  WARNING: Center-enriched SNP pattern (possible read-end concatenation artifact).")
            elif not np.isnan(evc) and evc >= 1.5:
                report_lines.append("  WARNING: Edge-enriched SNP pattern (possible end-effect or trimming artifact).")
            else:
                report_lines.append("  OK: No strong center/edge enrichment signal detected.")

# Reads in consensus per sample (from JSON)
sample_reads_consens = []
if os.path.exists(json_path):
    sample_reads_consens = get_reads_consens_from_json(json_path, replicate_suffix)

if sample_reads_consens:
    report_lines.append("")
    plot_lines = create_sorted_sample_plot(
        sample_reads_consens,
        "Reads in Consensus per Sample (sorted low to high):",
        max_bar_length=60,
    )
    for line in plot_lines:
        report_lines.append(line)

# RAD loci per sample (from step7 stats file)
sample_loci = []
if os.path.exists(s7_stats_path):
    sample_loci = get_rad_loci_per_sample_from_s7_stats(s7_stats_path, replicate_suffix)

if sample_loci:
    report_lines.append("")
    plot_lines = create_sorted_sample_plot(
        sample_loci,
        "RAD Loci per Sample (sorted low to high):",
        max_bar_length=60,
    )
    for line in plot_lines:
        report_lines.append(line)

report_lines.append("")
report_lines.append("=" * 80)

# Write report
output_report = snakemake.output.report if hasattr(snakemake.output, "report") else snakemake.output[0]
os.makedirs(os.path.dirname(output_report) or ".", exist_ok=True)
with open(output_report, 'w') as f:
    f.write('\n'.join(report_lines))

print(f"QC report written to {output_report}")

# Paper supplementary figure + metrics CSV
output_csv = snakemake.output.get("csv") if hasattr(snakemake.output, "get") else None
output_pdf = snakemake.output.get("pdf") if hasattr(snakemake.output, "get") else None
if not output_csv:
    output_csv = getattr(snakemake.output, "csv", None)
if not output_pdf:
    output_pdf = getattr(snakemake.output, "pdf", None)

if output_csv or output_pdf:
    metric_rows = []
    sample_rows = []
    extracted = {}
    if vcf_file and os.path.exists(vcf_file):
        extracted = extract_vcf_stats(vcf_file, json_path, replicate_suffix)
        for sname, svals in extracted.get("per_sample", {}).items():
            row = {"sample": sname}
            row.update(svals)
            sample_rows.append(row)
            for key, val in svals.items():
                metric_rows.append({"metric": key, "sample": sname, "value": val})
        for metric_name, values in [
            ("vcf_dp", extracted.get("dp_values", [])),
            ("vcf_maf", extracted.get("maf_values", [])),
            ("vcf_ns", extracted.get("ns_values", [])),
        ]:
            for val in values:
                metric_rows.append({"metric": metric_name, "sample": "", "value": val})
        for key in [
            "n_fragments", "n_variants", "n_biallelic", "n_triallelic",
            "n_quadriallelic", "n_gt4_allelic", "singleton_proportion",
            "doubleton_proportion",
        ]:
            metric_rows.append({"metric": f"vcf_{key}", "sample": "", "value": extracted.get(key, np.nan)})
        for key in ("hdplot_n_tested", "hdplot_n_flagged", "hdplot_prop_flagged"):
            metric_rows.append({"metric": key, "sample": "", "value": extracted.get(key, np.nan)})
        dist_miss = None
        sample_names = extracted.get("sample_names") or []
        dosage = extracted.get("dosage_matrix")
        per_sample = extracted.get("per_sample") or {}
        if sample_names and isinstance(dosage, np.ndarray) and dosage.size:
            missing_rates = [per_sample.get(s, {}).get("missing_rate", np.nan) for s in sample_names]
            dist_miss = calculate_genetic_distance_missingness_from_matrix(dosage, missing_rates)
            if dist_miss is not None:
                metric_rows.append({
                    "metric": "genetic_distance_missingness_corr",
                    "sample": "",
                    "value": dist_miss,
                })

    hwe_summary = {}
    if hardy_hwe and os.path.exists(hardy_hwe):
        try:
            position_map = build_vcf_position_map(vcf_file) if vcf_file else None
            hwe_loci = parse_hardy_file(hardy_hwe, position_map=position_map)
            hwe_summary = summarize_hwe(hwe_loci, alpha=0.05)
            valid = max(1, hwe_summary.get("valid_loci", 0))
            metric_rows.append({"metric": "hwe_prop_rejected", "sample": "", "value": hwe_summary.get("sig_deviations", 0) / valid if hwe_summary.get("valid_loci", 0) else np.nan})
            metric_rows.append({"metric": "hwe_prop_deficit", "sample": "", "value": hwe_summary.get("sig_deficit", 0) / valid if hwe_summary.get("valid_loci", 0) else np.nan})
            metric_rows.append({"metric": "hwe_prop_excess", "sample": "", "value": hwe_summary.get("sig_excess", 0) / valid if hwe_summary.get("valid_loci", 0) else np.nan})
            for pval in hwe_summary.get("p_values", []):
                metric_rows.append({"metric": "hwe_p", "sample": "", "value": pval})
        except Exception as e:
            print(f"Warning: HWE summary for figure failed: {e}", file=sys.stderr)
            hwe_summary = {}

    whoa_summary = parse_whoa_summary(whoa_report) if whoa_report else None
    if whoa_summary:
        metric_rows.append({
            "metric": "whoa_mean_miscall_rate",
            "sample": "",
            "value": whoa_summary.get("mean_miscall_rate", np.nan),
        })
        for bin_row in whoa_summary.get("depth_bins", []):
            metric_rows.append({
                "metric": "whoa_mean_m",
                "sample": "",
                "value": bin_row.get("mean_m", np.nan),
                "mean_dp": bin_row.get("mean_dp", np.nan),
            })

    rep = compute_replicate_discordance_details(vcf_file, replicate_suffix) if vcf_file else {"pair_rows": [], "summary": {}}
    replicate_df = pd.DataFrame(rep.get("pair_rows") or [])
    for pair in rep.get("pair_rows") or []:
        metric_rows.append({
            "metric": "replicate_discordance",
            "sample": pair.get("sample", ""),
            "value": pair.get("replicate_discordance", np.nan),
        })

    pca_details = None
    if vcf_file and os.path.exists(vcf_file):
        pca_details = calculate_pca_variance(
            vcf_file, n_components=8, exclude_replicates=True,
            replicate_suffix=replicate_suffix, return_details=True
        )
        if isinstance(pca_details, dict):
            metric_rows.append({
                "metric": "pca_cumulative_variance",
                "sample": "",
                "value": pca_details.get("cumulative_variance", np.nan),
            })

    snps_per_locus = []
    missing_per_locus = []
    if vcf_file and os.path.exists(vcf_file):
        if vcf_stats.get("variants_per_locus"):
            snps_per_locus = [len(v) for v in vcf_stats["variants_per_locus"].values()]
        for locus_id, total in vcf_stats.get("total_calls_per_locus", {}).items():
            miss = vcf_stats.get("missing_data_per_locus", {}).get(locus_id, 0)
            if total:
                missing_per_locus.append(miss / total)
        for val in snps_per_locus:
            metric_rows.append({"metric": "snps_per_locus", "sample": "", "value": val})
        for val in missing_per_locus:
            metric_rows.append({"metric": "missing_per_locus", "sample": "", "value": val})

    for val in locus_lengths:
        metric_rows.append({"metric": "locus_length", "sample": "", "value": val})

    shared_samples = None
    shared_matrix = None
    if loci_file and os.path.exists(loci_file):
        try:
            shared_samples, presence = parse_loci_presence(loci_file)
            shared_matrix = shared_loci_matrix(presence)
        except Exception as e:
            print(f"Warning: shared-loci parse for figure failed: {e}", file=sys.stderr)

    sample_df = pd.DataFrame(sample_rows)
    if output_csv:
        os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
        pd.DataFrame(metric_rows).to_csv(output_csv, index=False)
        print(f"QC metrics CSV written to {output_csv}")

    if output_pdf:
        os.makedirs(os.path.dirname(output_pdf) or ".", exist_ok=True)
        plot_full_assembly_qc(output_pdf, {
            "sample_df": sample_df,
            "replicate_df": replicate_df,
            "pca": pca_details if isinstance(pca_details, dict) else None,
            "whoa": whoa_summary,
            "hwe": hwe_summary,
            "dp_values": extracted.get("dp_values", []),
            "maf_values": extracted.get("maf_values", []),
            "snps_per_locus": snps_per_locus,
            "locus_lengths": locus_lengths,
            "missing_per_locus": missing_per_locus,
            "shared_samples": shared_samples,
            "shared_matrix": shared_matrix,
            "hdplot_h": extracted.get("hdplot_h", []),
            "hdplot_d": extracted.get("hdplot_d", []),
            "hdplot_n_tested": extracted.get("hdplot_n_tested", 0),
            "hdplot_n_flagged": extracted.get("hdplot_n_flagged", 0),
        })


