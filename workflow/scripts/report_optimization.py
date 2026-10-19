#!/usr/bin/env python3
"""
Optimization threshold comparison.
Extracts per-sample statistics (excluding replicates) and generates plots.
Runs only in Snakemake mode.
"""
import sys
import os
import subprocess
from itertools import combinations

# Snakemake injects the 'snakemake' object when using script: directive
# Check if snakemake is available - if not, this script is being run incorrectly
if 'snakemake' not in globals():
    raise RuntimeError(
        "This script must be run within Snakemake using the 'script:' directive. "
        "The 'snakemake' object was not found. Make sure you're using 'script:' not 'shell:' in the Snakefile."
    )

# Redirect output to log
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

# Import the extraction functions
sys.path.insert(0, os.path.dirname(__file__))
from extract_sample_stats import (
    extract_sample_stats,
    extract_vcf_stats,
    extract_hwe_per_sample
)
from qc_common import (
    calculate_genetic_distance_missingness_from_matrix,
    compute_replicate_discordance_details,
    parse_whoa_summary,
)
from summarize_hardy import parse_hardy_file, summarize_hwe, build_vcf_position_map
from shared_loci_matrix import parse_loci_presence, shared_loci_matrix, mean_pairwise_shared
import ipyrad as ip
import pandas as pd
import numpy as np

# Parse parameters from snakemake
def _as_list(value):
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        return [value]
    try:
        return list(value)
    except TypeError:
        return [value]

thresholds = snakemake.params.thresholds
project_dir = snakemake.wildcards.project_dir
json_files = _as_list(snakemake.input.jsons)
vcf_files = _as_list(snakemake.input.get('vcfs', []))
hardy_files = _as_list(snakemake.input.get('hardy_reports', []))
loci_files = _as_list(snakemake.input.get('loci', []))
whoa_files = _as_list(snakemake.input.get('whoa_reports', []))
output_csv = snakemake.output.csv
output_vcf_csv = snakemake.output.get('vcf_csv', None)
output_hwe_csv = snakemake.output.get('hwe_csv', None)
output_replicate_csv = snakemake.output.get('replicate_csv', None)
output_pdf = snakemake.output.pdf
output_supp_pdf = snakemake.output.get('supp_pdf', None)
replicate_suffix = snakemake.params.get("replicate_suffix", "_repl")
script_dir = os.path.dirname(__file__)

# Step 1: Extract sample statistics
print("Extracting per-sample statistics (excluding replicates)...")

all_stats = []
vcf_stats_list = []  # Per-locus VCF statistics
hwe_stats = []  # Per-sample HWE deviation counts
replicate_stats = []  # Per-pair replicate stats across thresholds
matrix_by_threshold = {}
variant_ids_by_threshold = {}

def _append_metric(threshold, assembly_name, metric, value, **extra):
    row = {
        'threshold': threshold,
        'assembly_name': assembly_name,
        'metric': metric,
        'value': value,
    }
    row.update(extra)
    vcf_stats_list.append(row)

def _resolve_sidecar(explicit_list, index, json_file, *relative_parts):
    if index < len(explicit_list) and explicit_list[index]:
        candidate = explicit_list[index]
        if os.path.exists(candidate):
            return candidate
    if json_file:
        json_dir = os.path.dirname(json_file)
        assembly_name = os.path.basename(json_file).replace('.json', '')
        inferred = os.path.join(json_dir, *relative_parts).format(assembly_name=assembly_name)
        if os.path.exists(inferred):
            return inferred
    return None

for i, threshold in enumerate(thresholds):
    json_file = json_files[i] if i < len(json_files) else None
    
    # Get VCF file path
    vcf_file = _resolve_sidecar(
        vcf_files, i, json_file, "{assembly_name}_outfiles", "{assembly_name}.vcf"
    )

    # Get HWE file path
    hwe_file = _resolve_sidecar(
        hardy_files, i, json_file, "{assembly_name}_outfiles", "qc", "{assembly_name}.hardy.hwe"
    )

    loci_file = _resolve_sidecar(
        loci_files, i, json_file, "{assembly_name}_outfiles", "{assembly_name}.loci"
    )
    whoa_file = _resolve_sidecar(
        whoa_files, i, json_file, "{assembly_name}_outfiles", "qc", "{assembly_name}.whoa_report.txt"
    )
    
    # Extract sample statistics from JSON
    threshold_sample_rows = []
    assembly_name = f'optimization_{threshold}'
    if json_file and os.path.exists(json_file):
        stats = extract_sample_stats(json_file, replicate_suffix)
        # Add threshold to each stat
        for stat in stats:
            stat['threshold'] = threshold
        threshold_sample_rows = stats
        try:
            data = ip.load_json(json_file)
            assembly_name = data.name if hasattr(data, 'name') else assembly_name
        except Exception:
            assembly_name = os.path.basename(json_file).replace('.json', '')

    # Extract VCF statistics (per-locus)
    if vcf_file and os.path.exists(vcf_file):
        vcf_stats = extract_vcf_stats(vcf_file, json_file, replicate_suffix)
        # Create per-locus rows for NS, DP, and MAF
        for ns in vcf_stats['ns_values']:
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'vcf_ns',
                'value': ns
            })
        
        for dp in vcf_stats['dp_values']:
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'vcf_dp',
                'value': dp
            })

        for maf in vcf_stats.get('maf_values', []):
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'vcf_maf',
                'value': maf
            })
        
        # Add summary statistics
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'vcf_n_fragments',
            'value': vcf_stats['n_fragments']
        })
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'vcf_n_variants',
            'value': vcf_stats['n_variants']
        })
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'vcf_n_biallelic',
            'value': vcf_stats['n_biallelic']
        })
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'vcf_n_triallelic',
            'value': vcf_stats['n_triallelic']
        })
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'vcf_n_quadriallelic',
            'value': vcf_stats['n_quadriallelic']
        })
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'vcf_n_gt4_allelic',
            'value': vcf_stats['n_gt4_allelic']
        })
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'vcf_singleton_proportion',
            'value': vcf_stats.get('singleton_proportion', np.nan)
        })
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'vcf_doubleton_proportion',
            'value': vcf_stats.get('doubleton_proportion', np.nan)
        })

        for cutoff, retained in vcf_stats.get('missingness_retention', {}).items():
            pct = int(round(float(cutoff) * 100))
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': f'loci_retained_ge{pct}pct',
                'value': retained
            })

        _append_metric(
            threshold, assembly_name, 'hdplot_n_tested',
            vcf_stats.get('hdplot_n_tested', 0)
        )
        _append_metric(
            threshold, assembly_name, 'hdplot_n_flagged',
            vcf_stats.get('hdplot_n_flagged', 0)
        )
        _append_metric(
            threshold, assembly_name, 'hdplot_prop_flagged',
            vcf_stats.get('hdplot_prop_flagged', np.nan)
        )

        # Merge per-sample VCF-derived metrics into sample rows
        if threshold_sample_rows:
            per_sample = vcf_stats.get('per_sample', {})
            for row in threshold_sample_rows:
                sname = row.get('sample')
                svals = per_sample.get(sname, {})
                row['vcf_call_rate'] = svals.get('call_rate', np.nan)
                row['vcf_missing_rate'] = svals.get('missing_rate', np.nan)
                row['vcf_observed_het'] = svals.get('observed_het', np.nan)
                row['vcf_expected_het'] = svals.get('expected_het', np.nan)
                row['vcf_fis'] = svals.get('fis', np.nan)
                # Legacy/compat metric alias requested in interpretation notes
                row['singleton_proportion'] = vcf_stats.get('singleton_proportion', np.nan)

        # Collect dosage matrix for cross-threshold stability analyses
        sample_names = vcf_stats.get('sample_names', [])
        dosage_matrix = vcf_stats.get('dosage_matrix', np.empty((0, 0), dtype=float))
        if isinstance(dosage_matrix, np.ndarray) and dosage_matrix.size > 0 and sample_names:
            matrix_by_threshold[threshold] = {
                'assembly_name': assembly_name,
                'sample_names': sample_names,
                'dosage_matrix': dosage_matrix,
            }
            missing_rates = [
                vcf_stats.get('per_sample', {}).get(sname, {}).get('missing_rate', np.nan)
                for sname in sample_names
            ]
            dist_miss = calculate_genetic_distance_missingness_from_matrix(
                dosage_matrix, missing_rates
            )
            if dist_miss is not None:
                _append_metric(
                    threshold, assembly_name,
                    'genetic_distance_missingness_corr',
                    dist_miss,
                )
        variant_ids_by_threshold[threshold] = set(vcf_stats.get('variant_ids') or [])

        # Detailed replicate concordance / technical reproducibility
        rep = compute_replicate_discordance_details(vcf_file, replicate_suffix)
        rep_summary = rep.get("summary", {})
        rep_pairs = rep.get("pair_rows", [])
        if rep_pairs:
            for row in rep_pairs:
                replicate_stats.append({
                    "threshold": threshold,
                    "assembly_name": assembly_name,
                    **row,
                })
        if rep_summary.get("n_pairs", 0) > 0:
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'replicate_discordance_mean',
                'value': rep_summary.get("mean_discordance", np.nan)
            })
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'replicate_discordance_sd',
                'value': rep_summary.get("sd_discordance", np.nan)
            })
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'replicate_loci_dropout_mean',
                'value': rep_summary.get("mean_loci_dropout", np.nan)
            })
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'replicate_allele_dropout_hom_het_mean',
                'value': rep_summary.get("mean_allele_dropout_hom_het", np.nan)
            })
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'replicate_allele_dropout_hom_hom_mean',
                'value': rep_summary.get("mean_allele_dropout_hom_hom", np.nan)
            })
            vcf_stats_list.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'metric': 'replicate_snp_error_mean',
                'value': rep_summary.get("mean_snp_error", np.nan)
            })

    if loci_file and os.path.exists(loci_file):
        try:
            _samples, presence = parse_loci_presence(loci_file)
            shared = shared_loci_matrix(presence)
            _append_metric(
                threshold, assembly_name, 'mean_pairwise_shared_loci',
                mean_pairwise_shared(shared)
            )
            if shared.size:
                totals = np.diag(shared)
                _append_metric(
                    threshold, assembly_name, 'mean_loci_per_sample',
                    float(np.mean(totals))
                )
        except Exception as e:
            print(f"Warning: shared-loci summary failed at threshold {threshold}: {e}", file=sys.stderr)

    if whoa_file and os.path.exists(whoa_file):
        whoa = parse_whoa_summary(whoa_file)
        if whoa:
            if whoa.get('mean_miscall_rate') is not None:
                _append_metric(
                    threshold, assembly_name, 'whoa_mean_miscall_rate',
                    whoa['mean_miscall_rate']
                )
            for bin_row in whoa.get('depth_bins', []):
                _append_metric(
                    threshold, assembly_name, 'whoa_mean_m',
                    bin_row.get('mean_m', np.nan),
                    mean_dp=bin_row.get('mean_dp', np.nan),
                    lo95=bin_row.get('lo95', np.nan),
                    hi95=bin_row.get('hi95', np.nan),
                    n_bin=bin_row.get('n_bin', np.nan),
                )

    if threshold_sample_rows:
        all_stats.extend(threshold_sample_rows)
    
    # Extract HWE per-sample statistics
    if hwe_file and vcf_file and os.path.exists(hwe_file) and os.path.exists(vcf_file):
        hwe_sample_data = extract_hwe_per_sample(hwe_file, vcf_file, json_file, replicate_suffix)
        # Add threshold to each HWE stat
        for stat in hwe_sample_data:
            stat['threshold'] = threshold
        hwe_stats.extend(hwe_sample_data)

        try:
            position_map = build_vcf_position_map(vcf_file)
            loci = parse_hardy_file(hwe_file, position_map=position_map)
            hwe_summary = summarize_hwe(loci, alpha=0.05)
            valid = max(1, hwe_summary.get('valid_loci', 0))
            _append_metric(
                threshold, assembly_name, 'hwe_prop_rejected',
                hwe_summary.get('sig_deviations', 0) / valid if hwe_summary.get('valid_loci', 0) else np.nan
            )
            _append_metric(
                threshold, assembly_name, 'hwe_prop_deficit',
                hwe_summary.get('sig_deficit', 0) / valid if hwe_summary.get('valid_loci', 0) else np.nan
            )
            _append_metric(
                threshold, assembly_name, 'hwe_prop_excess',
                hwe_summary.get('sig_excess', 0) / valid if hwe_summary.get('valid_loci', 0) else np.nan
            )
        except Exception as e:
            print(f"Warning: global HWE summary failed at threshold {threshold}: {e}", file=sys.stderr)

# Cross-threshold PCA/relatedness stability summaries
def _impute_by_col_mean(mat):
    if mat.size == 0:
        return mat
    col_means = np.nanmean(mat, axis=0)
    col_means = np.where(np.isnan(col_means), 0.0, col_means)
    inds = np.where(np.isnan(mat))
    out = mat.copy()
    out[inds] = np.take(col_means, inds[1])
    return out

for threshold, payload in matrix_by_threshold.items():
    mat = payload['dosage_matrix']
    assembly_name = payload['assembly_name']
    if mat.ndim != 2 or mat.shape[0] < 3 or mat.shape[1] < 2:
        continue
    try:
        x = _impute_by_col_mean(mat)
        x = x - np.mean(x, axis=0, keepdims=True)
        _, svals, _ = np.linalg.svd(x, full_matrices=False)
        var = (svals ** 2)
        var_ratio = var / np.sum(var) if np.sum(var) > 0 else np.zeros_like(var)
        pc1 = float(var_ratio[0]) if var_ratio.size > 0 else np.nan
        pc2 = float(var_ratio[1]) if var_ratio.size > 1 else np.nan
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'pca_pc1_variance_explained',
            'value': pc1
        })
        vcf_stats_list.append({
            'threshold': threshold,
            'assembly_name': assembly_name,
            'metric': 'pca_pc2_variance_explained',
            'value': pc2
        })
    except Exception as e:
        print(f"Warning: PCA summary failed at threshold {threshold}: {e}", file=sys.stderr)

for t1, t2 in combinations(sorted(matrix_by_threshold.keys()), 2):
    d1 = matrix_by_threshold[t1]
    d2 = matrix_by_threshold[t2]
    a1 = d1['assembly_name']
    a2 = d2['assembly_name']
    s1 = d1['sample_names']
    s2 = d2['sample_names']
    common = [s for s in s1 if s in set(s2)]
    if len(common) < 3:
        continue
    idx1 = [s1.index(s) for s in common]
    idx2 = [s2.index(s) for s in common]
    m1 = d1['dosage_matrix'][:, idx1].T  # samples x loci
    m2 = d2['dosage_matrix'][:, idx2].T
    try:
        x1 = _impute_by_col_mean(m1)
        x2 = _impute_by_col_mean(m2)
        # Use sample-wise Euclidean distance matrices as relatedness proxy
        dist1 = np.sqrt(((x1[:, None, :] - x1[None, :, :]) ** 2).mean(axis=2))
        dist2 = np.sqrt(((x2[:, None, :] - x2[None, :, :]) ** 2).mean(axis=2))
        triu = np.triu_indices(len(common), k=1)
        v1 = dist1[triu]
        v2 = dist2[triu]
        if v1.size > 2 and np.std(v1) > 0 and np.std(v2) > 0:
            corr = float(np.corrcoef(v1, v2)[0, 1])
            # assign to higher threshold for plotting trend
            vcf_stats_list.append({
                'threshold': t2,
                'assembly_name': a2,
                'metric': 'relatedness_matrix_correlation_prev',
                'value': corr
            })
    except Exception as e:
        print(f"Warning: relatedness stability failed for {t1} vs {t2}: {e}", file=sys.stderr)

sorted_thresholds = sorted(variant_ids_by_threshold.keys())
for t1, t2 in zip(sorted_thresholds, sorted_thresholds[1:]):
    ids1 = variant_ids_by_threshold.get(t1) or set()
    ids2 = variant_ids_by_threshold.get(t2) or set()
    union = ids1 | ids2
    jaccard = (len(ids1 & ids2) / len(union)) if union else np.nan
    assembly_name = matrix_by_threshold.get(t2, {}).get('assembly_name', f'optimization_{t2}')
    _append_metric(t2, assembly_name, 'variant_jaccard_prev', jaccard)

r80_by_threshold = {}
assembly_by_threshold = {}
for row in vcf_stats_list:
    metric = row.get('metric')
    t = row.get('threshold')
    if t is None:
        continue
    if row.get('assembly_name'):
        assembly_by_threshold.setdefault(t, row['assembly_name'])
    if metric == 'loci_retained_ge80pct':
        r80_by_threshold[t] = row.get('value', np.nan)
sorted_r80 = sorted(r80_by_threshold)
for t1, t2 in zip(sorted_r80, sorted_r80[1:]):
    v1 = r80_by_threshold.get(t1)
    v2 = r80_by_threshold.get(t2)
    if v1 is None or v2 is None or (isinstance(v1, float) and np.isnan(v1)) or (isinstance(v2, float) and np.isnan(v2)):
        continue
    assembly_name = assembly_by_threshold.get(t2, f'optimization_{t2}')
    _append_metric(t2, assembly_name, 'r80_delta', float(v2) - float(v1))

# Save to CSV files
if all_stats:
    df = pd.DataFrame(all_stats)
    df.to_csv(output_csv, index=False)
    print(f"Extracted {len(df)} sample records from {len(thresholds)} thresholds")
    print(f"Saved to {output_csv}")
else:
    # Create empty CSV
    empty_df = pd.DataFrame(columns=[
        'threshold', 'assembly_name', 'sample',
        's3_merged_pairs', 's3_clusters_total', 's3_clusters_hidepth',
        's3_avg_depth_total', 's3_avg_depth_mj', 's3_avg_depth_stat',
        's3_sd_depth_total', 's3_sd_depth_mj', 's3_sd_depth_stat',
        's3_filtered_bad_align',
        's4_hetero_est', 's4_error_est',
        's5_clusters_total', 's5_filtered_by_depth', 's5_filtered_by_maxH',
        's5_filtered_by_maxAlleles', 's5_filtered_by_maxN', 's5_reads_consens',
        's5_nsites', 's5_nhetero', 's5_heterozygosity', 'hwe_deviated_loci',
        'vcf_call_rate', 'vcf_missing_rate', 'vcf_observed_het', 'vcf_expected_het',
        'vcf_fis', 'singleton_proportion'
    ])
    empty_df.to_csv(output_csv, index=False)
    print(f"No data found. Created empty CSV: {output_csv}")

# Save VCF statistics
if vcf_stats_list and output_vcf_csv:
    vcf_df = pd.DataFrame(vcf_stats_list)
    vcf_df.to_csv(output_vcf_csv, index=False)
    print(f"Extracted {len(vcf_df)} VCF locus records from {len(thresholds)} thresholds")
    print(f"VCF stats saved to {output_vcf_csv}")
elif output_vcf_csv:
    # Create empty VCF CSV
    empty_vcf_df = pd.DataFrame(columns=['threshold', 'assembly_name', 'metric', 'value'])
    empty_vcf_df.to_csv(output_vcf_csv, index=False)

# Save HWE statistics
if hwe_stats and output_hwe_csv:
    hwe_df = pd.DataFrame(hwe_stats)
    hwe_df.to_csv(output_hwe_csv, index=False)
    print(f"Extracted {len(hwe_df)} HWE per-sample records from {len(thresholds)} thresholds")
    print(f"HWE stats saved to {output_hwe_csv}")
elif output_hwe_csv:
    # Create empty HWE CSV
    empty_hwe_df = pd.DataFrame(columns=['threshold', 'assembly_name', 'sample', 'metric', 'value'])
    empty_hwe_df.to_csv(output_hwe_csv, index=False)

# Save detailed replicate statistics
if output_replicate_csv:
    if replicate_stats:
        rep_df = pd.DataFrame(replicate_stats)
        rep_df.to_csv(output_replicate_csv, index=False)
        print(f"Detailed replicate stats saved to {output_replicate_csv}")
    else:
        empty_rep_df = pd.DataFrame(columns=[
            'threshold', 'assembly_name', 'sample', 'replicate',
            'compared_loci', 'discordant_loci', 'replicate_discordance',
            'loci_dropout', 'allele_dropout_hom_het', 'allele_dropout_hom_hom',
            'snp_error'
        ])
        empty_rep_df.to_csv(output_replicate_csv, index=False)

# Step 2: Generate plots using plot_from_csv.py
print("Generating plots using plot_from_csv.py...")
plot_script = os.path.join(script_dir, 'plot_from_csv.py')
if os.path.exists(plot_script):
    try:
        result = subprocess.run(
            [
                sys.executable, plot_script,
                "--layout", "diagnostic",
                output_csv, output_pdf,
                "--vcf", output_vcf_csv or "",
                "--hwe", output_hwe_csv or "",
            ],
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        print(f"Plots saved to {output_pdf}")
        if output_supp_pdf:
            supp_cmd = [
                sys.executable, plot_script,
                "--layout", "supplementary",
                output_csv, output_supp_pdf,
                "--vcf", output_vcf_csv or "",
                "--hwe", output_hwe_csv or "",
            ]
            if output_replicate_csv:
                supp_cmd.extend(["--replicate", output_replicate_csv])
            result = subprocess.run(supp_cmd, capture_output=True, text=True, check=True)
            print(result.stdout)
            if result.stderr:
                print(result.stderr, file=sys.stderr)
            print(f"Supplementary figure saved to {output_supp_pdf}")
    except subprocess.CalledProcessError as e:
        print(f"Error generating plots: {e}", file=sys.stderr)
        print(f"STDOUT: {e.stdout}", file=sys.stderr)
        print(f"STDERR: {e.stderr}", file=sys.stderr)
        raise
    except Exception as e:
        print(f"Error calling plot_from_csv.py: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        raise
else:
    print(f"Warning: plot_from_csv.py not found at {plot_script}", file=sys.stderr)
    raise FileNotFoundError(f"Missing plotting script: {plot_script}")

print(f"\nOptimization comparison complete!")
print(f"  Statistics CSV: {output_csv}")
if output_vcf_csv:
    print(f"  VCF Statistics CSV: {output_vcf_csv}")
if output_hwe_csv:
    print(f"  HWE Statistics CSV: {output_hwe_csv}")
if output_replicate_csv:
    print(f"  Replicate Statistics CSV: {output_replicate_csv}")
print(f"  Plots PDF: {output_pdf}")
if output_supp_pdf:
    print(f"  Supplementary figure: {output_supp_pdf}")
