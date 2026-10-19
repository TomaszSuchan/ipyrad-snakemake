#!/usr/bin/env python3
"""
Extract per-sample statistics from optimization assembly JSON files.
Excludes replicate samples and saves to CSV for plotting.

This module provides functions for extracting statistics and is only used
as an importable module by report_optimization.py.
"""
import sys
import os
import json
import ipyrad as ip
import pandas as pd
import numpy as np
from collections import defaultdict

from qc_common import (
    hdplot_is_paralog,
    parse_allele_depths,
)

def is_replicate(sample_name, replicate_suffix="_repl"):
    """Check if sample is a replicate."""
    return replicate_suffix in sample_name

def extract_hwe_per_sample(hwe_file, vcf_file, json_file, replicate_suffix="_repl", alpha=0.05):
    """Extract per-sample HWE deviation rates.
    
    For each sample, counts how many HWE-deviated loci they are present in.
    Returns a list of dictionaries with threshold, assembly_name, sample, and hwe_deviated_loci.
    """
    hwe_data = []
    
    if not os.path.exists(hwe_file) or not os.path.exists(vcf_file):
        return hwe_data
    
    try:
        # Parse HWE file to get deviated loci (tab-separated format)
        deviated_loci = set()
        with open(hwe_file, 'r') as f:
            header = f.readline().strip()
            # Expected columns: CHR, POS, OBS(HOM1/HET/HOM2), E(HOM1/HET/HOM2), ChiSq_HWE, P_HWE, P_HET_DEFICIT, P_HET_EXCESS
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        chrom = parts[0]
                        pos = parts[1]
                        try:
                            p_hwe = float(parts[5]) if parts[5] != 'nan' else np.nan
                            if not np.isnan(p_hwe) and p_hwe < alpha:
                                deviated_loci.add((chrom, pos))
                        except (ValueError, IndexError):
                            continue
        
        if len(deviated_loci) == 0:
            return hwe_data
        
        # Get threshold and assembly name from JSON
        threshold = None
        assembly_name = 'unknown'
        if json_file and os.path.exists(json_file):
            try:
                data = ip.load_json(json_file)
                if hasattr(data, 'name'):
                    assembly_name = data.name
                    if 'optimization_' in assembly_name:
                        try:
                            threshold = int(assembly_name.split('_')[-1])
                        except (ValueError, IndexError):
                            pass
                if threshold is None and hasattr(data, 'params'):
                    if hasattr(data.params, 'clust_threshold'):
                        threshold = int(data.params.clust_threshold * 100)
            except:
                pass
        
        # Parse VCF to count deviated loci per sample
        import gzip
        opener = gzip.open if vcf_file.endswith('.gz') else open
        
        sample_names = []
        sample_deviated_counts = {}
        
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    if len(fields) > 9:
                        sample_names = fields[9:]
                        # Initialize counters, excluding replicates
                        for sample_name in sample_names:
                            if not is_replicate(sample_name, replicate_suffix):
                                sample_deviated_counts[sample_name] = 0
                    continue
                
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                
                chrom = fields[0]
                pos = fields[1]
                
                # Check if this locus deviates from HWE
                if (chrom, pos) in deviated_loci:
                    # Count samples present at this deviated locus
                    samples = fields[9:]
                    for i, sample_gt in enumerate(samples):
                        if i < len(sample_names):
                            sample_name = sample_names[i]
                            if not is_replicate(sample_name, replicate_suffix):
                                # Check if sample has data (not missing)
                                gt_fields = sample_gt.split(':')
                                genotype = gt_fields[0]
                                is_missing = genotype.startswith('./.') or genotype.startswith('.|.')
                                if not is_missing:
                                    sample_deviated_counts[sample_name] += 1
        
        # Create data rows
        for sample_name, count in sample_deviated_counts.items():
            hwe_data.append({
                'threshold': threshold,
                'assembly_name': assembly_name,
                'sample': sample_name,
                'hwe_deviated_loci': count
            })
    
    except Exception as e:
        print(f"Warning: Error extracting HWE per-sample stats: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    
    return hwe_data

def extract_vcf_stats(vcf_file, json_file, replicate_suffix="_repl"):
    """Extract per-locus statistics from VCF file.
    
    Returns a dict with summary, per-locus, and per-sample metrics:
    - summary counts: fragments, variants, allelic categories
    - per-locus values: NS, DP, MAF
    - per-sample: call_rate, missing_rate, observed_het, expected_het, fis
    - missingness retention at cutoffs 50/70/80/90% (80% is Paris r80)
    - singleton / doubleton proportions
    - genotype matrix (dosage) for downstream PCA/relatedness checks
    - McKinney HDplot H/D and the proportion of biallelic SNPs flagged as paralogs
    """
    stats = {
        'n_fragments': 0,
        'n_variants': 0,
        'n_biallelic': 0,
        'n_triallelic': 0,
        'n_quadriallelic': 0,
        'n_gt4_allelic': 0,
        'ns_values': [],
        'dp_values': [],
        'maf_values': [],
        'singleton_proportion': np.nan,
        'doubleton_proportion': np.nan,
        'missingness_retention': {},
        'per_sample': {},
        'sample_names': [],
        'dosage_matrix': np.empty((0, 0), dtype=float),
        'variant_ids': set(),
        'hdplot_h': [],
        'hdplot_d': [],
        'hdplot_n_tested': 0,
        'hdplot_n_flagged': 0,
        'hdplot_prop_flagged': np.nan,
    }
    
    if not os.path.exists(vcf_file):
        return stats
    
    try:
        import gzip
        
        # Determine if file is gzipped
        opener = gzip.open if vcf_file.endswith('.gz') else open
        
        fragments = set()
        variant_ids = set()
        variants = []
        missingness_cutoffs = [0.5, 0.7, 0.8, 0.9]
        retained_by_cutoff = {c: 0 for c in missingness_cutoffs}
        sample_names = []
        all_sample_names = []
        hdplot_h = []
        hdplot_d = []
        hdplot_n_tested = 0
        hdplot_n_flagged = 0
        dose_rows = []
        per_sample_called = defaultdict(int)
        per_sample_nonmissing_hets = defaultdict(int)
        per_sample_nonmissing_loci = defaultdict(int)
        total_nonmissing_alleles = 0
        singleton_sites = 0
        doubleton_sites = 0
        
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    header_fields = line.strip().split('\t')
                    if len(header_fields) > 9:
                        all_sample_names = header_fields[9:]
                        sample_names = [
                            s for s in all_sample_names
                            if not is_replicate(s, replicate_suffix)
                        ]
                    continue
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                
                chrom = fields[0]  # Column 1: CHROM
                pos = fields[1]    # Column 2: POS
                ref = fields[3]    # Column 4: REF
                alt = fields[4]    # Column 5: ALT
                info = fields[7]   # Column 8: INFO
                
                fragments.add(chrom)
                variant_ids.add(f"{chrom}:{pos}")
                
                # Count alleles: REF + ALT alleles
                # ALT can be "." (missing), single allele, or comma-separated multiple alleles
                if alt == '.':
                    n_alleles = 1  # Only REF
                else:
                    n_alt_alleles = alt.count(',') + 1  # Count commas + 1
                    n_alleles = 1 + n_alt_alleles  # REF + ALT alleles
                
                # Categorize by number of alleles
                if n_alleles == 2:
                    stats['n_biallelic'] += 1
                elif n_alleles == 3:
                    stats['n_triallelic'] += 1
                elif n_alleles == 4:
                    stats['n_quadriallelic'] += 1
                elif n_alleles > 4:
                    stats['n_gt4_allelic'] += 1
                
                # Parse INFO field for NS and DP
                ns = None
                dp = None
                for info_field in info.split(';'):
                    if '=' in info_field:
                        key, value = info_field.split('=', 1)
                        if key == 'NS':
                            try:
                                ns = int(value)
                            except ValueError:
                                pass
                        elif key == 'DP':
                            try:
                                dp = int(value)
                            except ValueError:
                                pass
                
                if ns is not None:
                    stats['ns_values'].append(ns)
                    # Treat NS as number of samples with data at this site.
                    if sample_names:
                        call_rate = ns / max(1, len(sample_names))
                        for cutoff in missingness_cutoffs:
                            if call_rate >= cutoff:
                                retained_by_cutoff[cutoff] += 1
                if dp is not None:
                    stats['dp_values'].append(dp)

                # Parse genotypes for non-replicate samples (keep header order)
                gt_fields = fields[9:] if len(fields) > 9 else []
                format_keys = fields[8].split(':') if len(fields) > 8 else []
                dosage_row = []
                nonmissing_dosages = []
                n_called_locus = 0
                n_het_locus = 0
                het_ref_fractions = []
                is_biallelic_snp = (
                    n_alleles == 2
                    and len(ref) == 1
                    and len(alt) == 1
                    and ref != '.'
                    and alt != '.'
                )
                for sname, sample_gt in zip(all_sample_names, gt_fields):
                    if is_replicate(sname, replicate_suffix):
                        continue
                    gt_main = sample_gt.split(':')[0] if sample_gt else './.'
                    alleles = gt_to_alleles(gt_main)
                    if -1 in alleles:
                        dosage_row.append(np.nan)
                        continue
                    dosage = float(alleles[0] + alleles[1])
                    dosage_row.append(dosage)
                    nonmissing_dosages.append(dosage)
                    per_sample_called[sname] += 1
                    per_sample_nonmissing_loci[sname] += 1
                    n_called_locus += 1
                    if dosage == 1.0:
                        per_sample_nonmissing_hets[sname] += 1
                        n_het_locus += 1
                        depths = parse_allele_depths(sample_gt, format_keys, ref, alt)
                        if depths is not None:
                            ref_c, alt_c = depths
                            tot = ref_c + alt_c
                            if tot > 0:
                                het_ref_fractions.append(ref_c / tot)

                if sample_names:
                    if len(dosage_row) < len(sample_names):
                        dosage_row.extend([np.nan] * (len(sample_names) - len(dosage_row)))
                    elif len(dosage_row) > len(sample_names):
                        dosage_row = dosage_row[:len(sample_names)]
                    dose_rows.append(dosage_row)

                if is_biallelic_snp and n_called_locus >= 5 and n_het_locus >= 2 and het_ref_fractions:
                    h_val = n_het_locus / n_called_locus
                    d_val = (float(np.mean(het_ref_fractions)) - 0.5) * 100.0
                    hdplot_h.append(h_val)
                    hdplot_d.append(d_val)
                    hdplot_n_tested += 1
                    if hdplot_is_paralog(h_val, d_val):
                        hdplot_n_flagged += 1

                # MAF and singleton/doubleton stats from dosage
                if nonmissing_dosages:
                    alt_sum = float(np.sum(nonmissing_dosages))
                    called_samples = len(nonmissing_dosages)
                    alt_count = int(round(alt_sum))
                    total_alleles = 2 * called_samples
                    ref_count = total_alleles - alt_count
                    minor_count = min(alt_count, ref_count)
                    if total_alleles > 0:
                        maf = minor_count / total_alleles
                        stats['maf_values'].append(maf)
                        total_nonmissing_alleles += total_alleles
                        if minor_count == 1:
                            singleton_sites += 1
                        elif minor_count == 2:
                            doubleton_sites += 1
                
                variants.append({
                    'chrom': chrom,
                    'n_alleles': n_alleles,
                    'ns': ns,
                    'dp': dp
                })
        
        stats['n_fragments'] = len(fragments)
        stats['n_variants'] = len(variants)
        stats['variant_ids'] = variant_ids
        stats['missingness_retention'] = retained_by_cutoff
        if stats['n_variants'] > 0:
            stats['singleton_proportion'] = singleton_sites / stats['n_variants']
            stats['doubleton_proportion'] = doubleton_sites / stats['n_variants']
        else:
            stats['singleton_proportion'] = np.nan
            stats['doubleton_proportion'] = np.nan

        # Per-sample metrics
        stats['sample_names'] = sample_names
        for sname in sample_names:
            n_called = per_sample_called.get(sname, 0)
            n_nonmissing = per_sample_nonmissing_loci.get(sname, 0)
            call_rate = n_called / max(1, stats['n_variants'])
            missing_rate = 1.0 - call_rate
            observed_het = (
                per_sample_nonmissing_hets.get(sname, 0) / max(1, n_nonmissing)
                if n_nonmissing > 0 else np.nan
            )
            # Approximate expected heterozygosity from sample mean allele frequency.
            # dosage/2 gives allele frequency at each non-missing locus for diploids.
            if sample_names and dose_rows:
                col_idx = sample_names.index(sname)
                sample_dosages = np.array(
                    [row[col_idx] for row in dose_rows if not np.isnan(row[col_idx])],
                    dtype=float
                )
                if sample_dosages.size > 0:
                    p_hat = np.mean(sample_dosages / 2.0)
                    expected_het = 2 * p_hat * (1 - p_hat)
                else:
                    expected_het = np.nan
            else:
                expected_het = np.nan
            if expected_het is not None and not np.isnan(expected_het) and expected_het > 0:
                fis = 1.0 - (observed_het / expected_het)
            else:
                fis = np.nan
            stats['per_sample'][sname] = {
                'call_rate': call_rate,
                'missing_rate': missing_rate,
                'observed_het': observed_het,
                'expected_het': expected_het,
                'fis': fis,
            }

        if dose_rows and sample_names:
            stats['dosage_matrix'] = np.array(dose_rows, dtype=float)
        else:
            stats['dosage_matrix'] = np.empty((0, 0), dtype=float)

        stats['hdplot_h'] = hdplot_h
        stats['hdplot_d'] = hdplot_d
        stats['hdplot_n_tested'] = hdplot_n_tested
        stats['hdplot_n_flagged'] = hdplot_n_flagged
        stats['hdplot_prop_flagged'] = (
            hdplot_n_flagged / hdplot_n_tested if hdplot_n_tested > 0 else np.nan
        )
        
    except Exception as e:
        print(f"Warning: Error reading VCF {vcf_file}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    
    return stats


def gt_to_alleles(gt_str):
    """Convert GT string to allele indices, handling missing data."""
    if gt_str is None or gt_str in ['./.', '.|.', '.']:
        return [-1, -1]

    if '|' in gt_str:
        alleles = gt_str.split('|')
    elif '/' in gt_str:
        alleles = gt_str.split('/')
    else:
        return [int(gt_str) if gt_str != '.' else -1, -1]

    result = []
    for allele in alleles:
        if allele == '.':
            result.append(-1)
        else:
            result.append(int(allele))

    if len(result) == 1:
        result.append(result[0])
    elif len(result) > 2:
        result = result[:2]

    return result

def extract_sample_stats(json_file, replicate_suffix="_repl"):
    """Extract per-sample statistics from JSON file, excluding replicates."""
    all_data = []
    
    if not os.path.exists(json_file):
        return all_data
    
    try:
        data = ip.load_json(json_file)
        
        # Get threshold from assembly name or params
        threshold = None
        if hasattr(data, 'name'):
            # Extract threshold from name like "optimization_93"
            if 'optimization_' in data.name:
                try:
                    threshold = int(data.name.split('_')[-1])
                except (ValueError, IndexError):
                    pass
        if threshold is None and hasattr(data, 'params'):
            if hasattr(data.params, 'clust_threshold'):
                threshold = int(data.params.clust_threshold * 100)
        
        # Get assembly name
        assembly_name = data.name if hasattr(data, 'name') else 'unknown'
        
        # Extract stats for each sample
        if hasattr(data, 'samples'):
            for sample_name, sample_obj in data.samples.items():
                # Skip replicates
                if is_replicate(sample_name, replicate_suffix):
                    continue
                
                # Get stats_dfs - can be accessed as attribute or dict key
                stats_dfs = None
                if hasattr(sample_obj, 'stats_dfs'):
                    stats_dfs = sample_obj.stats_dfs
                elif isinstance(sample_obj, dict) and 'stats_dfs' in sample_obj:
                    stats_dfs = sample_obj['stats_dfs']
                
                if stats_dfs is None:
                    continue
                
                # Initialize row with metadata
                row = {
                    'threshold': threshold,
                    'assembly_name': assembly_name,
                    'sample': sample_name
                }
                
                # Extract s3 stats (except hidepth_min)
                # stats_dfs can be a dict (from JSON) or an object (from ipyrad)
                if 's3' in stats_dfs:
                    s3 = stats_dfs['s3']
                    # Handle both dict and object access
                    def get_val(obj, key):
                        if isinstance(obj, dict):
                            return obj.get(key)
                        return getattr(obj, key, None)
                    
                    val = get_val(s3, 'merged_pairs')
                    row['s3_merged_pairs'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s3, 'clusters_total')
                    row['s3_clusters_total'] = val if val is not None and pd.notna(val) else None
                    # Skip hidepth_min as requested
                    
                    val = get_val(s3, 'clusters_hidepth')
                    row['s3_clusters_hidepth'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s3, 'avg_depth_total')
                    row['s3_avg_depth_total'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s3, 'avg_depth_mj')
                    row['s3_avg_depth_mj'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s3, 'avg_depth_stat')
                    row['s3_avg_depth_stat'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s3, 'sd_depth_total')
                    row['s3_sd_depth_total'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s3, 'sd_depth_mj')
                    row['s3_sd_depth_mj'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s3, 'sd_depth_stat')
                    row['s3_sd_depth_stat'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s3, 'filtered_bad_align')
                    row['s3_filtered_bad_align'] = val if val is not None and pd.notna(val) else None
                
                # Extract s4 stats
                if 's4' in stats_dfs:
                    s4 = stats_dfs['s4']
                    def get_val(obj, key):
                        if isinstance(obj, dict):
                            return obj.get(key)
                        return getattr(obj, key, None)
                    
                    val = get_val(s4, 'hetero_est')
                    row['s4_hetero_est'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s4, 'error_est')
                    row['s4_error_est'] = val if val is not None and pd.notna(val) else None
                
                # Extract s5 stats
                if 's5' in stats_dfs:
                    s5 = stats_dfs['s5']
                    def get_val(obj, key):
                        if isinstance(obj, dict):
                            return obj.get(key)
                        return getattr(obj, key, None)
                    
                    val = get_val(s5, 'clusters_total')
                    row['s5_clusters_total'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s5, 'filtered_by_depth')
                    row['s5_filtered_by_depth'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s5, 'filtered_by_maxH')
                    row['s5_filtered_by_maxH'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s5, 'filtered_by_maxAlleles')
                    row['s5_filtered_by_maxAlleles'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s5, 'filtered_by_maxN')
                    row['s5_filtered_by_maxN'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s5, 'reads_consens')
                    row['s5_reads_consens'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s5, 'nsites')
                    row['s5_nsites'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s5, 'nhetero')
                    row['s5_nhetero'] = val if val is not None and pd.notna(val) else None
                    
                    val = get_val(s5, 'heterozygosity')
                    row['s5_heterozygosity'] = val if val is not None and pd.notna(val) else None
                
                all_data.append(row)
    
    except Exception as e:
        print(f"Warning: Error reading JSON {json_file}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    
    return all_data

