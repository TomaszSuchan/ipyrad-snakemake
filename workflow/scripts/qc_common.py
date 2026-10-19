#!/usr/bin/env python3
"""
Common QC analysis functions shared between optimization and final reports.
"""
import os
import re
import gzip
from collections import defaultdict
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
try:
    from sklearn.decomposition import PCA
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    # Fallback: use numpy SVD for PCA

def create_ascii_histogram(values, bins=20, width=60, label=""):
    """Create an ASCII histogram with vertical layout: value then horizontal bar."""
    if not values:
        return f"{label}: No data\n"
    
    values = [v for v in values if v is not None and not np.isnan(v)]
    if not values:
        return f"{label}: No valid data\n"
    
    hist, bin_edges = np.histogram(values, bins=bins)
    max_count = max(hist) if max(hist) > 0 else 1
    
    lines = []
    if label:
        lines.append(f"{label}:")
    
    for i in range(len(hist)):
        bin_start = bin_edges[i]
        bin_end = bin_edges[i + 1]
        count = hist[i]
        bar_length = int((count / max_count) * width)
        bar = '█' * bar_length
        # Format: value range | count | bar
        lines.append(f"  {count:6d} │ {bin_start:8.2f} - {bin_end:8.2f} │{bar}")
    
    lines.append(f"  Min: {min(values):.2f}, Max: {max(values):.2f}, Mean: {np.mean(values):.2f}, Median: {np.median(values):.2f}")
    return '\n'.join(lines) + '\n'

def create_ascii_boxplot(values_list, labels, width=60):
    """Create ASCII boxplots with vertical layout: stats then horizontal bar."""
    if not values_list or not labels:
        return ""
    
    # Filter and validate data
    filtered_data = []
    filtered_labels = []
    for vals, lbl in zip(values_list, labels):
        valid_vals = [v for v in vals if v is not None and not np.isnan(v)]
        if valid_vals:
            filtered_data.append(valid_vals)
            filtered_labels.append(lbl)
    
    if not filtered_data:
        return "No valid data for boxplots\n"
    
    lines = []
    
    for vals, label in zip(filtered_data, filtered_labels):
        q1 = np.percentile(vals, 25)
        median = np.percentile(vals, 50)
        q3 = np.percentile(vals, 75)
        iqr = q3 - q1
        lower_whisker = max(min(vals), q1 - 1.5 * iqr)
        upper_whisker = min(max(vals), q3 + 1.5 * iqr)
        min_val = min(vals)
        max_val = max(vals)
        
        # Scale to width
        if max_val == min_val:
            scale = 1
        else:
            scale = width / (max_val - min_val)
        
        def scale_pos(val):
            return int((val - min_val) * scale)
        
        # Create boxplot line
        line = [' '] * (width + 1)
        
        # Whiskers
        lw_pos = scale_pos(lower_whisker)
        uw_pos = scale_pos(upper_whisker)
        q1_pos = scale_pos(q1)
        med_pos = scale_pos(median)
        q3_pos = scale_pos(q3)
        
        # Draw whiskers
        if lw_pos >= 0 and lw_pos <= width:
            line[lw_pos] = '|'
        if uw_pos >= 0 and uw_pos <= width:
            line[uw_pos] = '|'
        
        # Draw box
        for i in range(q1_pos, q3_pos + 1):
            if 0 <= i <= width:
                line[i] = '─'
        
        # Draw median
        if 0 <= med_pos <= width:
            line[med_pos] = '│'
        
        # Draw box boundaries
        if 0 <= q1_pos <= width:
            line[q1_pos] = '│'
        if 0 <= q3_pos <= width:
            line[q3_pos] = '│'
        
        line_str = ''.join(line)
        # Format: label | stats | bar
        stats_str = f"Min={min_val:.2f} Q1={q1:.2f} Med={median:.2f} Q3={q3:.2f} Max={max_val:.2f}"
        lines.append(f"{label:15s} │ {stats_str} │{line_str}")
    
    lines.append("")
    return '\n'.join(lines)

def count_loci_from_file(loci_file):
    """Count loci from .loci file."""
    if not os.path.exists(loci_file):
        return 0
    try:
        with open(loci_file, 'r') as f:
            count = sum(1 for line in f if line.startswith('//'))
        return count
    except Exception as e:
        print(f"Warning: Error counting loci in {loci_file}: {e}", file=__import__('sys').stderr)
        return 0

def extract_rad_locus(variant_id, pattern=None):
    """
    Extract RAD locus ID from variant ID.
    Supports both patterns: loc(\d+)_pos\d+ and loc(\d+)_
    """
    if pattern is None:
        # Try both patterns
        match = re.match(r'(loc\d+)_pos\d+', variant_id)
        if match:
            return match.group(1)
        match = re.match(r'(loc\d+)_', variant_id)
        if match:
            return match.group(1)
        return None
    else:
        match = re.match(pattern, variant_id)
        if match:
            return match.group(1)
        return None

def identify_replicate_pairs(sample_names, replicate_suffix="_repl"):
    """Identify replicate pairs from sample names."""
    replicate_pairs = {}
    replicate_samples = set()
    
    for sample in sample_names:
        if replicate_suffix in sample:
            base_name = sample.split(replicate_suffix)[0]
            if base_name in sample_names:
                if base_name not in replicate_pairs:
                    replicate_pairs[base_name] = []
                replicate_pairs[base_name].append(sample)
                replicate_samples.add(sample)
                replicate_samples.add(base_name)
    
    return replicate_pairs, replicate_samples

def get_non_replicate_samples(sample_names, replicate_samples):
    """Get list of samples that are not replicates."""
    return [s for s in sample_names if s not in replicate_samples]


def gt_to_alleles(gt_str):
    """Convert GT string to two allele indices, with -1 for missing."""
    if gt_str is None or gt_str in ['./.', '.|.', '.']:
        return [-1, -1]
    if '|' in gt_str:
        alleles = gt_str.split('|')
    elif '/' in gt_str:
        alleles = gt_str.split('/')
    else:
        return [int(gt_str) if gt_str != '.' else -1, -1]
    out = []
    for a in alleles:
        if a == '.':
            out.append(-1)
        else:
            out.append(int(a))
    if len(out) == 1:
        out.append(out[0])
    elif len(out) > 2:
        out = out[:2]
    return out


CATG_INDEX = {"C": 0, "A": 1, "T": 2, "G": 3}
HDPLOT_H_CUTOFF = 0.5
HDPLOT_D_CUTOFF = 5.0  # |mean(ref fraction) - 0.5| in percentage points (SNPfiltR-style)


def parse_allele_depths(sample_gt, format_keys, ref, alt):
    """Return (ref_count, alt_count) from FORMAT AD or ipyrad CATG, else None."""
    if not sample_gt or not format_keys:
        return None
    parts = sample_gt.split(":")
    fmt = {key: parts[i] if i < len(parts) else "" for i, key in enumerate(format_keys)}
    ad = fmt.get("AD", "")
    if ad and ad != ".":
        ads = ad.split(",")
        if len(ads) >= 2:
            try:
                return int(float(ads[0])), int(float(ads[1]))
            except ValueError:
                pass
    catg = fmt.get("CATG", "")
    if catg and catg != "." and len(ref) == 1 and len(alt) == 1:
        try:
            counts = [int(float(x)) for x in catg.split(",")]
        except ValueError:
            return None
        if len(counts) < 4:
            return None
        ri = CATG_INDEX.get(ref[0].upper())
        ai = CATG_INDEX.get(alt[0].upper())
        if ri is None or ai is None:
            return None
        return counts[ri], counts[ai]
    return None


def hdplot_is_paralog(h, d, h_cut=HDPLOT_H_CUTOFF, d_cut=HDPLOT_D_CUTOFF):
    """McKinney et al. 2017-style flag: excess heterozygotes + read-ratio deviation."""
    if h is None or d is None or np.isnan(h) or np.isnan(d):
        return False
    return h > h_cut and abs(d) > d_cut


def calculate_genetic_distance_missingness_from_matrix(dosage_matrix, missing_rates):
    """
    McCartney-Melstad-style correlation: pairwise genetic distance vs pairwise missingness.
    dosage_matrix: variants x samples. missing_rates: length n_samples.
    """
    mat = np.asarray(dosage_matrix, dtype=float)
    miss = np.asarray(missing_rates, dtype=float).reshape(-1)
    if mat.ndim != 2 or mat.shape[1] != miss.size or miss.size < 3:
        return None
    if not np.all(np.isfinite(miss)):
        return None
    x = mat.T  # samples x variants
    missing_per_variant = np.isnan(x).sum(axis=0) / max(1, x.shape[0])
    x = x[:, missing_per_variant < 0.5]
    if x.shape[0] < 3 or x.shape[1] < 50:
        return None
    imputed = x.copy()
    col_means = np.nanmean(imputed, axis=0)
    col_means = np.where(np.isnan(col_means), 0.0, col_means)
    inds = np.where(np.isnan(imputed))
    imputed[inds] = np.take(col_means, inds[1])
    genetic_distances = pdist(imputed, metric="euclidean")
    missingness_diff = pdist(miss.reshape(-1, 1), metric="euclidean")
    if genetic_distances.size < 3 or np.std(genetic_distances) == 0 or np.std(missingness_diff) == 0:
        return None
    corr, _ = pearsonr(genetic_distances, missingness_diff)
    return float(corr)


def compute_replicate_discordance_details(vcf_file, replicate_suffix="_repl"):
    """
    Compute detailed replicate concordance/discordance from VCF.
    Returns:
      {
        "pair_rows": [per-pair dicts],
        "summary": dict
      }
    """
    result = {
        "pair_rows": [],
        "summary": {
            "n_pairs": 0,
            "mean_discordance": np.nan,
            "sd_discordance": np.nan,
            "mean_loci_dropout": np.nan,
            "mean_allele_dropout_hom_het": np.nan,
            "mean_allele_dropout_hom_hom": np.nan,
            "mean_snp_error": np.nan,
            "total_loci": 0,
        }
    }
    if not os.path.exists(vcf_file):
        return result

    opener = gzip.open if vcf_file.endswith('.gz') else open
    try:
        with opener(vcf_file, 'rt') as f:
            sample_names = []
            for line in f:
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    sample_names = fields[9:] if len(fields) > 9 else []
                    break

        if not sample_names:
            return result

        pairs = []
        for idx, sample in enumerate(sample_names):
            if replicate_suffix in sample:
                base = sample.split(replicate_suffix)[0]
                if base in sample_names:
                    pairs.append((base, sample, sample_names.index(base), idx))

        if not pairs:
            return result

        counters = []
        for base, repl, base_idx, repl_idx in pairs:
            counters.append({
                "base": base,
                "replicate": repl,
                "base_idx": base_idx,
                "repl_idx": repl_idx,
                "compared_loci": 0,
                "discordant_loci": 0,
                "loci_dropout": 0,
                "allele_dropout_hom_het": 0,
                "allele_dropout_hom_hom": 0,
                "snp_error": 0,
            })

        total_loci = 0
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                total_loci += 1
                calls = fields[9:]
                for ctr in counters:
                    if ctr["base_idx"] >= len(calls) or ctr["repl_idx"] >= len(calls):
                        continue
                    gt1 = calls[ctr["base_idx"]].split(':')[0]
                    gt2 = calls[ctr["repl_idx"]].split(':')[0]
                    a1 = gt_to_alleles(gt1)
                    a2 = gt_to_alleles(gt2)

                    m1 = (a1[0] == -1 and a1[1] == -1)
                    m2 = (a2[0] == -1 and a2[1] == -1)
                    if m1 and m2:
                        continue
                    if m1 or m2:
                        ctr["loci_dropout"] += 1
                        continue

                    ctr["compared_loci"] += 1
                    if sorted(a1) != sorted(a2):
                        ctr["discordant_loci"] += 1

                    if max(a1 + a2) > 1:
                        ctr["snp_error"] += 1
                    elif (a1[0] == a1[1]) and (a2[0] == a2[1]) and (a1[0] != a2[0]):
                        ctr["allele_dropout_hom_hom"] += 1
                    elif ((a1[0] == a1[1]) or (a2[0] == a2[1])) and len(set(a1) | set(a2)) == 2:
                        ctr["allele_dropout_hom_het"] += 1
                    elif len(set(a1) | set(a2)) > 2:
                        ctr["snp_error"] += 1

        pair_rows = []
        discordances = []
        dropouts = []
        ahh = []
        ahh2 = []
        snperr = []
        for ctr in counters:
            if ctr["compared_loci"] > 0:
                disc = ctr["discordant_loci"] / ctr["compared_loci"]
            else:
                disc = np.nan
            ld = ctr["loci_dropout"] / max(1, total_loci)
            ad_hh = ctr["allele_dropout_hom_hom"] / max(1, total_loci)
            ad_hh2 = ctr["allele_dropout_hom_het"] / max(1, total_loci)
            se = ctr["snp_error"] / max(1, total_loci)
            pair_rows.append({
                "sample": ctr["base"],
                "replicate": ctr["replicate"],
                "compared_loci": ctr["compared_loci"],
                "discordant_loci": ctr["discordant_loci"],
                "replicate_discordance": disc,
                "loci_dropout": ld,
                "allele_dropout_hom_het": ad_hh2,
                "allele_dropout_hom_hom": ad_hh,
                "snp_error": se,
            })
            if not np.isnan(disc):
                discordances.append(disc)
            dropouts.append(ld)
            ahh.append(ad_hh2)
            ahh2.append(ad_hh)
            snperr.append(se)

        result["pair_rows"] = pair_rows
        result["summary"] = {
            "n_pairs": len(pair_rows),
            "mean_discordance": float(np.mean(discordances)) if discordances else np.nan,
            "sd_discordance": float(np.std(discordances)) if discordances else np.nan,
            "mean_loci_dropout": float(np.mean(dropouts)) if dropouts else np.nan,
            "mean_allele_dropout_hom_het": float(np.mean(ahh)) if ahh else np.nan,
            "mean_allele_dropout_hom_hom": float(np.mean(ahh2)) if ahh2 else np.nan,
            "mean_snp_error": float(np.mean(snperr)) if snperr else np.nan,
            "total_loci": total_loci,
        }
        return result
    except Exception as e:
        print(f"Warning: Error computing replicate discordance for {vcf_file}: {e}", file=__import__('sys').stderr)
        return result

def analyze_vcf(vcf_file, id_pattern=r'loc(\d+)_', replicate_suffix="_repl"):
    """
    VCF analysis for QC metrics.
    Returns statistics for loci, variants, missing data, and heterozygosity.
    """
    stats = {
        'total_variants': 0,
        'variants_per_locus': defaultdict(list),
        'missing_data_per_variant': [],
        'missing_data_per_variant_no_repl': [],  # Without replicates
        'missing_data_per_sample': defaultdict(int),
        'total_calls_per_sample': defaultdict(int),
        'missing_data_per_locus': defaultdict(int),
        'total_calls_per_locus': defaultdict(int),
        'sample_names': [],
        'sample_het_sites': defaultdict(int),
        'sample_hom_sites': defaultdict(int),
        'replicate_pairs': {},
        'replicate_samples': set(),
        'non_replicate_samples': [],
    }
    
    if not os.path.exists(vcf_file):
        return stats
    
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    try:
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    if len(fields) > 9:
                        stats['sample_names'] = fields[9:]
                        # Identify replicates
                        replicate_pairs, replicate_samples = identify_replicate_pairs(
                            stats['sample_names'], replicate_suffix
                        )
                        stats['replicate_pairs'] = replicate_pairs
                        stats['replicate_samples'] = replicate_samples
                        stats['non_replicate_samples'] = get_non_replicate_samples(
                            stats['sample_names'], replicate_samples
                        )
                    continue
                
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                
                stats['total_variants'] += 1
                variant_id = fields[2]
                samples = fields[9:]
                
                # Extract locus ID
                locus_id = extract_rad_locus(variant_id, id_pattern)
                if locus_id:
                    stats['variants_per_locus'][locus_id].append(variant_id)
                    stats['total_calls_per_locus'][locus_id] += len(samples)
                
                # Per-variant missing data
                missing_count = 0
                missing_count_no_repl = 0
                non_repl_total = 0
                
                for i, sample_gt in enumerate(samples):
                    sample_name = stats['sample_names'][i] if i < len(stats['sample_names']) else f"Sample_{i}"
                    stats['total_calls_per_sample'][sample_name] += 1
                    
                    gt_fields = sample_gt.split(':')
                    genotype = gt_fields[0]
                    is_missing = genotype.startswith('./.') or genotype.startswith('.|.')
                    
                    if is_missing:
                        stats['missing_data_per_sample'][sample_name] += 1
                        missing_count += 1
                        if locus_id:
                            stats['missing_data_per_locus'][locus_id] += 1
                    else:
                        # Track heterozygosity
                        alleles = genotype.replace('|', '/').split('/')
                        if len(alleles) == 2:
                            if alleles[0] != alleles[1]:
                                stats['sample_het_sites'][sample_name] += 1
                            else:
                                stats['sample_hom_sites'][sample_name] += 1
                    
                    # Count missing for non-replicate samples only
                    if sample_name in stats['non_replicate_samples']:
                        non_repl_total += 1
                        if is_missing:
                            missing_count_no_repl += 1
                
                # Missing data per variant (all samples)
                if len(samples) > 0:
                    missing_prop = missing_count / len(samples)
                    stats['missing_data_per_variant'].append(missing_prop)
                
                # Missing data per variant (excluding replicates)
                if non_repl_total > 0:
                    missing_prop_no_repl = missing_count_no_repl / non_repl_total
                    stats['missing_data_per_variant_no_repl'].append(missing_prop_no_repl)
    
    except Exception as e:
        print(f"Warning: Error parsing VCF {vcf_file}: {e}", file=__import__('sys').stderr)
    
    return stats

def parse_hardy_summary(hardy_file):
    """Parse Hardy-Weinberg summary from summarize_hardy output."""
    if not os.path.exists(hardy_file):
        return None
    
    summary = {
        'total_loci': 0,
        'hwe_deviations': 0,
        'p_value_distribution': []
    }
    
    try:
        with open(hardy_file, 'r') as f:
            for line in f:
                if 'Total loci tested' in line:
                    match = re.search(r'(\d+)', line)
                    if match:
                        summary['total_loci'] = int(match.group(1))
                elif 'Loci with significant HWE deviation' in line:
                    match = re.search(r'(\d+)', line)
                    if match:
                        summary['hwe_deviations'] = int(match.group(1))
                elif 'Mean p-value' in line:
                    match = re.search(r'([\d.]+)', line)
                    if match:
                        summary['mean_p_value'] = float(match.group(1))
    except Exception as e:
        print(f"Warning: Error parsing Hardy file {hardy_file}: {e}", file=__import__('sys').stderr)
    
    return summary if summary['total_loci'] > 0 else None

def parse_whoa_summary(whoa_file):
    """Extract key summary from whoa report."""
    if not os.path.exists(whoa_file):
        return None

    summary = {
        'mean_miscall_rate': None,
        'depth_bins': []
    }

    try:
        with open(whoa_file, 'r') as f:
            lines = f.readlines()

        in_tsv = False
        for line in lines:
            stripped = line.strip()
            if stripped == "WHOA_POSTERIORS_TSV":
                in_tsv = True
                continue
            if in_tsv:
                if not stripped or stripped.startswith("=") or stripped.startswith("-"):
                    if summary['depth_bins']:
                        break
                    continue
                if stripped.startswith("mean_m"):
                    continue
                parts = stripped.split()
                if len(parts) >= 4:
                    try:
                        mean_m = float(parts[0])
                        lo95 = float(parts[1])
                        hi95 = float(parts[2])
                        n_bin = int(float(parts[3]))
                        mean_dp = float(parts[4]) if len(parts) > 4 else None
                        summary['depth_bins'].append({
                            'mean_m': mean_m,
                            'lo95': lo95,
                            'hi95': hi95,
                            'n_bin': n_bin,
                            'mean_dp': mean_dp
                        })
                    except (ValueError, IndexError):
                        if summary['depth_bins']:
                            break

        if not summary['depth_bins']:
            in_table = False
            for line in lines:
                if 'mean_m' in line and 'lo95' in line:
                    in_table = True
                    continue
                if in_table and line.strip() and not line.startswith('-'):
                    parts = line.strip().split()
                    if parts and parts[0].rstrip('.').isdigit():
                        parts = parts[1:]
                    if len(parts) >= 4:
                        try:
                            mean_m = float(parts[0])
                            lo95 = float(parts[1])
                            hi95 = float(parts[2])
                            n_bin = int(float(parts[3]))
                            mean_dp = float(parts[4]) if len(parts) > 4 else None
                            summary['depth_bins'].append({
                                'mean_m': mean_m,
                                'lo95': lo95,
                                'hi95': hi95,
                                'n_bin': n_bin,
                                'mean_dp': mean_dp
                            })
                        except (ValueError, IndexError):
                            pass

        if summary['depth_bins']:
            summary['mean_miscall_rate'] = float(np.mean(
                [b['mean_m'] for b in summary['depth_bins'] if b.get('mean_m') is not None]
            ))
    except Exception as e:
        print(f"Warning: Error parsing whoa file {whoa_file}: {e}", file=__import__('sys').stderr)

    return summary if summary['depth_bins'] else None

def parse_replicate_summary(replicate_file):
    """Extract summary from replicate analysis."""
    if not os.path.exists(replicate_file):
        return None
    
    summary = {
        'n_pairs': 0,
        'mean_loci_dropout': 0.0,
        'mean_allele_dropout': 0.0,
        'mean_snp_error': 0.0
    }
    
    try:
        with open(replicate_file, 'r') as f:
            for line in f:
                if 'Number of replicate pairs found' in line:
                    match = re.search(r'(\d+)', line)
                    if match:
                        summary['n_pairs'] = int(match.group(1))
                elif 'Mean stats for all' in line:
                    continue
                elif 'Loci dropout (false hom):' in line:
                    match = re.search(r'([\d.]+)', line)
                    if match:
                        summary['mean_loci_dropout'] = float(match.group(1))
                elif 'Allele dropout (hom-het):' in line:
                    match = re.search(r'([\d.]+)', line)
                    if match:
                        summary['mean_allele_dropout'] = float(match.group(1))
                elif 'SNP error (>2 alleles):' in line:
                    match = re.search(r'([\d.]+)', line)
                    if match:
                        summary['mean_snp_error'] = float(match.group(1))
    except Exception as e:
        print(f"Warning: Error parsing replicate file {replicate_file}: {e}", file=__import__('sys').stderr)
    
    return summary if summary['n_pairs'] > 0 else None

def calculate_pca_variance(vcf_file, n_components=8, exclude_replicates=True, replicate_suffix="_repl", return_details=False):
    """
    Calculate PCA from VCF and return cumulative variance explained by first n_components axes.
    Returns cumulative variance explained (0-1) or None if calculation fails.
    If return_details is True, return a dict with scores and sample names.
    """
    if not os.path.exists(vcf_file):
        return None
    
    try:
        import vcfpy
        
        # Read VCF and build genotype matrix
        reader = vcfpy.Reader.from_path(vcf_file)
        sample_names = reader.header.samples.names
        
        # Filter samples if excluding replicates
        if exclude_replicates:
            replicate_samples = set()
            for name in sample_names:
                if replicate_suffix in name:
                    replicate_samples.add(name)
            sample_indices = [i for i, name in enumerate(sample_names) if name not in replicate_samples]
            sample_names = [sample_names[i] for i in sample_indices]
        else:
            sample_indices = list(range(len(sample_names)))
        
        if len(sample_names) < 2:
            return None
        
        # Build genotype matrix: rows = variants, cols = samples
        genotype_matrix = []
        variant_count = 0
        
        for record in reader:
            variant_count += 1
            if variant_count > 10000:  # Limit to first 10k variants for performance
                break
            
            row = []
            for i in sample_indices:
                call = record.calls[i]
                if call.is_variant:
                    # Encode genotype: 0=hom ref, 1=het, 2=hom alt
                    gt = call.gt_type
                    if gt is None:
                        row.append(np.nan)
                    else:
                        row.append(float(gt))
                else:
                    row.append(np.nan)
            genotype_matrix.append(row)
        
        reader.close()
        
        if len(genotype_matrix) < n_components:
            return None
        
        # Convert to numpy array
        X = np.array(genotype_matrix).T  # Transpose: rows = samples, cols = variants
        
        # Remove variants with >50% missing data
        missing_per_variant = np.isnan(X).sum(axis=0) / X.shape[0]
        keep_variants = missing_per_variant < 0.5
        X = X[:, keep_variants]
        
        # Remove samples with >50% missing data
        missing_per_sample = np.isnan(X).sum(axis=1) / X.shape[1]
        keep_samples = missing_per_sample < 0.5
        kept_names = [name for name, keep in zip(sample_names, keep_samples) if keep]
        X = X[keep_samples, :]
        
        if X.shape[0] < 2 or X.shape[1] < n_components:
            return None
        
        # Impute missing values with mean
        X_imputed = X.copy()
        for i in range(X.shape[0]):
            row = X[i, :]
            mask = ~np.isnan(row)
            if mask.sum() > 0:
                X_imputed[i, ~mask] = np.nanmean(row[mask])
        
        # Standardize
        X_std = (X_imputed - np.nanmean(X_imputed, axis=0)) / (np.nanstd(X_imputed, axis=0) + 1e-10)
        
        n_comp = min(n_components, X_std.shape[0] - 1, X_std.shape[1])
        # Perform PCA
        if HAS_SKLEARN:
            pca = PCA(n_components=n_comp)
            scores = pca.fit_transform(X_std)
            explained_var = pca.explained_variance_ratio_
        else:
            # Manual PCA using SVD
            U, s, Vt = np.linalg.svd(X_std, full_matrices=False)
            explained_var = (s ** 2) / np.sum(s ** 2)
            explained_var = explained_var[:n_comp]
            scores = U[:, :n_comp] * s[:n_comp]
        
        # Cumulative variance explained
        cum_var = np.sum(explained_var[:n_components])
        if return_details:
            return {
                'cumulative_variance': float(cum_var),
                'explained_variance_ratio': np.asarray(explained_var),
                'sample_names': kept_names,
                'scores': np.asarray(scores),
            }
        return float(cum_var)
        
    except Exception as e:
        import sys
        print(f"Warning: Error calculating PCA for {vcf_file}: {e}", file=sys.stderr)
        return None

def calculate_genetic_distance_missingness_correlation(vcf_file, exclude_replicates=True, replicate_suffix="_repl"):
    """
    Calculate correlation between pairwise genetic distance and missingness.
    Returns correlation coefficient (float) or None if calculation fails.
    """
    if not os.path.exists(vcf_file):
        return None
    
    try:
        import vcfpy
        
        # Read VCF and build genotype matrix
        reader = vcfpy.Reader.from_path(vcf_file)
        sample_names = reader.header.samples.names
        
        # Filter samples if excluding replicates
        if exclude_replicates:
            replicate_samples = set()
            for name in sample_names:
                if replicate_suffix in name:
                    replicate_samples.add(name)
            sample_indices = [i for i, name in enumerate(sample_names) if name not in replicate_samples]
            sample_names = [sample_names[i] for i in sample_indices]
        else:
            sample_indices = list(range(len(sample_names)))
        
        if len(sample_names) < 3:  # Need at least 3 samples for meaningful correlation
            return None
        
        # Build genotype matrix and missingness per sample
        genotype_matrix = []
        missing_per_sample = np.zeros(len(sample_indices))
        total_calls_per_sample = np.zeros(len(sample_indices))
        variant_count = 0
        
        for record in reader:
            variant_count += 1
            if variant_count > 10000:  # Limit to first 10k variants for performance
                break
            
            row = []
            for idx, i in enumerate(sample_indices):
                call = record.calls[i]
                total_calls_per_sample[idx] += 1
                if call.is_variant:
                    gt = call.gt_type
                    if gt is None:
                        row.append(np.nan)
                        missing_per_sample[idx] += 1
                    else:
                        row.append(float(gt))
                else:
                    row.append(0.0)  # Reference allele
            genotype_matrix.append(row)
        
        reader.close()
        
        if len(genotype_matrix) < 100:  # Need sufficient variants
            return None
        
        # Calculate missingness proportion per sample
        missingness = missing_per_sample / total_calls_per_sample
        
        # Convert to numpy array
        X = np.array(genotype_matrix).T  # Transpose: rows = samples, cols = variants
        
        # Remove variants with >50% missing data
        missing_per_variant = np.isnan(X).sum(axis=0) / X.shape[0]
        keep_variants = missing_per_variant < 0.5
        X = X[:, keep_variants]
        
        if X.shape[0] < 3 or X.shape[1] < 100:
            return None
        
        # Impute missing values with mean
        X_imputed = X.copy()
        for i in range(X.shape[0]):
            row = X[i, :]
            mask = ~np.isnan(row)
            if mask.sum() > 0:
                X_imputed[i, ~mask] = np.nanmean(row[mask])
        
        # Calculate pairwise genetic distances (Euclidean distance)
        genetic_distances = pdist(X_imputed, metric='euclidean')
        genetic_dist_matrix = squareform(genetic_distances)
        
        # Calculate pairwise missingness differences
        missingness_diff = pdist(missingness.reshape(-1, 1), metric='euclidean')
        
        # Calculate correlation
        if len(genetic_distances) > 0 and len(missingness_diff) > 0:
            corr, p_value = pearsonr(genetic_distances, missingness_diff)
            return float(corr)
        else:
            return None
        
    except Exception as e:
        import sys
        print(f"Warning: Error calculating genetic distance-missingness correlation for {vcf_file}: {e}", file=sys.stderr)
        return None

