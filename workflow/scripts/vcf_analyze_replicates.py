#!/usr/bin/env python3
import argparse
import vcfpy
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(
        description='Compare VCF replicates using vcfpy.')
    parser.add_argument(
        '--suffix',
        type=str,
        required=True,
        help='Suffix pattern in replicate names (e.g., "_repl" for samples like sample_repl04a).')
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Input VCF file.')
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output report file.')
    parser.add_argument(
        '--list-samples',
        action='store_true',
        help='List all samples in VCF and exit.')
    return parser.parse_args()

def validate_sample_pairs(samples, suffix):
    """
    Find and validate replicate pairs based on suffix pattern.
    For example, if suffix is "_repl", this finds pairs like:
    - sample_repl04a (replicate) matches with sample (original)

    Returns:
        original_indices: list of indices for original samples
        replicate_indices: list of indices for replicate samples
        replicate_names: list of replicate sample names
    """
    original_indices = []
    replicate_indices = []
    replicate_names = []

    # Find samples with the suffix (these are replicates)
    for idx, sample in enumerate(samples):
        if suffix in sample:
            # Extract the original sample name by removing suffix and library ID
            # e.g., "sample_repl04a" -> "sample"
            base_name = sample.split(suffix)[0]

            # Find matching original sample
            if base_name in samples:
                original_idx = samples.index(base_name)
                original_indices.append(original_idx)
                replicate_indices.append(idx)
                replicate_names.append(sample)
                print(f"Found replicate pair: {base_name} <-> {sample}")

    if not replicate_indices:
        print(f"Warning: No samples found with suffix '{suffix}'")

    return original_indices, replicate_indices, replicate_names

def gt_to_alleles(gt_str):
    """Convert GT string to allele indices, handling missing data."""
    if gt_str is None or gt_str in ['./.', '.|.', '.']:
        return [-1, -1]

    # Handle phased and unphased genotypes
    if '|' in gt_str:
        alleles = gt_str.split('|')
    elif '/' in gt_str:
        alleles = gt_str.split('/')
    else:
        # Single allele case
        return [int(gt_str) if gt_str != '.' else -1, -1]

    # Convert to integers, handling missing alleles
    result = []
    for allele in alleles:
        if allele == '.':
            result.append(-1)
        else:
            result.append(int(allele))

    # Ensure we always return exactly 2 alleles
    if len(result) == 1:
        result.append(result[0])  # Haploid to diploid
    elif len(result) > 2:
        result = result[:2]  # Take first 2 if more than diploid

    return result

def main():
    args = parse_args()
    
    # Read VCF file
    reader = vcfpy.Reader.from_path(args.input)
    samples = reader.header.samples.names
    
    # List samples if requested
    if args.list_samples:
        print(f"All {len(samples)} samples in VCF file:")
        for i, sample in enumerate(samples, 1):
            print(f"  {i:3d}. {sample}")
        reader.close()
        return
    
    # Find and validate replicate pairs
    original_indices, replicate_indices, replicate_names = validate_sample_pairs(samples, args.suffix)
    
    if not replicate_indices:
        print("\nNo valid replicate pairs found. Exiting.")
        reader.close()
        return
    
    # Initialize counters
    num_pairs = len(replicate_indices)
    identical = [0] * num_pairs
    loci_dropped = [0] * num_pairs
    allele_dropped_hom_hom = [0] * num_pairs
    allele_dropped_hom_het = [0] * num_pairs
    snp_error = [0] * num_pairs
    total_loci = 0
    
    # Process each variant
    for record in reader:
        total_loci += 1
        
        # Get genotypes for original and replicate samples
        original_gts = []
        replicate_gts = []
        
        for idx in original_indices:
            call = record.calls[idx]
            gt_str = call.data.GT if hasattr(call.data, 'GT') else None
            original_gts.append(gt_to_alleles(gt_str))
        
        for idx in replicate_indices:
            call = record.calls[idx]
            gt_str = call.data.GT if hasattr(call.data, 'GT') else None
            replicate_gts.append(gt_to_alleles(gt_str))
        
        # Compare each original-replicate pair
        for s in range(min(len(original_gts), len(replicate_gts))):
            l1 = original_gts[s]  # original
            l2 = replicate_gts[s]  # replicate
            
            # Check for missing data first
            l1_missing = (l1[0] == -1 and l1[1] == -1)
            l2_missing = (l2[0] == -1 and l2[1] == -1)
            
            if l1_missing and l2_missing:
                # Both missing - skip (not a dropout)
                continue
            elif l1_missing or l2_missing:
                loci_dropped[s] += 1
            # Check for identical genotypes (order-independent for unphased)
            elif (set(l1) == set(l2) and l1[0] != -1 and l2[0] != -1):
                identical[s] += 1
            # Check for multiallelic variants (allele index > 1 indicates multiallelic)
            elif max(l1 + l2) > 1:
                snp_error[s] += 1
            elif (l1[0] == l1[1]) and (l2[0] == l2[1]) and (l1[0] != l2[0]):
                allele_dropped_hom_hom[s] += 1
            elif ((l1[0] == l1[1]) or (l2[0] == l2[1])) and len(set(l1) | set(l2)) == 2:
                allele_dropped_hom_het[s] += 1
            elif len(set(l1) | set(l2)) > 2:
                snp_error[s] += 1
            else:
                # Fallback: any other mismatch
                snp_error[s] += 1
    
    reader.close()
    
    # Calculate proportions
    loci_dropped_proportion = [round(x / total_loci, 4) for x in loci_dropped]
    allele_dropped_hom_hom_proportion = [round(x / total_loci, 4) for x in allele_dropped_hom_hom]
    allele_dropped_hom_het_proportion = [round(x / total_loci, 4) for x in allele_dropped_hom_het]
    snp_error_proportion = [round(x / total_loci, 4) for x in snp_error]
    
    # Get original sample names for output
    original_names = [samples[idx] for idx in original_indices]

    # Write results to output file
    with open(args.output, 'w') as out:
        out.write("=" * 80 + "\n")
        out.write("VCF Replicate Analysis Report\n")
        out.write("=" * 80 + "\n\n")
        out.write(f"Input VCF: {args.input}\n")
        out.write(f"Suffix pattern: {args.suffix}\n")
        out.write(f"Total loci analyzed: {total_loci}\n")
        out.write(f"Number of replicate pairs found: {len(replicate_names)}\n\n")

        out.write("=" * 80 + "\n")
        out.write("Stats per replicate pair:\n")
        out.write("=" * 80 + "\n\n")
        out.write("Name of the original sample:\t" + "\t".join(original_names) + "\n")
        out.write("Name of the replicated sample:\t" + "\t".join(replicate_names) + "\n")
        out.write("Loci dropout (false hom):\t" + "\t".join(map(str, loci_dropped_proportion)) + "\n")
        out.write("Allele dropout (hom-het):\t" + "\t".join(map(str, allele_dropped_hom_het_proportion)) + "\n")
        out.write("Allele dropout (hom-hom):\t" + "\t".join(map(str, allele_dropped_hom_hom_proportion)) + "\n")
        out.write("SNP error (>2 alleles):\t" + "\t".join(map(str, snp_error_proportion)) + "\n\n")

        if len(replicate_names) > 0:
            out.write("=" * 80 + "\n")
            out.write(f"Mean stats for all {len(replicate_names)} replicates:\n")
            out.write("=" * 80 + "\n\n")
            out.write(f"Loci dropout (false hom):\t{round(sum(loci_dropped_proportion) / len(loci_dropped_proportion), 4)}\n")
            out.write(f"Allele dropout (hom-het):\t{round(sum(allele_dropped_hom_het_proportion) / len(allele_dropped_hom_het_proportion), 4)}\n")
            out.write(f"Allele dropout (hom-hom):\t{round(sum(allele_dropped_hom_hom_proportion) / len(allele_dropped_hom_hom_proportion), 4)}\n")
            out.write(f"SNP error (>2 alleles):\t{round(sum(snp_error_proportion) / len(snp_error_proportion), 4)}\n\n")

            out.write("=" * 80 + "\n")
            out.write("Interpretation:\n")
            out.write("=" * 80 + "\n")
            out.write("- Loci dropout: One sample has data while the other is missing\n")
            out.write("- Allele dropout (het-hom): One sample is heterozygous, the other falsely homozygous\n")
            out.write("- Allele dropout (hom-hom): Both samples are falsely homozygous\n")
            out.write("- SNP error: More than 2 unique alleles across the replicate pair (genotyping error or multiallelic variant)\n")

    print(f"Analysis complete. Report written to: {args.output}")

if __name__ == "__main__":
    main()