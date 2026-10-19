import ipyrad as ip
import random

opt_config = snakemake.params.project_config["optimization"]
random.seed(opt_config.get("random_seed", 42))

# Load samples and read counts
all_samples = {}
for json_file in snakemake.input.jsons:
    data = ip.load_json(json_file)
    for sample_name in data.samples.keys():
        # Try to get reads_passed_filter (step 2), fall back to reads_raw (step 1)
        try:
            if 'reads_passed_filter' in data.stats.columns:
                read_count = int(data.stats.loc[sample_name, 'reads_passed_filter'])
            else:
                read_count = int(data.stats.loc[sample_name, 'reads_raw'])
        except Exception as e:
            read_count = 0
        all_samples[sample_name] = read_count

# Select samples
explicit = opt_config.get("samples", None)
if explicit:
    selected = [s for s in explicit if s in all_samples]
else:
    # Filter by replicate suffix and bottom percent
    replicate_suffix = opt_config.get("replicate_suffix", "_repl")
    use_replicates = opt_config.get("use_replicates", True)
    remove_bottom_percent = opt_config.get("remove_bottom_percent", 20)
    n_samples = opt_config["n_samples"]

    import re

    all_sample_names = list(all_samples.keys())

    # Step 1: Find all replicates using the configured suffix
    suffix_re = re.escape(replicate_suffix)
    all_replicates = [s for s in all_sample_names if re.search(suffix_re, s, re.IGNORECASE)]

    # Step 2: Find the original samples for each replicate
    replicate_originals = []
    for repl in all_replicates:
        base_name = re.split(suffix_re, repl, maxsplit=1, flags=re.IGNORECASE)[0]
        for sample in all_sample_names:
            if sample.startswith(base_name) and replicate_suffix.lower() not in sample.lower():
                if sample not in replicate_originals:
                    replicate_originals.append(sample)

    # Step 3: Always include replicates and their originals
    always_include = list(set(all_replicates + replicate_originals))

    # Step 4: Get remaining samples (exclude replicates and their originals)
    remaining_samples = [s for s in all_sample_names if s not in always_include]

    # Step 5: Apply bottom percent filtering to remaining samples
    sorted_remaining = sorted(remaining_samples, key=lambda s: all_samples[s])
    cutoff_idx = int(len(sorted_remaining) * remove_bottom_percent / 100)
    excluded_bottom = sorted_remaining[:cutoff_idx]
    candidates_for_selection = sorted_remaining[cutoff_idx:]

    # Step 6: Randomly select n_samples from the filtered remaining pool
    n_to_select = min(n_samples, len(candidates_for_selection))
    additional_selected = random.sample(candidates_for_selection, n_to_select)

    # Step 7: Combine always_include + additional selections
    selected = always_include + additional_selected

# Write outputs
with open(snakemake.output.sample_list, 'w') as f:
    for sample in sorted(selected):
        f.write(f"{sample}\n")

# Calculate statistics for the report
read_counts = list(all_samples.values())
selected_counts = [all_samples[s] for s in selected]

with open(snakemake.output.report, 'w') as f:
    f.write(f"Optimization Sample Selection\n")
    f.write(f"{'='*70}\n\n")

    # Overall statistics
    f.write(f"Total samples available: {len(all_samples)}\n")
    f.write(f"Selected samples: {len(selected)}\n\n")

    # Selection criteria
    f.write(f"Selection Criteria:\n")
    f.write(f"{'-'*70}\n")
    if explicit:
        f.write(f"  Method: Explicit sample list\n")
    else:
        f.write(f"  Method: Replicates + random selection\n")
        f.write(f"  Target additional samples: {opt_config['n_samples']}\n")
        f.write(f"  Remove bottom percent (from non-replicate pool): {remove_bottom_percent}%\n")
        f.write(f"  Random seed: {opt_config.get('random_seed', 42)}\n\n")
        f.write(f"Sample Selection Strategy:\n")
        f.write(f"  1. Always include all replicates: {len(all_replicates)}\n")
        f.write(f"  2. Always include replicate originals: {len(replicate_originals)}\n")
        f.write(f"  3. From remaining {len(remaining_samples)} samples:\n")
        f.write(f"     - Excluded bottom {remove_bottom_percent}%: {len(excluded_bottom)} samples\n")
        f.write(f"     - Selection pool: {len(candidates_for_selection)} samples\n")
        f.write(f"     - Randomly selected: {len(additional_selected)} samples\n")
        f.write(f"  Total selected: {len(selected)} samples\n")
    f.write("\n")

    # Read count statistics
    f.write(f"Read Count Statistics:\n")
    f.write(f"{'-'*70}\n")
    f.write(f"  All samples - Min: {min(read_counts)} | Max: {max(read_counts)} | Mean: {sum(read_counts)//len(read_counts)}\n")
    f.write(f"  Selected    - Min: {min(selected_counts)} | Max: {max(selected_counts)} | Mean: {sum(selected_counts)//len(selected_counts)}\n")
    f.write("\n")

    # Selected samples table
    f.write(f"Selected Samples:\n")
    f.write(f"{'-'*70}\n")
    f.write(f"  {'Sample Name':<40} {'Read Count':>15}\n")
    f.write(f"  {'-'*40} {'-'*15}\n")
    for sample in sorted(selected):
        f.write(f"  {sample:<40} {all_samples[sample]:>15}\n")

    # Show breakdown of selected samples
    if not explicit:
        # Show replicate samples
        if all_replicates:
            f.write(f"\n")
            f.write(f"Replicate Samples (always included):\n")
            f.write(f"{'-'*70}\n")
            f.write(f"  {'Sample Name':<40} {'Read Count':>15}\n")
            f.write(f"  {'-'*40} {'-'*15}\n")
            for sample in sorted(all_replicates):
                f.write(f"  {sample:<40} {all_samples[sample]:>15}\n")

        # Show replicate originals
        if replicate_originals:
            f.write(f"\n")
            f.write(f"Replicate Original Samples (always included):\n")
            f.write(f"{'-'*70}\n")
            f.write(f"  {'Sample Name':<40} {'Read Count':>15}\n")
            f.write(f"  {'-'*40} {'-'*15}\n")
            for sample in sorted(replicate_originals):
                f.write(f"  {sample:<40} {all_samples[sample]:>15}\n")

        # Show additional randomly selected samples
        if additional_selected:
            f.write(f"\n")
            f.write(f"Additional Randomly Selected Samples:\n")
            f.write(f"{'-'*70}\n")
            f.write(f"  {'Sample Name':<40} {'Read Count':>15}\n")
            f.write(f"  {'-'*40} {'-'*15}\n")
            for sample in sorted(additional_selected):
                f.write(f"  {sample:<40} {all_samples[sample]:>15}\n")

        # Show excluded bottom percent samples
        if excluded_bottom:
            f.write(f"\n")
            f.write(f"Excluded Bottom {remove_bottom_percent}% Samples (from non-replicate pool):\n")
            f.write(f"{'-'*70}\n")
            f.write(f"  {'Sample Name':<40} {'Read Count':>15}\n")
            f.write(f"  {'-'*40} {'-'*15}\n")
            for sample in sorted(excluded_bottom, key=lambda s: all_samples[s]):
                f.write(f"  {sample:<40} {all_samples[sample]:>15}\n")
