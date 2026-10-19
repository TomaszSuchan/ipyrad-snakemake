import ipyrad as ip
import pandas as pd

# Collect stats from all libraries
all_stats = []

for json_file in snakemake.input.jsons:
    data = ip.load_json(json_file)
    if hasattr(data, 'stats') and data.stats is not None and len(data.stats) > 0:
        df = data.stats.copy()
        # Add assembly name and raw fastq path columns
        df.insert(0, 'assembly_name', data.params.assembly_name)
        df.insert(1, 'raw_fastq_path', data.params.raw_fastq_path)
        all_stats.append(df)

# Combine all stats (keep all occurrences)
if all_stats:
    combined = pd.concat(all_stats)

    # Export sorted by sample name
    by_name = combined.sort_index()
    by_name.to_csv(snakemake.output.by_name, sep='\t')

    # Export sorted by reads_passed_filter (low to high)
    if 'reads_passed_filter' in combined.columns:
        by_reads = combined.sort_values('reads_passed_filter', ascending=True)
    else:
        by_reads = combined.sort_values('reads_raw', ascending=True)
    by_reads.to_csv(snakemake.output.by_reads, sep='\t')
else:
    # Create empty files if no stats found
    with open(snakemake.output.by_name, 'w') as f:
        f.write("# No samples found\n")
    with open(snakemake.output.by_reads, 'w') as f:
        f.write("# No samples found\n")
