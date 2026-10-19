import os
import sys
import ipyrad as ip

# Redirect stdout to log file
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

data = ip.load_json(snakemake.params.json_path)
# Use the assembly name from wildcards instead of deriving from output path
# to avoid double _outfiles suffix
data.name = snakemake.wildcards.assembly_name
proj_config = snakemake.params.project_config
ipyrad_params = proj_config["ipyrad"]

# Set number of cores for parallel processing
data.ipcluster["cores"] = snakemake.threads

# Filter failed samples (those with state < 6 have zero loci after previous steps)
initial_count = len(data.samples)
samples_to_keep = [s for s in data.samples if data.samples[s].stats.state >= 6]
failed_samples = [s for s in data.samples if data.samples[s].stats.state < 6]

# Report filtered samples
if failed_samples:
    print(f"WARNING: {len(failed_samples)} sample(s) filtered out (zero loci after step 6):")
    for sample in failed_samples:
        state = data.samples[sample].stats.state
        print(f"  - {sample} (state: {state})")

# Handle case when all samples are filtered out
if len(samples_to_keep) == 0:
    max_state = max((data.samples[s].stats.state for s in data.samples), default=float("nan"))
    if max_state <= 2:
        raise Exception(
            f"ERROR: All {initial_count} samples are at state <= 2 in {snakemake.params.json_path}. "
            "This usually means the merged assembly JSON was recreated and now reflects pre-assembly state. "
            "Rerun steps 3-6 for this assembly before running step 7."
        )
    raise Exception(f"ERROR: All {initial_count} samples were filtered out (zero loci). Cannot proceed with step 7.")

# Branch to keep only samples with loci (if any were filtered)
if len(samples_to_keep) < initial_count:
    # Branch with a temporary name, then rename to keep the original assembly name
    # Note: branch() requires newname as first positional argument, not a keyword
    temp_name = f"{data.name}_filtered_temp"
    data = data.branch(temp_name, subsamples=samples_to_keep)
    data.name = snakemake.wildcards.assembly_name  # Restore original name

# Set parameters
n_samples = len(data.samples)
min_samples_locus = max(1, int(n_samples * ipyrad_params.get("min_samples_percent", 50) / 100))
data.params.min_samples_locus = min_samples_locus

# For optimization assemblies, only output VCF (minimal)
# For full assemblies, use all configured formats
is_optimization = "optimization" in snakemake.output.loci
if is_optimization:
    data.params.output_formats = "v"
else:
    data.params.output_formats = ipyrad_params.get("output_formats", "p,v")

# Run step 7 - ipyrad automatically creates stats files
data.run("7", force=True, auto=True)

# Note: Subsets are now handled by separate s7_filter_export_subset rules
# This allows parallel execution of subset creation
