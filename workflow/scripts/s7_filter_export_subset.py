import os
import sys
import ipyrad as ip

# Redirect stdout to log file
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

data = ip.load_json(snakemake.params.json_path)
proj_config = snakemake.params.project_config
ipyrad_params = proj_config["ipyrad"]
subset_name = snakemake.wildcards.subset_name

# Set number of cores for parallel processing
data.ipcluster["cores"] = snakemake.threads

# Identify failed samples (state < 6 have zero loci after previous steps)
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
            "Rerun steps 3-6 for this assembly to refresh {assembly}.json before running subset export."
        )
    raise Exception(f"ERROR: All {initial_count} samples were filtered out (zero loci). Cannot proceed with step 7.")

# Get the sample list for this subset
subsets = proj_config.get("step7_subsets", {})
if subset_name not in subsets:
    raise ValueError(f"Subset '{subset_name}' not found in step7_subsets config")

sample_list = subsets[subset_name]

# Handle both comma-separated and space-separated formats
if isinstance(sample_list, str):
    # Try comma first, then space
    if ',' in sample_list:
        samples = [s.strip() for s in sample_list.split(',') if s.strip()]
    else:
        samples = sample_list.split()
else:
    samples = sample_list

# Filter to only existing samples
existing = [s for s in samples if s in data.samples]
missing = [s for s in samples if s not in data.samples]
eligible = [s for s in existing if s in samples_to_keep]
ineligible = [s for s in existing if s not in samples_to_keep]

if missing:
    print(f"WARNING: {len(missing)} sample(s) from subset '{subset_name}' not found in assembly:")
    for sample in missing[:10]:  # Show first 10 missing samples
        print(f"  - {sample}")
    if len(missing) > 10:
        print(f"  ... and {len(missing) - 10} more")

if not existing:
    raise Exception(f"ERROR: No valid samples found for subset '{subset_name}'. Cannot proceed.")

if ineligible:
    print(
        f"WARNING: {len(ineligible)} sample(s) in subset '{subset_name}' are below step 6 and will be excluded."
    )
    for sample in ineligible[:10]:
        state = data.samples[sample].stats.state
        print(f"  - {sample} (state: {state})")
    if len(ineligible) > 10:
        print(f"  ... and {len(ineligible) - 10} more")

if not eligible:
    raise Exception(
        f"ERROR: Subset '{subset_name}' has no samples at step >= 6. "
        "Cannot run subset step 7 without consensus loci."
    )

print(
    f"Creating subset '{subset_name}' with {len(eligible)} eligible samples "
    f"(out of {len(samples)} specified)"
)

# Branch to create subset
subset_assembly_name = f"{snakemake.wildcards.assembly_name}_{subset_name}"
subset_data = data.branch(subset_assembly_name, subsamples=eligible)
subset_data.params.min_samples_locus = max(1, int(len(eligible) * ipyrad_params.get("min_samples_percent", 50) / 100))
subset_data.ipcluster["cores"] = snakemake.threads

# Set output formats
subset_data.params.output_formats = ipyrad_params.get("output_formats", "p,v")

# Run step 7
subset_data.run("7", force=True, auto=True)

