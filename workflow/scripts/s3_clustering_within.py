import os
import sys
import ipyrad as ip

# Redirect stdout to log file
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

data = ip.load_json(snakemake.input.json)

# Set parameters
proj_config = snakemake.params.project_config
ipyrad_params = proj_config["ipyrad"]
threshold = snakemake.params.threshold if snakemake.params.threshold else ipyrad_params["clust_threshold"]
data.params.clust_threshold = threshold
data.params.max_alleles_consens = ipyrad_params.get("max_alleles_consens", 2)

# Set number of cores for parallel processing
data.ipcluster["cores"] = snakemake.threads

# Set the working directories to match Snakefile's expected output structure
# Output dir is like: {project_dir}/{assembly_dir}/{assembly_name}/{assembly_name}_clust_{threshold}/
# Stats file should be written directly to the output directory
output_dir = os.path.dirname(snakemake.output.stats)
data.dirs.clusts = output_dir

# Run step 3 - ipyrad automatically creates stats files in the output_dir
data.run("3", force=True, auto=True)

# Verify the output file was created
expected_output = snakemake.output.stats
if not os.path.exists(expected_output):
    raise FileNotFoundError(f"Expected output file not created: {expected_output}")
