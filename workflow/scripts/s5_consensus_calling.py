import os
import sys
import ipyrad as ip

# Redirect stdout to log file
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

data = ip.load_json(snakemake.params.json_path)

# Set parameters
proj_config = snakemake.params.project_config
ipyrad_params = proj_config["ipyrad"]
data.params.mindepth_statistical = ipyrad_params.get("mindepth_statistical", 6)
data.params.mindepth_majrule = ipyrad_params.get("mindepth_majrule", 6)
data.params.max_Ns_consens = ipyrad_params.get("max_Ns_consens", 0.05)
data.params.max_Hs_consens = ipyrad_params.get("max_Hs_consens", 0.05)

# Set number of cores for parallel processing
data.ipcluster["cores"] = snakemake.threads

# Set working directories to match Snakefile's expected output structure
# Stats file should be written directly to the output directory
output_dir = os.path.dirname(snakemake.output.stats)
data.dirs.consens = output_dir

# Run step 5 - ipyrad automatically creates stats files
data.run("5", force=True, auto=True)
