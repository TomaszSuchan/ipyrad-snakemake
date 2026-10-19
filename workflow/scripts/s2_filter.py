import os
import sys
import ipyrad as ip

# Redirect stdout to log file
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

library_name = snakemake.params.library_name
library_id = snakemake.params.library_id

# Load assembly from step 1
data = ip.load_json(snakemake.params.json_path)

# Set filtering parameters
data.params.datatype = snakemake.params.demux_datatype
data.params.filter_adapters = snakemake.params.filter_adapters
data.params.filter_min_trim_len = snakemake.params.filter_min_trim_len

# Set number of cores for parallel processing
data.ipcluster["cores"] = snakemake.threads

# Run step 2 only
data.run("2", force=True, auto=True)
data.save()
