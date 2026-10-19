import os
import sys
import ipyrad as ip

# Redirect stdout to log file
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

library_name = snakemake.params.library_name
library_id = snakemake.params.library_id
raw_dir = snakemake.params.raw_dir

# Setup
library_dir = f"{snakemake.params.project_dir}/demultiplexed/{library_id}"
os.makedirs(library_dir, exist_ok=True)

data = ip.Assembly(library_name)
data.params.project_dir = library_dir
data.params.raw_fastq_path = f"{raw_dir}/{library_name}_R*_.fastq.gz"
data.params.barcodes_path = snakemake.input.barcode
data.params.datatype = snakemake.params.demux_datatype
data.params.restriction_overhang = tuple(snakemake.params.restriction_overhang)

# Set number of cores for parallel processing
data.ipcluster["cores"] = snakemake.threads

# Run step 1 only
data.run("1", force=True, auto=True)
data.save()
