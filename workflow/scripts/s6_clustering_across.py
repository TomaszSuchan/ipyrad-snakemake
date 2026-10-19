import os
import sys
import ipyrad as ip

# Redirect stdout to log file
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

data = ip.load_json(snakemake.params.json_path)
data.name = os.path.basename(os.path.dirname(snakemake.params.json_path))
data.params.project_dir = os.path.dirname(snakemake.params.json_path)

# Set number of cores for parallel processing
data.ipcluster["cores"] = snakemake.threads

# Run step 6
data.run("6", force=True, auto=True)
data.save()

# Create dummy stats file if ipyrad doesn't create one
stats_path = snakemake.output.stats
if not os.path.exists(stats_path):
    os.makedirs(os.path.dirname(stats_path), exist_ok=True)
    with open(stats_path, 'w') as f:
        f.write(f"Step 6 completed: {data.name}\n")
        f.write(data.stats.to_string() + "\n")
