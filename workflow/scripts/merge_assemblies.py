import os
import ipyrad as ip

proj_config = snakemake.params.project_config
assembly_name = snakemake.params.assembly_name
ipyrad_params = proj_config["ipyrad"]

# Load assemblies
assemblies = []
for json_file in snakemake.params.jsons:
    data = ip.load_json(json_file)

    # Filter to selected samples if sample_list provided
    sample_list_file = snakemake.input.get("sample_list", None)
    if sample_list_file and os.path.exists(sample_list_file):
        with open(sample_list_file) as f:
            selected = set(line.strip() for line in f if line.strip())
        data.samples = {name: sample for name, sample in data.samples.items() if name in selected}
        if not data.samples:
            continue

    assemblies.append(data)

if not assemblies:
    raise Exception("No assemblies to merge!")

# Merge
merged = ip.merge(assembly_name, assemblies)
merged.params.project_dir = os.path.dirname(snakemake.output.json)
merged.params.datatype = ipyrad_params.get("datatype", "pairddrad")
merged.params.restriction_overhang = tuple(ipyrad_params["restriction_overhang"])

# Set reference-based assembly parameters if configured
assembly_method = ipyrad_params.get("assembly_method", "denovo")
if assembly_method != "denovo":
    merged.params.assembly_method = assembly_method

    # Set reference sequence path
    reference_path = None
    if "reference_sequence" in ipyrad_params and ipyrad_params["reference_sequence"]:
        reference_path = ipyrad_params["reference_sequence"]
    elif "reference_assembly" in ipyrad_params and ipyrad_params["reference_assembly"]:
        assembly_id = ipyrad_params["reference_assembly"]
        project_dir = snakemake.wildcards.project_dir
        profile = snakemake.params.project_config.get("profile", "default")
        reference_path = f"{project_dir}/reference/{profile}_{assembly_id}_genomic.fna"

    if reference_path:
        merged.params.reference_sequence = reference_path
        print(f"Using reference-based assembly: {assembly_method}")
        print(f"Reference genome: {reference_path}")

merged.save()

# Marker file used by Snakefile to trigger s3 when merge_assemblies reruns.
with open(snakemake.output.merge_done, "w") as f:
    f.write("merge_complete\n")
