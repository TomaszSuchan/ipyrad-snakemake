## Quick orientation for AI contributors

This repository is a Snakemake workflow that orchestrates ipyrad assemblies for multiple independent "projects" (assemblies). The goal of these instructions is to help an AI coding agent be immediately productive here by highlighting the important structure, conventions, and common tasks.

- Primary entry points:
  - `workflow/Snakefile` — core orchestration. Read this first.
  - `config/config.yaml` — per-project configuration, resources, and defaults.
  - `workflow/scripts/*.py` — implementation of each Snakemake rule (s1..s7, select_optimization_samples, final_report).
  - `README.md` — high-level usage examples and run commands.

- Run pattern and examples:
  - Typical run (uses conda environments defined under `workflow/envs`):
    - `snakemake --use-conda --cores 48 --resources demux_lock=1`
  - Run a single project by passing config on the command-line:
    - `snakemake --use-conda --cores 48 --config project=Pmi_test`
  - The Snakefile defines `FINAL_OUTPUTS` per-project; avoid changing that layout without updating rules that build `JSON_TO_LIBRARY`.

- Important config conventions (documented in `config/config.yaml`):
  - Projects live under the `projects:` key. Each project must define `library_prefixes` and `project_dir`.
  - A library prefix corresponds to three files:
    - `{prefix}_R1_.fastq.gz` (forward reads)
    - `{prefix}_R2_.fastq.gz` (reverse reads)
    - `{prefix}.txt` (barcode file)
  - Example project key: `Pmi_test` is defined in the default `config.yaml` — use it for local testing.

- File/path naming and mapping rules the code relies on:
  - Snakefile builds a safe library id by replacing `/` with `__` in the prefix (safe_id = prefix.replace("/", "__")). That safe id is used for demultiplexed output directories.
  - Demultiplexed JSONs are expected at: `{project_dir}/demultiplexed/{safe_id}/{lib_basename}.json` (see `JSON_TO_LIBRARY` mapping).
  - Optimization assemblies are named like `clust_85` for threshold 0.85 — many rules parse these names (e.g., `threshold = int(assembly_name.split("_")[-1]) / 100.0`).

- Workflow toggles and hotspots:
  - `workflow.run_optimization` and `workflow.run_full_assembly` control the two main modes.
  - Heavy/long steps (CPU/memory/RT-sensitive): `s3_clustering_within` and `s6_clustering_across` (see `config/resources` for default thread counts and mem_per_cpu). Treat resource changes carefully.
  - `scripts/*` are executed via Snakemake `script:` directives and expect specific input/output paths — editing their I/O contract requires updating corresponding rule blocks.

- Conda and environments:
  - The Snakefile references `envs/ipyrad.yaml` for rule environments. If you modify dependencies, update that YAML and test locally with `--use-conda`.

- Logging, outputs and where to look when things fail:
  - Rule logs are written under `logs/` with predictable names like `logs/step2_{json_path}.log`.
  - Final assembled outputs are placed under each project's `project_dir` in subfolders `optimization/` or `full/` depending on the mode.

- Tests / validation approach for changes:
  - Use the included `Pmi_test` project in `config/config.yaml` as a minimal integration example (it points to local `~/Pmi-test/*` paths in the example config). For unit-style checks, prefer running the Snakefile targeting a small project and a single rule.
  - Re-run `snakemake --use-conda --cores 2 --config project=Pmi_test` to validate end-to-end behavior for small edits.

- Patterns and pitfalls for contributors (do not change without care):
  - Wildcard constraints: the Snakefile builds large regex constraints from `JSON_TO_LIBRARY` keys — adding new wildcard-based rules may require regenerating or extending these constraints.
  - Library merging logic: libraries that share the same basename are intentionally merged (see prints in Snakefile that group by `lib_info['name']`). Changing that behavior affects downstream merge and ipyrad steps.
  - Don't refactor `JSON_TO_LIBRARY` shape lightly — many lambda closures in rule params depend on its keys and structure.

- Example snippets to reference while coding:
  - Path check and safe id creation (Snakefile):
    - `safe_id = prefix.replace("/", "__")`
    - `r1_path = f"{prefix}_R1_.fastq.gz"`
  - Common run command used in README: `snakemake --use-conda --cores 48 --resources demux_lock=1`

If any section is unclear or you want more examples (e.g., how to add a new rule, how to add a project entry, or how to extend `envs/ipyrad.yaml`), tell me which part to expand and I will update this file.

---
Files to inspect for more details: `workflow/Snakefile`, `config/config.yaml`, `workflow/scripts/*.py`, `workflow/envs/ipyrad.yaml`, `README.md`.
