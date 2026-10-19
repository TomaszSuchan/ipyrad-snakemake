# ipyrad-Snakemake Workflow Guide

## Overview

This pipeline supports two workflows per project and per assembly profile:

1. **Optimization workflow** - test multiple clustering thresholds on a subset of samples
2. **Full assembly workflow** - run complete assembly on all samples

Demultiplexing/filtering (ipyrad steps 1-2) is shared and reused across profiles to avoid redundant processing.

## Config Structure

Edit `config/config.yaml` at these keys (Snakemake-style paths):

- `config["projects"][project]["library_prefixes"]`
- `config["projects"][project]["project_dir"]`
- `config["projects"][project]["ipyrad_step12"]` (shared settings for steps 1-2)
- `config["projects"][project]["assemblies"][profile]["ipyrad"]` (settings for steps 3-7)
- `config["projects"][project]["assemblies"][profile]["workflow"]`
- `config["projects"][project]["assemblies"][profile]["optimization"]`
- `config["projects"][project]["assemblies"][profile]["step7_subsets"]`

Use full profile names like `denovo` and `reference`.

### Minimal example

```yaml
projects:
  Pmi_test:
    library_prefixes:
      - "/path/run1/Pmi_RAD02a"
      - "/path/run2/Pmi_RAD02a"
    project_dir: "Pmi_test"

    # Shared step1-2 settings (must be identical for all profiles)
    ipyrad_step12:
      datatype: "pair3rad"
      restriction_overhang: ["TGCAG", "TAA"]
      filter_adapters: 2
      filter_min_trim_len: 35

    assemblies:
      denovo:
        ipyrad:
          datatype: "pairddrad"
          assembly_method: "denovo"
          clust_threshold: 0.85
          mindepth_statistical: 6
          mindepth_majrule: 6
          max_Ns_consens: 0.05
          max_Hs_consens: 0.05
          max_alleles_consens: 2
          min_samples_percent: 50
          output_formats: "pv"
        workflow:
          run_optimization: true
          run_full_assembly: false
        optimization:
          n_samples: 10
          thresholds: [0.80, 0.85]  # Proportions (preferred). Legacy percent ints still supported.
          replicate_suffix: "_repl"
          remove_bottom_percent: 20
          random_seed: 42
          samples: null
          use_replicates: true
        step7_subsets: {}

      reference:
        ipyrad:
          datatype: "pairddrad"
          assembly_method: "reference"
          # reference_sequence: "/path/reference.fasta"
          # or:
          # reference_assembly: "GCA_900128785.1"
          clust_threshold: 0.85
          mindepth_statistical: 6
          mindepth_majrule: 6
          max_Ns_consens: 0.05
          max_Hs_consens: 0.05
          max_alleles_consens: 2
          min_samples_percent: 50
          output_formats: "pv"
        workflow:
          run_optimization: true
          run_full_assembly: false
        optimization:
          n_samples: 10
          thresholds: [0.80, 0.85]
          replicate_suffix: "_repl"
          remove_bottom_percent: 20
          random_seed: 42
          samples: null
          use_replicates: true
        step7_subsets: {}
```

## Library Prefix Requirements

Each `library_prefixes` entry must have:

- `{prefix}_R1_.fastq.gz`
- `{prefix}_R2_.fastq.gz`
- `{prefix}.txt` (barcode file)

If the same library base name appears in multiple directories, each is demultiplexed independently and merged later.

## Output Layout (current)

The pipeline uses a flat project layout (no `assemblies/<profile>/...` directory).
Profile separation is encoded in assembly names and summary filenames.

```text
<project_dir>/
├── demultiplexed/
│   ├── <library_id>/<library_name>.json
│   ├── <library_id>/<library_name>_fastqs/s1_demultiplex_stats.txt
│   ├── <library_id>/<library_name>_edits/s2_rawedit_stats.txt
│   ├── s2_stats_by_sample.tsv
│   └── s2_stats_by_reads.tsv
├── optimization/
│   ├── selected_samples_denovo.txt
│   ├── sample_selection_report_denovo.txt
│   ├── optimization_denovo_80/
│   │   ├── optimization_denovo_80.json
│   │   ├── optimization_denovo_80_clust_.../
│   │   └── optimization_denovo_80_outfiles/
│   ├── optimization_reference_80/
│   │   └── ...
│   ├── optimization_sample_stats_denovo.csv
│   ├── optimization_sample_stats_reference.csv
│   ├── optimization_comparison_denovo.pdf
│   ├── optimization_comparison_reference.pdf
│   ├── supp_fig1_threshold_selection_denovo.pdf
│   └── supp_fig1_threshold_selection_reference.pdf
├── full/
│   ├── <project>_denovo/
│   │   ├── qc_report.txt
│   │   ├── qc/qc_metrics.csv
│   │   └── qc/supp_fig2_assembly_qc.pdf
│   └── <project>_reference/
├── reference/
│   └── reference_<assembly_id>_genomic.fna*
├── logs/
└── benchmarks/
```

## Usage

Run all projects:

```bash
snakemake --use-conda --cores 48 --resources demux_lock=1
```

Run one project:

```bash
snakemake --use-conda --cores 48 --resources demux_lock=1 --config project=Pmi_test
```

Dry-run:

```bash
snakemake -n --use-conda --cores 1 --configfile config/config.yaml
```

## Workflow toggles

Control per profile:

```yaml
assemblies:
  denovo:
    workflow:
      run_optimization: true
      run_full_assembly: false
  reference:
    workflow:
      run_optimization: false
      run_full_assembly: true
```

## Optimization sample selection

Automatic mode:

```yaml
optimization:
  n_samples: 20
  replicate_suffix: "_repl"
  remove_bottom_percent: 20
  random_seed: 42
  samples: null
  use_replicates: true
```

Explicit mode:

```yaml
optimization:
  samples:
    - SampleA
    - SampleB
  use_replicates: true
```

Selection logic:

1. Identify replicate samples using `config["projects"][project]["assemblies"][profile]["optimization"]["replicate_suffix"]`.
2. Keep replicates available for technical concordance checks (`use_replicates: true`).
3. In automatic mode, optionally remove low-read samples using `remove_bottom_percent`.
4. Randomly select up to `n_samples` with fixed `random_seed` for reproducibility.
5. Write selected sample IDs to `optimization/selected_samples_<profile>.txt`.

## Optimization methods (manuscript-ready)

This section describes exactly what the optimization workflow does and what is compared across thresholds.

### Study design and scope

For each profile (`denovo`, `reference`), the pipeline runs a threshold sweep defined in:

- `config["projects"][project]["assemblies"][profile]["optimization"]["thresholds"]`

Thresholds are configured as proportions (for example, `0.85`). Legacy percent-style values (for example, `85`) are still accepted and normalized to proportions internally. Output directory names use percent-like tokens (for example, `optimization_denovo_85`). Shared demultiplexed and filtered reads (steps 1-2) are reused from:

- `<project_dir>/demultiplexed/`

Threshold-specific assemblies are written to:

- `<project_dir>/optimization/optimization_<profile>_<threshold>/`

### Per-threshold assembly procedure

For every threshold in the list, the workflow performs:

1. **Library merge / subset preparation**  
   Builds threshold-specific optimization assemblies from selected samples.
2. **ipyrad step 3 (within-sample clustering)**  
   Uses threshold-specific `clust_threshold`.
3. **ipyrad step 4 (joint estimation)**  
   Estimates heterozygosity and sequencing error parameters.
4. **ipyrad step 5 (consensus calling)**  
   Applies depth, heterozygosity, and allele filters.
5. **ipyrad step 6 (across-sample clustering)**  
   Produces across-sample loci set.
6. **ipyrad step 7 (filter/export)**  
   Exports `vcf` and requested formats.
7. **Hardy-Weinberg scan**  
   Runs HWE test on the exported VCF and stores `.hardy.hwe`.
8. **WHOA genotyping-error scan**  
   Estimates heterozygote miscall rate vs read depth (`qc_whoa`).

### Comparison metrics generated

The `optimization_comparison` rule combines JSON, VCF, HWE, `.loci`, and WHOA outputs across thresholds and writes CSV summaries, a diagnostic PDF, and Supplementary Figure 1.

#### 1) Core ipyrad sample metrics (from JSON stats)

Per-sample metrics are extracted from ipyrad step summaries (excluding replicate-only rows from biological summaries), including:

- Step 3 depth and clustering summaries
- Step 4 heterozygosity and error estimates
- Step 5 consensus filtering and heterozygosity summaries

These are saved in:

- `optimization/optimization_sample_stats_<profile>.csv`

#### 2) Locus retention vs missing-data tradeoff

From VCF per-locus call-rate (`NS`), the pipeline reports loci retained at:

- `>=50%`, `>=70%`, `>=80%` (Paris **r80**), and `>=90%` sample presence
- consecutive-threshold change in r80: `r80_delta` (annotation only; not used as an automatic picker)

Metrics:

- `loci_retained_ge50pct`
- `loci_retained_ge70pct`
- `loci_retained_ge80pct`
- `loci_retained_ge90pct`
- `r80_delta`

#### 3) Per-sample call-rate distribution

Per sample:

- `vcf_call_rate` (called loci / total variant loci)
- `vcf_missing_rate` (1 - call rate)

These are merged into the sample stats CSV and plotted across thresholds.

#### 4) MAF spectrum and singleton burden

From diploid genotype dosages:

- per-locus minor allele frequency: `vcf_maf`
- singleton proportion: `vcf_singleton_proportion`
- doubleton proportion: `vcf_doubleton_proportion`

The pipeline also reports allele-count class totals:

- `vcf_n_biallelic`
- `vcf_n_triallelic`
- `vcf_n_quadriallelic`
- `vcf_n_gt4_allelic`

#### 5) Heterozygosity sanity checks

Per sample, derived from non-missing genotypes:

- observed heterozygosity: `vcf_observed_het`
- expected heterozygosity (approximate): `vcf_expected_het`
- inbreeding coefficient proxy: `vcf_fis = 1 - (H_O / H_E)`

Note: `vcf_expected_het` is calculated as a per-sample QC proxy (`2p(1-p)`) from each sample's mean allele frequency across non-missing loci. This is appropriate for comparing assembly settings and detecting implausible heterozygosity shifts, but it is not a formal population-level expected heterozygosity estimator. Therefore, `vcf_fis` should be interpreted as a heuristic QC indicator rather than a definitive population-genetic parameter.

In addition, HWE deviation burden is computed per sample:

- `hwe_deviated_loci`

and saved in:

- `optimization/optimization_hwe_stats_<profile>.csv`

#### 6) Replicate concordance (technical reproducibility)

Replicate pairing is determined by `replicate_suffix` using the same sample-pair logic as `workflow/scripts/vcf_analyze_replicates.py`.  
For each pair, non-missing genotypes are compared locus-by-locus and discordance rates are summarized as:

- `replicate_discordance_mean`
- `replicate_discordance_sd`
- `replicate_loci_dropout_mean`
- `replicate_allele_dropout_hom_het_mean`
- `replicate_snp_error_mean`

#### 7) Population-structure robustness and relatedness stability (proxy)

Using per-threshold genotype dosage matrices:

- PCA variance explained:
  - `pca_pc1_variance_explained`
  - `pca_pc2_variance_explained`
- Cross-threshold pairwise distance-matrix correlation (assigned to the later threshold):
  - `relatedness_matrix_correlation_prev`

These provide a practical stability check for whether individual-level structure is preserved as clustering thresholds change.

#### 8) McCartney-Melstad distance–missingness correlation

From the same per-threshold dosage matrix and per-sample missingness:

- `genetic_distance_missingness_corr`

A rise at high `clust_threshold` is the oversplitting signal (genetic distance starts tracking missing data).

#### 9) McKinney HDplot (paralog / collapsed-duplicate flag)

Biallelic SNPs with allele depths from FORMAT `AD` or ipyrad `CATG` are scored as heterozygote proportion `H` vs read-ratio deviation `D`. Sites with `H > 0.5` and `|D| > 5` are flagged (SNPfiltR-style cutoffs).

- `hdplot_n_tested`, `hdplot_n_flagged`, `hdplot_prop_flagged`

### VCF-level summary outputs

All VCF-derived long-format metrics are stored in:

- `optimization/optimization_sample_stats_vcf_<profile>.csv`

Included metrics:

- `vcf_ns`, `vcf_dp`, `vcf_maf`
- `vcf_n_fragments`, `vcf_n_variants`
- `vcf_n_biallelic`, `vcf_n_triallelic`, `vcf_n_quadriallelic`, `vcf_n_gt4_allelic`
- `vcf_singleton_proportion`, `vcf_doubleton_proportion`
- `loci_retained_ge50pct`, `loci_retained_ge70pct`, `loci_retained_ge80pct`, `loci_retained_ge90pct`, `r80_delta`
- `replicate_discordance_mean`, `replicate_discordance_sd`
- `pca_pc1_variance_explained`, `pca_pc2_variance_explained`
- `relatedness_matrix_correlation_prev`
- `variant_jaccard_prev`
- `hwe_prop_rejected`, `hwe_prop_deficit`, `hwe_prop_excess`
- `mean_pairwise_shared_loci`
- `whoa_mean_miscall_rate`, `whoa_mean_m` (with `mean_dp`, `lo95`, `hi95`)
- `hdplot_prop_flagged`, `hdplot_n_flagged`, `hdplot_n_tested`
- `genetic_distance_missingness_corr`

### Comparison figure (diagnostic dump)

The analyst dump:

- `optimization/optimization_comparison_<profile>.pdf`

contains multi-panel distributions and summary plots across thresholds for all available metrics listed above.

### Supplementary Figure 1 — threshold selection

Paper-ready figure (one PDF per profile):

- `optimization/supp_fig1_threshold_selection_<profile>.pdf`

Panels:

- **A** SNPs and RAD fragments vs threshold
- **B** loci retained at ≥50/70/80/90% sample presence; dotted line at the r80 maximum, with last-step Δ r80 annotated
- **C** consecutive-threshold variant-set Jaccard
- **D** replicate discordance, SNP error, and Mastretta allele/locus dropout
- **E** step 4 error and heterozygosity estimates
- **F** singleton, multiallelic, and HDplot-flagged proportions
- **G** per-sample observed heterozygosity
- **H** global HWE proportion rejected (deficit vs excess)
- **I** relatedness-matrix correlation with the previous threshold
- **J** WHOA heterozygote-miscall rate vs depth (one curve per threshold)
- **K** mean pairwise shared loci
- **L** McCartney-Melstad genetic distance–missingness correlation

WHOA is run per optimization assembly via `qc_whoa`. Optimization VCFs use `qc.whoa_min_bin_optimization` (default 200) so small sample subsets still fill depth bins. High `mean_m` that does not fall at high depth is the clustering/paralog signal; high `mean_m` only at low depth is a min-depth issue.

### Supplementary Figure 2 — final assembly QC

When `run_full_assembly: true`, the DAG runs WHOA, HWE, replicate analysis, and shared-loci QC and writes:

- `full/<assembly>/qc_report.txt` — numeric appendix (text)
- `full/<assembly>/qc/qc_metrics.csv` — long-format metrics for the figure
- `full/<assembly>/qc/supp_fig2_assembly_qc.pdf` — paper figure
- `full/<assembly>/<assembly>_outfiles/qc/<assembly>.whoa_plot.pdf`
- `full/<assembly>/<assembly>_outfiles/qc/<assembly>.hardy_report.txt`
- `full/<assembly>/<assembly>_outfiles/qc/<assembly>.replicate_analysis.txt`
- `full/<assembly>/<assembly>_outfiles/qc/<assembly>.shared_loci_heatmap.pdf`

Fig. S2 panels:

- **A** per-sample call rate
- **B** per-sample observed heterozygosity
- **C** locus depth (VCF DP)
- **D** SNPs per RAD locus
- **E** locus length
- **F** MAF spectrum
- **G** missingness per locus
- **H** replicate discordance per pair
- **I** HWE P-value histogram
- **J** shared-loci clustered heatmap (second page)
- **K** WHOA miscall rate vs depth
- **L** PCA PC1–PC2 (unlabeled outlier check)
- **M** HDplot (`H` vs `D`; sites with `H > 0.5` and `|D| > 5` flagged)

`qc.replicate_suffix` defaults to `optimization.replicate_suffix` so replicate pairing is the same in both workflows.
