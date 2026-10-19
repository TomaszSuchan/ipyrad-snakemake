#!/usr/bin/env python3
"""Generate diagnostic and supplementary plots from extracted statistics."""
import argparse
import os
import sys
from io import BytesIO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import joypy
    HAS_JOYPY = True
except ImportError:
    HAS_JOYPY = False

try:
    from qc_common import HDPLOT_D_CUTOFF, HDPLOT_H_CUTOFF, hdplot_is_paralog
except ImportError:
    HDPLOT_H_CUTOFF = 0.5
    HDPLOT_D_CUTOFF = 5.0

    def hdplot_is_paralog(h, d, h_cut=HDPLOT_H_CUTOFF, d_cut=HDPLOT_D_CUTOFF):
        if h is None or d is None or np.isnan(h) or np.isnan(d):
            return False
        return h > h_cut and abs(d) > d_cut


def threshold_label(value):
    """Format threshold values consistently as percentages."""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return str(value)
    pct = x * 100.0 if x <= 1.0 else x
    if abs(pct - round(pct)) < 1e-9:
        return f"{int(round(pct))}%"
    return f"{pct:.2f}%"


def threshold_numeric(value):
    """Return clustering threshold as a proportion in (0, 1]."""
    x = float(value)
    return x if x <= 1.0 else x / 100.0


def add_candidate(candidates, candidate_path, input_csv):
    """Add non-empty candidate paths excluding the input CSV itself."""
    if not candidate_path:
        return
    if os.path.normpath(candidate_path) == os.path.normpath(input_csv):
        return
    if candidate_path not in candidates:
        candidates.append(candidate_path)


def load_optional_csv(path):
    if path and os.path.exists(path):
        return pd.read_csv(path)
    return None


def first_by_threshold(vcf_df, metric):
    if vcf_df is None or vcf_df.empty or "metric" not in vcf_df.columns:
        return pd.DataFrame(columns=["threshold", "value"])
    sub = vcf_df[vcf_df["metric"] == metric].copy()
    if sub.empty:
        return pd.DataFrame(columns=["threshold", "value"])
    out = sub.groupby("threshold", as_index=False)["value"].first()
    out["threshold"] = out["threshold"].map(threshold_numeric)
    return out.sort_values("threshold")


def label_panel(ax, letter):
    ax.text(
        -0.12, 1.08, letter, transform=ax.transAxes,
        fontsize=12, fontweight="bold", va="top", ha="left"
    )


def style_ax(ax, ylabel, title=None):
    ax.set_ylabel(ylabel, fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, pad=4)
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
    ax.tick_params(labelsize=8)


def plot_empty(ax, message="No data"):
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, fontsize=9, color="0.4")
    ax.set_xticks([])
    ax.set_yticks([])


def plot_diagnostic(df, vcf_df, hwe_df, output_pdf):
    """Kitchen-sink ridgeline/bar grid used as the analyst dump."""
    if hwe_df is not None and "hwe_deviated_loci" in hwe_df.columns:
        merge_cols = ["threshold", "assembly_name", "sample"]
        available_merge_cols = [col for col in merge_cols if col in df.columns and col in hwe_df.columns]
        if available_merge_cols:
            df = df.merge(
                hwe_df[available_merge_cols + ["hwe_deviated_loci"]],
                on=available_merge_cols,
                how="left",
            )

    if len(df) == 0:
        print("No data to plot.", file=sys.stderr)
        sys.exit(1)

    metrics = {
        "s3": [
            "s3_merged_pairs", "s3_clusters_total", "s3_clusters_hidepth",
            "s3_avg_depth_total", "s3_avg_depth_mj", "s3_avg_depth_stat",
            "s3_sd_depth_total", "s3_sd_depth_mj", "s3_sd_depth_stat",
            "s3_filtered_bad_align",
        ],
        "s4": ["s4_hetero_est", "s4_error_est"],
        "s5": [
            "s5_clusters_total", "s5_filtered_by_depth", "s5_filtered_by_maxH",
            "s5_filtered_by_maxAlleles", "s5_filtered_by_maxN", "s5_reads_consens",
            "s5_nsites", "s5_nhetero", "s5_heterozygosity",
        ],
        "hwe": ["hwe_deviated_loci"],
        "sample_qc": [
            "vcf_call_rate", "vcf_missing_rate", "vcf_observed_het",
            "vcf_expected_het", "vcf_fis", "singleton_proportion",
        ],
        "vcf": [
            "vcf_ns", "vcf_dp", "vcf_maf", "vcf_n_fragments", "vcf_n_variants",
            "vcf_n_biallelic", "vcf_n_triallelic", "vcf_n_quadriallelic", "vcf_n_gt4_allelic",
            "vcf_singleton_proportion", "vcf_doubleton_proportion",
            "loci_retained_ge50pct", "loci_retained_ge70pct", "loci_retained_ge80pct",
            "loci_retained_ge90pct", "r80_delta",
            "replicate_discordance_mean", "replicate_discordance_sd",
            "replicate_loci_dropout_mean", "replicate_allele_dropout_hom_het_mean",
            "replicate_allele_dropout_hom_hom_mean", "replicate_snp_error_mean",
            "pca_pc1_variance_explained", "pca_pc2_variance_explained",
            "relatedness_matrix_correlation_prev",
            "hwe_prop_rejected", "variant_jaccard_prev",
            "mean_pairwise_shared_loci", "whoa_mean_miscall_rate",
            "hdplot_prop_flagged", "hdplot_n_flagged", "hdplot_n_tested",
            "genetic_distance_missingness_corr",
        ],
    }

    all_metrics = []
    for step_metrics in metrics.values():
        all_metrics.extend(step_metrics)
    available_metrics = [m for m in all_metrics if m in df.columns]
    vcf_metrics_in_data = set()

    if vcf_df is not None and len(vcf_df) > 0:
        vcf_metrics_in_data = set(vcf_df["metric"].unique())
        for vcf_metric in metrics["vcf"]:
            if vcf_metric in vcf_metrics_in_data and vcf_metric not in available_metrics:
                available_metrics.append(vcf_metric)

    if len(available_metrics) == 0:
        print("No matching metrics found in data.", file=sys.stderr)
        sys.exit(1)

    df = df.copy()
    df["threshold_str"] = df["threshold"].apply(threshold_label)

    n_metrics = len(available_metrics)
    n_cols = 4
    n_rows = (n_metrics + n_cols - 1) // n_cols
    fig = plt.figure(figsize=(20, 5 * n_rows))
    actual_idx = 0
    base_color = "#1f77b4"
    min_alpha = 0.4
    max_alpha = 0.8

    summary_stats = {
        "vcf_n_fragments", "vcf_n_variants", "vcf_n_biallelic",
        "vcf_n_triallelic", "vcf_n_quadriallelic", "vcf_n_gt4_allelic",
        "vcf_singleton_proportion", "vcf_doubleton_proportion",
        "loci_retained_ge50pct", "loci_retained_ge70pct", "loci_retained_ge80pct",
        "loci_retained_ge90pct", "r80_delta",
        "replicate_discordance_mean", "replicate_discordance_sd",
        "replicate_loci_dropout_mean", "replicate_allele_dropout_hom_het_mean",
        "replicate_allele_dropout_hom_hom_mean", "replicate_snp_error_mean",
        "pca_pc1_variance_explained", "pca_pc2_variance_explained",
        "relatedness_matrix_correlation_prev", "hwe_prop_rejected",
        "variant_jaccard_prev", "mean_pairwise_shared_loci", "whoa_mean_miscall_rate",
        "hdplot_prop_flagged", "hdplot_n_flagged", "hdplot_n_tested",
        "genetic_distance_missingness_corr",
    }

    for metric in available_metrics:
        is_summary_stat = False
        if vcf_df is not None and metric in vcf_metrics_in_data:
            metric_data = vcf_df[vcf_df["metric"] == metric].copy()
            if len(metric_data) == 0:
                continue
            is_summary_stat = metric in summary_stats
            if is_summary_stat:
                metric_data = metric_data.groupby("threshold")["value"].first().reset_index()
                metric_data["threshold_str"] = metric_data["threshold"].apply(threshold_label)
                plot_data = metric_data[["threshold_str", "value"]].rename(columns={"value": metric})
            else:
                metric_data["threshold_str"] = metric_data["threshold"].apply(threshold_label)
                plot_data = metric_data[["threshold_str", "value"]].rename(columns={"value": metric}).dropna()
        else:
            plot_data = df[["threshold_str", metric]].dropna()

        if len(plot_data) == 0:
            continue

        actual_idx += 1
        ax = plt.subplot(n_rows, n_cols, actual_idx)
        clean_name = metric.replace("_", " ").title()
        clean_name = clean_name.replace("S3 ", "Step 3: ")
        clean_name = clean_name.replace("S4 ", "Step 4: ")
        clean_name = clean_name.replace("S5 ", "Step 5: ")
        clean_name = clean_name.replace("Vcf ", "VCF: ")
        thresholds = sorted(plot_data["threshold_str"].unique())
        n_thresholds = len(thresholds)

        if is_summary_stat:
            plot_data_sorted = plot_data.sort_values("threshold_str")
            ax.bar(
                range(len(plot_data_sorted)), plot_data_sorted[metric],
                color=base_color, alpha=0.7, edgecolor="black", linewidth=0.5,
            )
            ax.set_xticks(range(len(plot_data_sorted)))
            ax.set_xticklabels(plot_data_sorted["threshold_str"], fontsize=8, rotation=0)
            ax.set_ylabel("Value", fontsize=9)
            ax.set_title(clean_name, fontsize=10, fontweight="bold", pad=8)
            ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
            continue

        if metric == "vcf_dp":
            plot_data_log = plot_data[plot_data[metric] > 0].copy()
            if len(plot_data_log) > 0:
                plot_data_log[metric] = np.log10(plot_data_log[metric])
                plot_data = plot_data_log
                clean_name = clean_name + " (log10)"

        if not HAS_JOYPY:
            ax.hist(plot_data[metric].dropna(), bins=30, color=base_color, alpha=0.7)
            ax.set_title(clean_name, fontsize=10, fontweight="bold", pad=8)
            continue

        fig_temp, axes_joy = joypy.joyplot(
            plot_data,
            by="threshold_str",
            column=metric,
            figsize=(8, 4),
            overlap=1.2,
            kind="kde",
            fade=False,
            linewidth=1.5,
        )
        for i, ax_joy in enumerate(axes_joy):
            alpha = min_alpha + (i / max(n_thresholds - 1, 1)) * (max_alpha - min_alpha)
            for collection in ax_joy.collections:
                collection.set_color(base_color)
                collection.set_alpha(alpha)
            for line in ax_joy.lines:
                line.set_color(base_color)
                line.set_alpha(alpha)

        buf = BytesIO()
        fig_temp.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0)
        img = plt.imread(buf)
        buf.close()
        ax.imshow(img, aspect="auto", extent=[0, 1, 0, 1])
        ax.axis("off")
        ax.set_title(clean_name, fontsize=10, fontweight="bold", pad=8)
        plt.close(fig_temp)

    plt.subplots_adjust(left=0.04, right=0.99, top=0.98, bottom=0.02, hspace=0.15, wspace=0.15)
    fig.savefig(output_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Plots saved to {output_pdf}")
    print(f"Plotted {len(available_metrics)} metrics across {len(df['threshold'].unique())} thresholds on one page")


def plot_threshold_selection_supp(df, vcf_df, hwe_df, replicate_df, output_pdf):
    """Journal-style Fig. S1: clustering-threshold selection."""
    fig, axes = plt.subplots(3, 4, figsize=(14.5, 10.5))
    axes = axes.ravel()
    cmap = plt.cm.viridis

    # A: SNPs and RAD fragments
    ax = axes[0]
    snps = first_by_threshold(vcf_df, "vcf_n_variants")
    frags = first_by_threshold(vcf_df, "vcf_n_fragments")
    if not snps.empty or not frags.empty:
        if not frags.empty:
            ax.plot(frags["threshold"], frags["value"], "o-", label="RAD fragments")
        if not snps.empty:
            ax.plot(snps["threshold"], snps["value"], "s--", label="SNPs")
        ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Count")
    else:
        plot_empty(ax)
    label_panel(ax, "A")

    # B: retention including Paris r80
    ax = axes[1]
    plotted = False
    for metric, label in [
        ("loci_retained_ge50pct", ">=50%"),
        ("loci_retained_ge70pct", ">=70%"),
        ("loci_retained_ge80pct", ">=80% (r80)"),
        ("loci_retained_ge90pct", ">=90%"),
    ]:
        series = first_by_threshold(vcf_df, metric)
        if not series.empty:
            ax.plot(series["threshold"], series["value"], "o-", label=label)
            plotted = True
    r80 = first_by_threshold(vcf_df, "loci_retained_ge80pct")
    if not r80.empty and r80["value"].notna().any():
        idx = r80["value"].idxmax()
        t_max = r80.loc[idx, "threshold"]
        ax.axvline(t_max, color="0.35", linestyle=":", linewidth=1.2, zorder=0)
        r80_sorted = r80.sort_values("threshold")
        delta = first_by_threshold(vcf_df, "r80_delta")
        note = f"max r80 at {threshold_label(t_max)}"
        if not delta.empty and delta["value"].notna().any():
            last = delta.sort_values("threshold").iloc[-1]
            note += f"\nΔ r80 last step: {last['value']:.0f}"
        elif len(r80_sorted) >= 2:
            last_delta = float(r80_sorted["value"].iloc[-1] - r80_sorted["value"].iloc[-2])
            note += f"\nΔ r80 last step: {last_delta:.0f}"
        ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=7, va="bottom")
    if plotted:
        ax.legend(fontsize=7, frameon=False, title="Presence")
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Loci retained")
    else:
        plot_empty(ax)
    label_panel(ax, "B")

    # C: Jaccard
    ax = axes[2]
    jac = first_by_threshold(vcf_df, "variant_jaccard_prev")
    if not jac.empty:
        ax.plot(jac["threshold"], jac["value"], "o-", color="#d62728")
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Jaccard vs previous")
    else:
        plot_empty(ax, "Need >=2 thresholds")
    label_panel(ax, "C")

    # D: replicate discordance / SNP error / Mastretta allele dropout
    ax = axes[3]
    disc = first_by_threshold(vcf_df, "replicate_discordance_mean")
    snp_err = first_by_threshold(vcf_df, "replicate_snp_error_mean")
    allele_do = first_by_threshold(vcf_df, "replicate_allele_dropout_hom_het_mean")
    loci_do = first_by_threshold(vcf_df, "replicate_loci_dropout_mean")
    plotted = False
    if replicate_df is not None and not replicate_df.empty and "replicate_discordance" in replicate_df.columns:
        r = replicate_df.copy()
        r["threshold"] = r["threshold"].map(threshold_numeric)
        ax.scatter(r["threshold"], r["replicate_discordance"], s=18, alpha=0.45, color="0.5", label="pairs")
        plotted = True
    if not disc.empty:
        ax.plot(disc["threshold"], disc["value"], "o-", color="#1f77b4", label="mean discordance")
        plotted = True
    if not snp_err.empty:
        ax.plot(snp_err["threshold"], snp_err["value"], "s--", color="#ff7f0e", label="mean SNP error")
        plotted = True
    if not allele_do.empty:
        ax.plot(
            allele_do["threshold"], allele_do["value"], "^-",
            color="#d62728", label="allele dropout",
        )
        plotted = True
    if not loci_do.empty:
        ax.plot(
            loci_do["threshold"], loci_do["value"], "x-.",
            color="#8c564b", label="locus dropout",
        )
        plotted = True
    if plotted:
        ax.legend(fontsize=6, frameon=False)
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Rate")
    else:
        plot_empty(ax, "No replicate pairs")
    label_panel(ax, "D")

    # E: s4 estimates
    ax = axes[4]
    plotted = False
    if df is not None and not df.empty:
        tmp = df.copy()
        tmp["threshold"] = tmp["threshold"].map(threshold_numeric)
        for col, label in [("s4_error_est", "error"), ("s4_hetero_est", "heterozygosity")]:
            if col in tmp.columns:
                g = tmp.groupby("threshold")[col].mean().reset_index()
                if g[col].notna().any():
                    ax.plot(g["threshold"], g[col], "o-", label=label)
                    plotted = True
    if plotted:
        ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Step 4 estimate")
    else:
        plot_empty(ax)
    label_panel(ax, "E")

    # F: singleton + multiallelic + HDplot paralog flag
    ax = axes[5]
    sing = first_by_threshold(vcf_df, "vcf_singleton_proportion")
    n_var = first_by_threshold(vcf_df, "vcf_n_variants")
    n_bi = first_by_threshold(vcf_df, "vcf_n_biallelic")
    hdplot = first_by_threshold(vcf_df, "hdplot_prop_flagged")
    plotted = False
    if not sing.empty:
        ax.plot(sing["threshold"], sing["value"], "o-", label="singletons")
        plotted = True
    if not n_var.empty and not n_bi.empty:
        merged = n_var.merge(n_bi, on="threshold", suffixes=("_n", "_bi"))
        merged["multi"] = 1.0 - (merged["value_bi"] / merged["value_n"].replace(0, np.nan))
        ax.plot(merged["threshold"], merged["multi"], "s--", label="multiallelic")
        plotted = True
    if not hdplot.empty:
        ax.plot(hdplot["threshold"], hdplot["value"], "^-", color="#d62728", label="HDplot flagged")
        plotted = True
    if plotted:
        ax.legend(fontsize=7, frameon=False)
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Proportion")
    else:
        plot_empty(ax)
    label_panel(ax, "F")

    # G: observed heterozygosity distributions
    ax = axes[6]
    if df is not None and not df.empty and "vcf_observed_het" in df.columns:
        tmp = df[["threshold", "vcf_observed_het"]].dropna().copy()
        tmp["threshold"] = tmp["threshold"].map(threshold_numeric)
        groups = list(tmp.groupby("threshold"))
        if groups:
            data = [g["vcf_observed_het"].values for _, g in groups]
            labels = [threshold_label(t) for t, _ in groups]
            ax.boxplot(data, showfliers=False)
            ax.set_xticks(range(1, len(labels) + 1))
            ax.set_xticklabels(labels)
            ax.set_xlabel("Clustering threshold", fontsize=9)
            style_ax(ax, "Observed heterozygosity")
        else:
            plot_empty(ax)
    else:
        plot_empty(ax)
    label_panel(ax, "G")

    # H: global HWE
    ax = axes[7]
    plotted = False
    for metric, label in [
        ("hwe_prop_rejected", "P < 0.05"),
        ("hwe_prop_deficit", "het deficit"),
        ("hwe_prop_excess", "het excess"),
    ]:
        series = first_by_threshold(vcf_df, metric)
        if not series.empty:
            ax.plot(series["threshold"], series["value"], "o-", label=label)
            plotted = True
    if plotted:
        ax.legend(fontsize=7, frameon=False)
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Proportion of loci")
    else:
        plot_empty(ax)
    label_panel(ax, "H")

    # I: relatedness stability
    ax = axes[8]
    rel = first_by_threshold(vcf_df, "relatedness_matrix_correlation_prev")
    if not rel.empty:
        ax.plot(rel["threshold"], rel["value"], "o-", color="#2ca02c")
        ax.set_ylim(min(0.0, float(rel["value"].min()) - 0.05), 1.05)
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Distance-matrix r")
    else:
        plot_empty(ax, "Need >=2 thresholds")
    label_panel(ax, "I")

    # J: WHOA mean_m vs depth
    ax = axes[9]
    if vcf_df is not None and not vcf_df.empty and "whoa_mean_m" in set(vcf_df["metric"].unique()):
        whoa = vcf_df[vcf_df["metric"] == "whoa_mean_m"].copy()
        whoa["threshold"] = whoa["threshold"].map(threshold_numeric)
        thresholds = sorted(whoa["threshold"].unique())
        colors = cmap(np.linspace(0.15, 0.9, max(len(thresholds), 1)))
        for color, t in zip(colors, thresholds):
            g = whoa[whoa["threshold"] == t].sort_values("mean_dp" if "mean_dp" in whoa.columns else "value")
            x = g["mean_dp"] if "mean_dp" in g.columns else np.arange(len(g))
            y = g["value"]
            if "lo95" in g.columns and "hi95" in g.columns and g["lo95"].notna().any():
                yerr = np.vstack([
                    np.clip(y - g["lo95"], 0, None),
                    np.clip(g["hi95"] - y, 0, None),
                ])
                ax.errorbar(x, y, yerr=yerr, fmt="o-", color=color, label=threshold_label(t), capsize=2, linewidth=1)
            else:
                ax.plot(x, y, "o-", color=color, label=threshold_label(t))
        ax.legend(fontsize=7, frameon=False, title="Threshold")
        ax.set_xlabel("Mean depth (bin)", fontsize=9)
        style_ax(ax, "Het. miscall rate (m)")
    else:
        plot_empty(ax, "No WHOA results")
    label_panel(ax, "J")

    # Extra useful panel: mean pairwise shared loci
    ax = axes[10]
    shared = first_by_threshold(vcf_df, "mean_pairwise_shared_loci")
    if not shared.empty:
        ax.plot(shared["threshold"], shared["value"], "o-", color="#9467bd")
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Mean pairwise shared loci")
    else:
        plot_empty(ax)
    label_panel(ax, "K")

    # L: McCartney-Melstad genetic distance vs missingness
    ax = axes[11]
    dist_miss = first_by_threshold(vcf_df, "genetic_distance_missingness_corr")
    if not dist_miss.empty:
        ax.plot(dist_miss["threshold"], dist_miss["value"], "o-", color="#17becf")
        ax.axhline(0.0, color="0.5", linestyle="--", linewidth=0.8)
        ax.set_xlabel("Clustering threshold", fontsize=9)
        style_ax(ax, "Distance–missingness r")
    else:
        plot_empty(ax, "Need >=3 samples")
    label_panel(ax, "L")
    fig.suptitle("Supplementary Figure 1. Clustering threshold selection", fontsize=13, y=0.995)
    fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.97])
    fig.savefig(output_pdf, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"Supplementary threshold-selection figure saved to {output_pdf}")


def plot_whoa_on_ax(ax, whoa_summary):
    if not whoa_summary or not whoa_summary.get("depth_bins"):
        plot_empty(ax, "No WHOA results")
        return
    bins = whoa_summary["depth_bins"]
    x = [b.get("mean_dp") if b.get("mean_dp") is not None else i for i, b in enumerate(bins)]
    y = [b["mean_m"] for b in bins]
    lo = [b.get("lo95", np.nan) for b in bins]
    hi = [b.get("hi95", np.nan) for b in bins]
    if np.isfinite(np.array(lo, dtype=float)).any():
        yerr = np.vstack([
            np.clip(np.array(y) - np.array(lo, dtype=float), 0, None),
            np.clip(np.array(hi, dtype=float) - np.array(y), 0, None),
        ])
        ax.errorbar(x, y, yerr=yerr, fmt="o-", capsize=3, color="#1f77b4")
    else:
        ax.plot(x, y, "o-")
    ax.set_xlabel("Mean depth (bin)", fontsize=9)
    style_ax(ax, "Het. miscall rate (m)")


def plot_full_assembly_qc(output_pdf, data):
    """Journal-style Fig. S2 for a finished assembly."""
    sample_df = data.get("sample_df")
    replicate_df = data.get("replicate_df")
    pca = data.get("pca")
    whoa = data.get("whoa")
    hwe = data.get("hwe")
    shared_samples = data.get("shared_samples")
    shared_matrix = data.get("shared_matrix")

    with PdfPages(output_pdf) as pdf:
        fig, axes = plt.subplots(3, 4, figsize=(14.5, 10.5))
        axes = axes.ravel()

        def sample_values(column):
            if sample_df is None or sample_df.empty or column not in sample_df.columns:
                return np.array([])
            return sample_df[column].dropna().to_numpy()

        # A call rate
        ax = axes[0]
        vals = sample_values("call_rate")
        if vals.size:
            ax.hist(vals, bins=20, color="#1f77b4", edgecolor="white")
            style_ax(ax, "Samples", None)
            ax.set_xlabel("Call rate", fontsize=9)
        else:
            plot_empty(ax)
        label_panel(ax, "A")

        # B heterozygosity
        ax = axes[1]
        vals = sample_values("observed_het")
        if vals.size:
            ax.hist(vals, bins=20, color="#2ca02c", edgecolor="white")
            ax.set_xlabel("Observed heterozygosity", fontsize=9)
            style_ax(ax, "Samples")
        else:
            plot_empty(ax)
        label_panel(ax, "B")

        # C depth
        ax = axes[2]
        dp = np.asarray(data.get("dp_values") or [], dtype=float)
        dp = dp[np.isfinite(dp) & (dp > 0)]
        if dp.size:
            ax.hist(np.log10(dp), bins=30, color="#ff7f0e", edgecolor="white")
            ax.set_xlabel("log10 locus DP", fontsize=9)
            style_ax(ax, "Loci")
        else:
            plot_empty(ax)
        label_panel(ax, "C")

        # D SNPs per locus
        ax = axes[3]
        vpl = np.asarray(data.get("snps_per_locus") or [], dtype=float)
        vpl = vpl[np.isfinite(vpl)]
        if vpl.size:
            ax.hist(vpl, bins=min(30, max(8, int(np.sqrt(vpl.size)))), color="#9467bd", edgecolor="white")
            ax.set_xlabel("SNPs per RAD locus", fontsize=9)
            style_ax(ax, "Loci")
        else:
            plot_empty(ax)
        label_panel(ax, "D")

        # E locus length
        ax = axes[4]
        lengths = np.asarray(data.get("locus_lengths") or [], dtype=float)
        lengths = lengths[np.isfinite(lengths)]
        if lengths.size:
            ax.hist(lengths, bins=20, color="#8c564b", edgecolor="white")
            ax.set_xlabel("Locus length (bp)", fontsize=9)
            style_ax(ax, "Loci")
        else:
            plot_empty(ax)
        label_panel(ax, "E")

        # F MAF
        ax = axes[5]
        maf = np.asarray(data.get("maf_values") or [], dtype=float)
        maf = maf[np.isfinite(maf)]
        if maf.size:
            ax.hist(maf, bins=25, color="#e377c2", edgecolor="white")
            ax.set_xlabel("Minor allele frequency", fontsize=9)
            style_ax(ax, "SNPs")
        else:
            plot_empty(ax)
        label_panel(ax, "F")

        # G missingness
        ax = axes[6]
        miss = np.asarray(data.get("missing_per_locus") or [], dtype=float)
        miss = miss[np.isfinite(miss)]
        if miss.size:
            ax.hist(miss, bins=20, color="#7f7f7f", edgecolor="white")
            ax.set_xlabel("Missingness per locus", fontsize=9)
            style_ax(ax, "Loci")
        else:
            plot_empty(ax)
        label_panel(ax, "G")

        # H replicates
        ax = axes[7]
        if replicate_df is not None and not replicate_df.empty and "replicate_discordance" in replicate_df.columns:
            r = replicate_df.dropna(subset=["replicate_discordance"]).copy()
            if not r.empty:
                labels = [f"{row.sample}" for row in r.itertuples()]
                ax.bar(range(len(r)), r["replicate_discordance"], color="#d62728", alpha=0.8)
                ax.set_xticks(range(len(r)))
                ax.set_xticklabels(labels, rotation=90, fontsize=6)
                style_ax(ax, "Discordance")
            else:
                plot_empty(ax, "No replicate pairs")
        else:
            plot_empty(ax, "No replicate pairs")
        label_panel(ax, "H")

        # I HWE
        ax = axes[8]
        pvals = np.asarray((hwe or {}).get("p_values") or [], dtype=float)
        pvals = pvals[np.isfinite(pvals)]
        if pvals.size:
            ax.hist(pvals, bins=20, color="#17becf", edgecolor="white")
            ax.axvline(0.05, color="black", linestyle="--", linewidth=1)
            valid = max(1, int((hwe or {}).get("valid_loci") or pvals.size))
            rejected = (hwe or {}).get("sig_deviations")
            deficit = (hwe or {}).get("sig_deficit")
            excess = (hwe or {}).get("sig_excess")
            txt = []
            if rejected is not None:
                txt.append(f"P<0.05: {100.0 * rejected / valid:.1f}%")
            if deficit is not None:
                txt.append(f"deficit: {100.0 * deficit / valid:.1f}%")
            if excess is not None:
                txt.append(f"excess: {100.0 * excess / valid:.1f}%")
            if txt:
                ax.text(0.98, 0.95, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=8)
            ax.set_xlabel("HWE P-value", fontsize=9)
            style_ax(ax, "Loci")
        else:
            plot_empty(ax)
        label_panel(ax, "I")

        # K WHOA (placed on page 1; J is heatmap on page 2)
        ax = axes[9]
        plot_whoa_on_ax(ax, whoa)
        label_panel(ax, "K")

        # L PCA
        ax = axes[10]
        if pca and pca.get("scores") is not None and np.asarray(pca["scores"]).shape[1] >= 2:
            scores = np.asarray(pca["scores"])
            ax.scatter(scores[:, 0], scores[:, 1], s=18, alpha=0.75, color="#1f77b4")
            ev = np.asarray(pca.get("explained_variance_ratio", [np.nan, np.nan]))
            xlab = "PC1"
            ylab = "PC2"
            if ev.size >= 1 and np.isfinite(ev[0]):
                xlab = f"PC1 ({100 * ev[0]:.1f}%)"
            if ev.size >= 2 and np.isfinite(ev[1]):
                ylab = f"PC2 ({100 * ev[1]:.1f}%)"
            ax.set_xlabel(xlab, fontsize=9)
            style_ax(ax, ylab)
        else:
            plot_empty(ax, "PCA unavailable")
        label_panel(ax, "L")

        # M HDplot
        ax = axes[11]
        h = np.asarray(data.get("hdplot_h") or [], dtype=float)
        d = np.asarray(data.get("hdplot_d") or [], dtype=float)
        if h.size and d.size and h.size == d.size:
            flagged = np.array([hdplot_is_paralog(hi, di) for hi, di in zip(h, d)])
            n_tested = int(data.get("hdplot_n_tested") or h.size)
            n_flagged = int(data.get("hdplot_n_flagged") if data.get("hdplot_n_flagged") is not None else flagged.sum())
            if h.size > 20000:
                rng = np.random.default_rng(0)
                keep = rng.choice(h.size, 20000, replace=False)
                h, d, flagged = h[keep], d[keep], flagged[keep]
            ax.scatter(d[~flagged], h[~flagged], s=6, alpha=0.2, color="0.55", rasterized=True, label="ok")
            if flagged.any():
                ax.scatter(d[flagged], h[flagged], s=8, alpha=0.5, color="#d62728", rasterized=True, label="flagged")
            ax.axhline(HDPLOT_H_CUTOFF, color="black", linestyle="--", linewidth=0.8)
            ax.axvline(HDPLOT_D_CUTOFF, color="black", linestyle="--", linewidth=0.8)
            ax.axvline(-HDPLOT_D_CUTOFF, color="black", linestyle="--", linewidth=0.8)
            ax.set_xlabel("D (read-ratio deviation)", fontsize=9)
            style_ax(ax, "H (het. proportion)")
            ax.legend(fontsize=7, frameon=False, loc="upper right")
            ax.text(
                0.02, 0.98,
                f"flagged {n_flagged}/{n_tested}"
                + (f" ({100.0 * n_flagged / n_tested:.1f}%)" if n_tested else ""),
                transform=ax.transAxes, ha="left", va="top", fontsize=8,
            )
        else:
            plot_empty(ax, "No HDplot depths")
        label_panel(ax, "M")
        fig.suptitle("Supplementary Figure 2. Final assembly QC", fontsize=13, y=0.995)
        fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.97])
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2: shared-loci heatmap (panel J)
        from shared_loci_matrix import plot_clustered_heatmap
        if shared_samples is not None and shared_matrix is not None:
            heat_fig = plot_clustered_heatmap(None, shared_samples, shared_matrix)
            heat_fig.suptitle("J. Pairwise shared loci", fontsize=13, y=0.995)
            pdf.savefig(heat_fig)
            plt.close(heat_fig)
        else:
            fig2, ax2 = plt.subplots(figsize=(8, 6))
            plot_empty(ax2, "Shared-loci matrix unavailable")
            label_panel(ax2, "J")
            pdf.savefig(fig2)
            plt.close(fig2)

    print(f"Supplementary assembly QC figure saved to {output_pdf}")


def resolve_sidecar(explicit, candidates, input_csv):
    if explicit and os.path.exists(explicit):
        print(f"Loaded sidecar from {explicit}")
        return pd.read_csv(explicit)
    for path in candidates:
        if os.path.exists(path):
            print(f"Loaded sidecar from {path}")
            return pd.read_csv(path)
    return None


def build_vcf_candidates(input_csv):
    input_dir = os.path.dirname(input_csv) if os.path.dirname(input_csv) else "."
    input_basename = os.path.basename(input_csv)
    input_stem = os.path.splitext(input_basename)[0]
    candidates = []
    add_candidate(candidates, os.path.join(input_dir, input_basename.replace(".csv", "_vcf.csv")), input_csv)
    rep1 = input_basename.replace("_sample_stats.csv", "_sample_stats_vcf.csv")
    if rep1 != input_basename:
        add_candidate(candidates, os.path.join(input_dir, rep1), input_csv)
    rep2 = input_basename.replace("_sample_stats_", "_sample_stats_vcf_")
    if rep2 != input_basename:
        add_candidate(candidates, os.path.join(input_dir, rep2), input_csv)
    if input_stem.startswith("optimization_sample_stats_"):
        profile_suffix = input_stem.replace("optimization_sample_stats_", "", 1)
        add_candidate(candidates, os.path.join(input_dir, f"optimization_sample_stats_vcf_{profile_suffix}.csv"), input_csv)
    add_candidate(candidates, os.path.join(input_dir, "optimization_sample_stats_vcf.csv"), input_csv)
    return candidates


def build_hwe_candidates(input_csv):
    input_dir = os.path.dirname(input_csv) if os.path.dirname(input_csv) else "."
    input_basename = os.path.basename(input_csv)
    input_stem = os.path.splitext(input_basename)[0]
    candidates = []
    add_candidate(candidates, os.path.join(input_dir, input_basename.replace(".csv", "_hwe.csv")), input_csv)
    hrep1 = input_basename.replace("_sample_stats.csv", "_hwe_stats.csv")
    if hrep1 != input_basename:
        add_candidate(candidates, os.path.join(input_dir, hrep1), input_csv)
    hrep2 = input_basename.replace("_sample_stats_", "_hwe_stats_")
    if hrep2 != input_basename:
        add_candidate(candidates, os.path.join(input_dir, hrep2), input_csv)
    add_candidate(candidates, os.path.join(input_dir, "optimization_hwe_stats.csv"), input_csv)
    if input_stem.startswith("optimization_sample_stats_"):
        profile_suffix = input_stem.replace("optimization_sample_stats_", "", 1)
        add_candidate(candidates, os.path.join(input_dir, f"optimization_hwe_stats_{profile_suffix}.csv"), input_csv)
    return candidates


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Plot optimization/QC statistics")
    parser.add_argument("--layout", choices=["diagnostic", "supplementary"], default="diagnostic")
    parser.add_argument("input_csv", nargs="?", default="optimization_sample_stats.csv")
    parser.add_argument("output_pdf", nargs="?", default="optimization_sample_stats.pdf")
    parser.add_argument("legacy_vcf_csv", nargs="?", default=None)
    parser.add_argument("legacy_hwe_csv", nargs="?", default=None)
    parser.add_argument("--vcf", dest="vcf_csv", default=None)
    parser.add_argument("--hwe", dest="hwe_csv", default=None)
    parser.add_argument("--replicate", dest="replicate_csv", default=None)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    input_csv = args.input_csv
    output_pdf = args.output_pdf
    vcf_csv = args.vcf_csv or args.legacy_vcf_csv
    hwe_csv = args.hwe_csv or args.legacy_hwe_csv

    if not os.path.exists(input_csv):
        print(f"No data to plot: missing {input_csv}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(input_csv)
    vcf_df = resolve_sidecar(vcf_csv, build_vcf_candidates(input_csv), input_csv)
    hwe_df = resolve_sidecar(hwe_csv, build_hwe_candidates(input_csv), input_csv)
    replicate_df = load_optional_csv(args.replicate_csv)
    if replicate_df is None and input_csv:
        input_dir = os.path.dirname(input_csv) if os.path.dirname(input_csv) else "."
        input_stem = os.path.splitext(os.path.basename(input_csv))[0]
        if input_stem.startswith("optimization_sample_stats_"):
            profile_suffix = input_stem.replace("optimization_sample_stats_", "", 1)
            replicate_df = load_optional_csv(
                os.path.join(input_dir, f"optimization_replicate_stats_{profile_suffix}.csv")
            )

    if args.layout == "supplementary":
        plot_threshold_selection_supp(df, vcf_df, hwe_df, replicate_df, output_pdf)
    else:
        plot_diagnostic(df, vcf_df, hwe_df, output_pdf)


if __name__ == "__main__":
    main()
