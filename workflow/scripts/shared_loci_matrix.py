#!/usr/bin/env python3
"""Build pairwise shared-loci matrix and clustered heatmap from an ipyrad .loci file."""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import squareform


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute sample x sample shared-locus counts from an ipyrad .loci file, "
            "write a TSV matrix, and plot a clustered heatmap with similarity trees."
        )
    )
    parser.add_argument("--loci", required=True, help="Input .loci file")
    parser.add_argument("--matrix", required=True, help="Output TSV matrix path")
    parser.add_argument("--heatmap", required=True, help="Output heatmap PDF path")
    return parser.parse_args()


def parse_loci_presence(loci_path):
    """
    Parse .loci into sample order and a binary presence matrix (samples x loci).
    Each block between // lines is one locus; sample names are the first token.
    """
    samples = []
    sample_index = {}
    loci_presence = []
    current = []

    with open(loci_path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("//"):
                if current:
                    loci_presence.append(current)
                    current = []
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name = parts[0]
            if name not in sample_index:
                sample_index[name] = len(samples)
                samples.append(name)
            current.append(sample_index[name])

    if current:
        loci_presence.append(current)

    n_samples = len(samples)
    n_loci = len(loci_presence)
    presence = np.zeros((n_samples, n_loci), dtype=np.uint8)
    for locus_idx, idxs in enumerate(loci_presence):
        presence[idxs, locus_idx] = 1
    return samples, presence


def shared_loci_matrix(presence):
    """Return integer shared-locus counts (samples x samples)."""
    if presence.size == 0:
        return np.zeros((presence.shape[0], presence.shape[0]), dtype=np.int32)
    return (presence.astype(np.int32) @ presence.T)


def mean_pairwise_shared(shared):
    """Mean off-diagonal pairwise shared-locus count."""
    n = shared.shape[0]
    if n < 2:
        return float("nan")
    triu = np.triu_indices(n, k=1)
    return float(np.mean(shared[triu]))


def jaccard_distance(shared):
    """
    Pairwise Jaccard distance from shared counts.
    Diagonal entries are sample totals; D_ii = 0.
    """
    n = shared.shape[0]
    dist = np.zeros((n, n), dtype=np.float64)
    totals = np.diag(shared).astype(np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            union = totals[i] + totals[j] - shared[i, j]
            if union <= 0:
                d = 0.0 if totals[i] == 0 and totals[j] == 0 else 1.0
            else:
                d = 1.0 - (float(shared[i, j]) / union)
            dist[i, j] = dist[j, i] = d
    return dist


def write_matrix_tsv(path, samples, shared):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as out:
        out.write("sample\t" + "\t".join(samples) + "\n")
        for i, name in enumerate(samples):
            row = "\t".join(str(int(v)) for v in shared[i])
            out.write(f"{name}\t{row}\n")


def plot_clustered_heatmap(path, samples, shared):
    """Heatmap of shared-locus counts with average-linkage trees (Jaccard distance).

    If path is None, the figure is returned without saving or closing.
    """
    if path:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    n = len(samples)
    if n == 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No samples found in .loci file", ha="center", va="center")
        ax.axis("off")
        if path:
            fig.savefig(path, bbox_inches="tight")
            plt.close(fig)
        return fig

    if n == 1:
        order = [0]
        Z = None
    else:
        dist = jaccard_distance(shared)
        # Guard against tiny negative values from floating point
        np.fill_diagonal(dist, 0.0)
        dist = np.clip(dist, 0.0, None)
        # If all distances are zero, skip meaningful clustering
        if np.allclose(dist, 0.0):
            order = list(range(n))
            Z = None
        else:
            condensed = squareform(dist, checks=False)
            Z = linkage(condensed, method="average")
            order = leaves_list(Z).tolist()

    ordered = shared[np.ix_(order, order)]
    ordered_names = [samples[i] for i in order]

    # Scale figure with sample count; keep labels readable for moderate n
    label_size = 8 if n <= 80 else max(3, int(200 / n))
    fig_w = max(8.0, min(28.0, 0.22 * n + 4.0))
    fig_h = max(8.0, min(28.0, 0.22 * n + 4.0))
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        2,
        2,
        width_ratios=[1.0, 4.0],
        height_ratios=[1.0, 4.0],
        wspace=0.02,
        hspace=0.02,
        figure=fig,
    )

    ax_col = fig.add_subplot(gs[0, 1])
    ax_row = fig.add_subplot(gs[1, 0])
    ax_mat = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_axes([0.92, 0.15, 0.02, 0.3])

    if Z is not None:
        dendrogram(Z, ax=ax_col, orientation="top", no_labels=True, color_threshold=0, above_threshold_color="0.3")
        dendrogram(Z, ax=ax_row, orientation="left", no_labels=True, color_threshold=0, above_threshold_color="0.3")
    ax_col.set_xticks([])
    ax_col.set_yticks([])
    ax_col.axis("off")
    ax_row.set_xticks([])
    ax_row.set_yticks([])
    ax_row.axis("off")

    im = ax_mat.imshow(ordered, aspect="auto", interpolation="nearest", cmap="viridis")
    ax_mat.set_xticks(range(n))
    ax_mat.set_yticks(range(n))
    ax_mat.set_xticklabels(ordered_names, rotation=90, fontsize=label_size)
    ax_mat.set_yticklabels(ordered_names, fontsize=label_size)
    ax_mat.tick_params(length=0)
    ax_mat.set_title("Pairwise shared loci (clustered by Jaccard similarity)", fontsize=11, pad=8)

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Shared loci", fontsize=9)

    if path:
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
    return fig


def main():
    args = parse_args()
    if not os.path.exists(args.loci):
        print(f"ERROR: loci file not found: {args.loci}", file=sys.stderr)
        sys.exit(1)

    samples, presence = parse_loci_presence(args.loci)
    shared = shared_loci_matrix(presence)
    write_matrix_tsv(args.matrix, samples, shared)
    plot_clustered_heatmap(args.heatmap, samples, shared)

    totals = np.diag(shared) if shared.size else np.array([])
    print(
        f"Wrote shared-loci matrix for {len(samples)} samples "
        f"({int(presence.shape[1])} loci); "
        f"per-sample loci min/median/max = "
        f"{int(totals.min()) if totals.size else 0}/"
        f"{int(np.median(totals)) if totals.size else 0}/"
        f"{int(totals.max()) if totals.size else 0}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
