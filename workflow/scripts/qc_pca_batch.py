"""PCA QC plots of an ipyrad VCF, grouped by barcode and by library.

Each sample is a point. Each group is drawn as a star at the centroid, with
lines from its samples and the group name beside the star.
Samples found in more than one library, or with more than one barcode, are
labeled "multiple" because ipyrad merges those reads into one sample.
"""
import csv
import os
import random
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

try:
    from sklearn.decomposition import PCA
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

MULTI = "multiple"
UNASSIGNED = "unassigned"


def as_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return list(value)


def dosage(gt):
    """Alt-allele count, or None if the call is missing or not 0/1."""
    if not gt or "." in gt:
        return None
    alt = 0
    for allele in gt.replace("|", "/").split("/"):
        if allele == "0":
            continue
        if allele == "1":
            alt += 1
            continue
        return None
    return float(alt)


def _info_ns(info):
    for field in info.split(";"):
        if field.startswith("NS="):
            try:
                return int(field[3:])
            except ValueError:
                return None
    return None


def parse_dosage(sample_fields, gt_index):
    """Return (dosage row, alt copies, allele copies, n missing), or None."""
    n = len(sample_fields)
    row = np.empty(n, dtype=np.float32)
    alt_copies = 0
    allele_copies = 0
    missing = 0
    for i, field in enumerate(sample_fields):
        if gt_index == 0:
            gt = field.split(":", 1)[0]
        else:
            bits = field.split(":")
            gt = bits[gt_index] if len(bits) > gt_index else "."
        val = dosage(gt)
        if val is None:
            missing += 1
            row[i] = np.nan
            continue
        row[i] = val
        alt_copies += int(val)
        allele_copies += gt.replace("|", "/").count("/") + 1
    return row, alt_copies, allele_copies, missing


def sample_sites(vcf_path, max_snps, max_missing, min_maf, seed):
    """One SNP per locus, then a random subset of loci that pass filters."""
    rng = random.Random(seed)
    samples = None
    min_called = None
    reservoir = []
    used_loci = set()
    n_seen = 0
    n_pass = 0

    opener = open
    if vcf_path.endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                min_called = int(np.ceil((1.0 - max_missing) * max(len(samples), 1)))
                continue
            if samples is None:
                continue

            n_seen += 1
            if n_seen % 200000 == 0:
                print(
                    f"  scanned {n_seen} variants, {n_pass} SNPs passed filters",
                    flush=True,
                )

            parts = line.rstrip("\n").split("\t", 9)
            if len(parts) < 10:
                continue
            ref, alt = parts[3], parts[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            locus = parts[2].split("_pos", 1)[0]
            if locus in used_loci:
                continue
            ns = _info_ns(parts[7])
            if ns is not None and min_called is not None and ns < min_called:
                continue

            # FORMAT is column 8; ipyrad writes GT first. Fall back to a scan.
            fmt = parts[8].split(":")
            try:
                gt_i = fmt.index("GT")
            except ValueError:
                continue

            fields = parts[9].split("\t")
            if len(fields) != len(samples):
                continue
            parsed = parse_dosage(fields, gt_i)
            row, alt_copies, allele_copies, missing = parsed
            if allele_copies == 0 or missing / len(samples) > max_missing:
                continue
            freq = alt_copies / allele_copies
            maf = freq if freq <= 0.5 else 1.0 - freq
            if maf < min_maf:
                continue

            used_loci.add(locus)
            n_pass += 1
            if len(reservoir) < max_snps:
                reservoir.append(row)
            else:
                j = rng.randrange(n_pass)
                if j < max_snps:
                    reservoir[j] = row

    if samples is None:
        raise ValueError(f"No #CHROM header in {vcf_path}")
    print(
        f"VCF scan: {n_seen} variants, {n_pass} SNPs passed, {len(reservoir)} used",
        flush=True,
    )
    return samples, reservoir, n_seen, n_pass


def load_sample_metadata(library_rows):
    """Map each barcode-file sample name to the libraries and barcodes it appears in."""
    meta = {}
    for row in as_list(library_rows):
        library, path = row.split("\t", 1)
        with open(path) as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                sample, barcode = parts[0], parts[1]
                rec = meta.setdefault(sample, {"libraries": [], "barcodes": []})
                if library not in rec["libraries"]:
                    rec["libraries"].append(library)
                if barcode not in rec["barcodes"]:
                    rec["barcodes"].append(barcode)
    return meta


def library_rows_from_barcode_files(paths):
    rows = []
    for path in as_list(paths):
        library = os.path.splitext(os.path.basename(path))[0]
        rows.append(f"{library}\t{path}")
    return rows


def group_label(values):
    if not values:
        return UNASSIGNED
    if len(values) == 1:
        return values[0]
    return MULTI


def run_pca(genotypes, max_sample_missing, seed):
    """
    genotypes: samples x sites, with NaN for missing calls.
    Returns scores (n_kept x 2), variance fractions, and the kept-sample mask.
    """
    call_rate = np.isfinite(genotypes).mean(axis=1)
    keep_samples = call_rate >= (1.0 - max_sample_missing)
    X = genotypes[keep_samples, :]
    if X.shape[0] < 3 or X.shape[1] < 2:
        return None, None, keep_samples

    site_obs = np.isfinite(X).sum(axis=0)
    X = X[:, site_obs >= 2]
    if X.shape[1] < 2:
        return None, None, keep_samples

    col_mean = np.nanmean(X, axis=0)
    inds = np.where(~np.isfinite(X))
    X[inds] = np.take(col_mean, inds[1])

    sd = X.std(axis=0)
    variable = sd > 1e-8
    X = X[:, variable]
    if X.shape[1] < 2:
        return None, None, keep_samples
    X = (X - X.mean(axis=0)) / X.std(axis=0)

    n_comp = min(2, X.shape[0] - 1, X.shape[1])
    if HAS_SKLEARN:
        pca = PCA(n_components=n_comp, random_state=seed)
        scores = pca.fit_transform(X)
        var = pca.explained_variance_ratio_
    else:
        _u, s, vt = np.linalg.svd(X, full_matrices=False)
        scores = X @ vt[:n_comp].T
        total = np.sum(s ** 2)
        var = (s[:n_comp] ** 2) / total if total else np.zeros(n_comp)

    if n_comp == 1:
        scores = np.column_stack([scores[:, 0], np.zeros(scores.shape[0])])
        var = np.array([var[0], 0.0])
    return scores, var, keep_samples


def group_color_map(groups):
    colors = {}
    i = 0
    for group in groups:
        if group == UNASSIGNED:
            colors[group] = (0.65, 0.65, 0.65, 1.0)
        elif group == MULTI:
            colors[group] = (0.15, 0.15, 0.15, 1.0)
        else:
            hue = (i * 0.61803398875) % 1.0
            colors[group] = plt.matplotlib.colors.hsv_to_rgb((hue, 0.62, 0.82))
            i += 1
    return colors


def _repel(xy, min_dist, max_shift, n_iter=300):
    pos = np.array(xy, dtype=float, copy=True)
    if len(pos) <= 1:
        return pos
    # Nudge exact ties so two groups that share a centroid can move apart.
    pos[:, 0] += (np.arange(len(pos)) - (len(pos) - 1) / 2.0) * (min_dist * 0.05)
    origin = pos.copy()
    n = len(pos)
    for _ in range(n_iter):
        delta = pos[:, None, :] - pos[None, :, :]
        dist = np.sqrt((delta * delta).sum(axis=2))
        np.fill_diagonal(dist, np.inf)
        overlap = dist < min_dist
        if not np.any(overlap):
            break
        safe = np.where(np.isfinite(dist) & (dist > 0), dist, 1.0)
        weight = np.where(overlap, (min_dist - safe) / safe, 0.0)
        shift = (delta * weight[:, :, None]).sum(axis=1) / n
        pos = pos + shift
        drift = pos - origin
        drift_norm = np.sqrt((drift * drift).sum(axis=1))
        too_far = drift_norm > max_shift
        if np.any(too_far):
            scale = np.ones(n)
            scale[too_far] = max_shift / drift_norm[too_far]
            pos = origin + drift * scale[:, None]
    return pos


def write_message_pdf(path, message):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_grouped(rows, group_key, var, out_pdf, title):
    plotted = [row for row in rows if row["in_pca"]]
    if len(plotted) < 3:
        write_message_pdf(out_pdf, f"{title}\nNot enough samples for PCA.")
        return

    groups = sorted(
        {row[group_key] for row in plotted},
        key=lambda name: (name == UNASSIGNED, name == MULTI, name),
    )
    colors = group_color_map(groups)
    fig, ax = plt.subplots(figsize=(10, 8))

    centroids = []
    for group in sorted(groups, key=lambda name: -sum(row[group_key] == name for row in plotted)):
        members = [row for row in plotted if row[group_key] == group]
        xy = np.array([[row["PC1"], row["PC2"]] for row in members], dtype=float)
        color = colors[group]
        centroid = xy.mean(axis=0)
        centroids.append((centroid, group, color))
        segments = np.stack([xy, np.repeat(centroid[None, :], len(xy), axis=0)], axis=1)
        ax.add_collection(LineCollection(
            segments, colors=[color], linewidths=0.4, alpha=0.28, zorder=1
        ))
        ax.scatter(xy[:, 0], xy[:, 1], s=16, color=color, alpha=0.8, linewidths=0, zorder=2)

    centers = np.array([c[0] for c in centroids], dtype=float)
    span = max(float(np.ptp(centers[:, 0])), float(np.ptp(centers[:, 1])), 1e-6)
    label_xy = _repel(centers, min_dist=0.055 * span, max_shift=0.18 * span)
    fontsize = 8 if len(groups) <= 20 else 6.5 if len(groups) <= 45 else 5.5

    for (centroid, group, color), (lx, ly) in zip(centroids, label_xy):
        ax.scatter(
            [centroid[0]], [centroid[1]], marker="*", s=220, color=color,
            edgecolors="black", linewidths=0.5, zorder=4,
        )
        if (lx - centroid[0]) ** 2 + (ly - centroid[1]) ** 2 > (0.015 * span) ** 2:
            ax.plot(
                [centroid[0], lx], [centroid[1], ly],
                color="0.45", lw=0.4, zorder=3,
            )
        ax.text(
            lx, ly, group, color=color, fontsize=fontsize,
            ha="center", va="center", zorder=5,
            path_effects=[pe.withStroke(linewidth=2.8, foreground="white")],
        )

    var1 = 100.0 * float(var[0])
    var2 = 100.0 * float(var[1])
    ax.set_xlabel(f"PC1 ({var1:.1f}%)")
    ax.set_ylabel(f"PC2 ({var2:.1f}%)")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    ax.autoscale_view()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def write_scores(path, rows, var, n_snps):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fields = [
        "sample", "in_pca", "PC1", "PC2", "var_PC1", "var_PC2", "n_snps",
        "barcode", "library", "barcodes_all", "libraries_all",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["var_PC1"] = f"{float(var[0]):.6f}" if var is not None else ""
            out["var_PC2"] = f"{float(var[1]):.6f}" if var is not None else ""
            out["n_snps"] = n_snps
            if not row["in_pca"]:
                out["PC1"] = ""
                out["PC2"] = ""
            writer.writerow(out)


def run(
    vcf,
    library_rows,
    library_pdf,
    barcode_pdf,
    scores_path,
    max_snps=5000,
    max_missing=0.5,
    min_maf=0.05,
    max_sample_missing=0.5,
    seed=1,
):
    for path in (library_pdf, barcode_pdf, scores_path):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    print(f"Reading {vcf}", flush=True)
    samples, sites, n_seen, n_pass = sample_sites(
        vcf, max_snps, max_missing, min_maf, seed
    )
    print(f"Variants seen: {n_seen}; passing filters: {n_pass}; used: {len(sites)}")

    meta = load_sample_metadata(library_rows)
    print(f"Barcode files assign {len(meta)} sample names")

    rows = []
    for sample in samples:
        rec = meta.get(sample, {"libraries": [], "barcodes": []})
        rows.append({
            "sample": sample,
            "in_pca": 0,
            "PC1": "",
            "PC2": "",
            "barcode": group_label(rec["barcodes"]),
            "library": group_label(rec["libraries"]),
            "barcodes_all": ",".join(rec["barcodes"]),
            "libraries_all": ",".join(rec["libraries"]),
        })

    var = None
    if not sites:
        message = "No biallelic SNPs passed the missingness and MAF filters."
        print(message)
        write_message_pdf(library_pdf, message)
        write_message_pdf(barcode_pdf, message)
        write_scores(scores_path, rows, None, 0)
        return

    genotypes = np.vstack(sites).T
    if genotypes.shape[0] != len(samples):
        raise ValueError("Genotype matrix does not match the VCF sample count")
    scores, var, keep = run_pca(genotypes, max_sample_missing, seed)
    n_snps = int(genotypes.shape[1])

    if scores is None:
        message = "Too few samples or variable sites for PCA."
        print(message)
        write_message_pdf(library_pdf, message)
        write_message_pdf(barcode_pdf, message)
        write_scores(scores_path, rows, None, n_snps)
        return

    kept_i = 0
    for i, row in enumerate(rows):
        if not keep[i]:
            continue
        row["in_pca"] = 1
        row["PC1"] = f"{scores[kept_i, 0]:.6f}"
        row["PC2"] = f"{scores[kept_i, 1]:.6f}"
        kept_i += 1

    plot_rows = []
    for row in rows:
        plot_row = dict(row)
        if row["in_pca"]:
            plot_row["PC1"] = float(row["PC1"])
            plot_row["PC2"] = float(row["PC2"])
        plot_rows.append(plot_row)

    n_in = sum(row["in_pca"] for row in rows)
    print(f"Samples in PCA: {n_in} / {len(rows)}")
    print(f"PC1 {100 * var[0]:.1f}%, PC2 {100 * var[1]:.1f}%")
    print(f"Library groups: {len({row['library'] for row in plot_rows if row['in_pca']})}")
    print(f"Barcode groups: {len({row['barcode'] for row in plot_rows if row['in_pca']})}")

    plot_grouped(plot_rows, "library", var, library_pdf, "PCA colored by library")
    plot_grouped(plot_rows, "barcode", var, barcode_pdf, "PCA colored by barcode")
    write_scores(scores_path, rows, var, n_snps)
    print(f"Wrote {library_pdf}")
    print(f"Wrote {barcode_pdf}")
    print(f"Wrote {scores_path}")


def _snakemake_main():
    try:
        sm = snakemake
    except NameError:
        return
    log_path = sm.log[0]
    os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
    sys.stdout = open(log_path, "w")
    sys.stderr = sys.stdout
    run(
        vcf=sm.input.vcf,
        library_rows=library_rows_from_barcode_files(sm.input.barcodes),
        library_pdf=sm.output.library_pdf,
        barcode_pdf=sm.output.barcode_pdf,
        scores_path=sm.output.scores,
        max_snps=int(sm.params.max_snps),
        max_missing=float(sm.params.max_missing),
        min_maf=float(sm.params.min_maf),
        max_sample_missing=float(sm.params.max_sample_missing),
        seed=int(sm.params.seed),
    )


_snakemake_main()
