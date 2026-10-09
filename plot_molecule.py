#!/usr/bin/env python3
"""
plot_molecule.py - find and plot the largest TruPath 'molecule' (proximity
constellation) in a genomic region, as a 2-panel figO-style figure:

    LEFT  : the molecule's nanowells on the flow cell (coloured by position
            within the molecule)
    RIGHT : the same reads laid out on the genome = the reconstructed molecule

A 'molecule' is recovered by spatially+genomically clustering reads within each
flow-cell tile (DBSCAN on x, y, and scaled genomic position): a group of
neighbouring nanowells that also map to a contiguous genomic segment is one
original long DNA molecule. Meaningful for TruPath / proximity-mapped data;
on standard WGS there are no real constellations.

Usage:
    python3 plot_molecule.py --cram sample.cram --reference genome.fa \\
        --region chr2:178525989-178830800 [--out fig.png] [--by span|reads]

Requires: samtools on PATH; numpy, matplotlib, scikit-learn.
All processing is local.
"""
import argparse, subprocess, sys, shutil, collections
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Plot the largest proximity molecule in a region (figO-style).")
    p.add_argument("--cram", required=True, help="CRAM/BAM file (TruPath / proximity data).")
    p.add_argument("--reference", "-T", required=True, help="Reference FASTA (for CRAM decoding).")
    p.add_argument("--region", required=True, help="Region, e.g. chr2:178525989-178830800.")
    p.add_argument("--out", help="Output PNG (default: molecule_<region>.png).")
    p.add_argument("--by", choices=["span", "reads", "density"], default="span",
                   help="Rank molecules by: span = longest (default); reads = most read-pairs; "
                        "density = sparsest (lowest read-pairs per 100 kb = biggest read-pair gaps).")
    p.add_argument("--top", type=int, default=1,
                   help="Plot the top N molecules as a stacked gallery (default 1 = single figO plot).")
    p.add_argument("--min-reads", type=int, default=5, help="Min read-pairs per molecule (default 5).")
    p.add_argument("--eps", type=float, default=350, help="DBSCAN spatial radius in px (default 350).")
    p.add_argument("--genomic-eps", type=float, default=60000,
                   help="Genomic distance (bp) mapped to one --eps in clustering (default 60000).")
    p.add_argument("--mapq", type=int, default=1, help="Minimum mapping quality (default 1).")
    p.add_argument("--max-span-kb", type=float, default=None,
                   help="Optional cap on molecule span (kb) to reject merged/outlier clusters.")
    return p.parse_args()


def load_reads(cram, reference, region, mapq):
    """Return mate-collapsed reads as arrays: lane, tile(str), x, y, pos."""
    if shutil.which("samtools") is None:
        sys.exit("error: samtools not found on PATH.")
    cmd = ["samtools", "view", "-F", "3328", "-q", str(mapq), "-T", reference, cram, region]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.exit(f"error: samtools failed:\n{proc.stderr.strip()}")
    seen = set(); lane=[]; tile=[]; x=[]; y=[]; pos=[]
    for ln in proc.stdout.splitlines():
        f = ln.split("\t")
        name = f[0]
        if name in seen:                      # collapse read mates (share name + x/y)
            continue
        seen.add(name)
        n = name.split(":")
        if len(n) < 7:
            continue
        try:
            lane.append(int(n[3])); tile.append(n[4])
            x.append(int(n[5])); y.append(int(n[6])); pos.append(int(f[3]))
        except ValueError:
            continue
    if not pos:
        sys.exit(f"error: no usable reads in {region} (check region, reference, and that read names "
                 "carry Illumina LANE:TILE:X:Y coordinates).")
    return (np.array(lane), np.array(tile), np.array(x, float), np.array(y, float), np.array(pos, float))


def find_molecules(reads, eps, genomic_eps, min_reads, max_span_kb):
    from sklearn.cluster import DBSCAN
    lane, tile, x, y, pos = reads
    scale = eps / genomic_eps
    mols = []
    keys = collections.defaultdict(list)
    for i in range(len(pos)):
        keys[(lane[i], tile[i])].append(i)
    for idx in keys.values():
        if len(idx) < min_reads:
            continue
        idx = np.array(idx)
        feat = np.column_stack([x[idx], y[idx], pos[idx] * scale])
        labels = DBSCAN(eps=eps, min_samples=min_reads).fit_predict(feat)
        for cl in set(labels):
            if cl == -1:
                continue
            sel = idx[labels == cl]
            span = (pos[sel].max() - pos[sel].min()) / 1e3
            if max_span_kb and span > max_span_kb:
                continue
            gaps = np.diff(np.sort(pos[sel])) / 1e3            # genomic gaps between read-pairs (kb)
            mols.append(dict(n=len(sel), span=span,
                             x=x[sel], y=y[sel], pos=pos[sel],
                             gap_med=float(np.median(gaps)), gap_max=float(gaps.max()),
                             reads_per_100kb=len(sel) / max(span, 1e-9) * 100,
                             lane=int(lane[sel][0]), tile=tile[sel][0]))
    return mols


def plot_molecule(m, region, cram, out, by="span"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    gp = (m["pos"] - m["pos"].min()) / 1e3            # kb within molecule
    xext = m["x"].max() - m["x"].min(); yext = m["y"].max() - m["y"].min()
    cx, cy = m["x"].mean(), m["y"].mean(); half = max(xext, yext) * 0.8 + 220

    fig, ax = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.5]})
    # --- left: constellation on the flow cell ---
    a = ax[0]
    sc = a.scatter(m["x"], m["y"], c=gp, cmap="turbo", s=220, ec="k", lw=.7, zorder=3)
    a.add_patch(Ellipse((cx, cy), xext + 280, yext + 280, fill=False, ec="#b00", lw=1.4, ls="--", zorder=2))
    for i in range(m["n"]):
        a.annotate(f"{gp[i]:.0f}", (m["x"][i], m["y"][i]), fontsize=7, ha="center", va="center", zorder=4)
    a.set_xlim(cx - half, cx + half); a.set_ylim(cy - half, cy + half)
    sb = 500
    a.plot([cx - half + 150, cx - half + 150 + sb], [cy - half + 120] * 2, "k-", lw=3)
    a.text(cx - half + 150, cy - half + 200, f"{sb} px", fontsize=8)
    a.set_xlabel("flow-cell X (px)"); a.set_ylabel("flow-cell Y (px)")
    a.set_title(f"ONE molecule on the flow cell\nlane {m['lane']}, tile {m['tile']}: "
                f"{m['n']} nanowells in {xext:.0f}x{yext:.0f} px", fontsize=10)
    cb = fig.colorbar(sc, ax=a); cb.set_label("genomic position within molecule (kb)")
    # --- right: reconstructed genomic segment ---
    b = ax[1]
    b.hlines(0, gp.min(), gp.max(), color="#b00", lw=2, zorder=1)
    b.scatter(gp, np.zeros(m["n"]), c=gp, cmap="turbo", s=220, ec="k", lw=.7, zorder=3)
    b.vlines(gp, -0.04, 0.04, color="#888", lw=.6, zorder=2)
    b.annotate("", (gp.min(), 0.2), (gp.max(), 0.2), arrowprops=dict(arrowstyle="<->", color="#333"))
    b.text((gp.min() + gp.max()) / 2, 0.25,
           f"reconstructed molecule ≈ {m['span']:.0f} kb, {m['n']} read-pairs", ha="center", fontsize=10)
    b.text((gp.min() + gp.max()) / 2, -0.32,
           f"read-pair gap: median {m['gap_med']:.1f} kb, max {m['gap_max']:.0f} kb  "
           f"({m['reads_per_100kb']:.1f} read-pairs/100 kb)", ha="center", fontsize=9, color="#555")
    b.set_ylim(-0.6, 0.5); b.set_yticks([])
    chrom = region.split(":")[0]
    b.set_xlabel(f"genomic position within molecule (kb)   [{chrom}:{m['pos'].min()/1e6:.3f} Mb start]")
    b.set_title("...the same reads tile a contiguous genomic segment.", fontsize=10)
    noun = {"span": "Largest molecule (longest span)",
            "reads": "Best-supported molecule (most read-pairs)",
            "density": "Sparsest molecule (fewest read-pairs per 100 kb)"}[by]
    fig.suptitle(f"{noun} in {region}   ({cram.split('/')[-1]})", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


def plot_gallery(mols, region, cram, out, by="span"):
    """figR-style stacked gallery: one row per molecule (constellation + genome)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    n = len(mols); chrom = region.split(":")[0]
    fig, ax = plt.subplots(n, 2, figsize=(13, 2.5 * n),
                           gridspec_kw={"width_ratios": [1, 1.6]}, squeeze=False)
    for k, m in enumerate(mols):
        gp = (m["pos"] - m["pos"].min()) / 1e3
        xext = m["x"].max() - m["x"].min(); yext = m["y"].max() - m["y"].min()
        cx, cy = m["x"].mean(), m["y"].mean(); half = max(xext, yext) * 0.7 + 200
        a = ax[k, 0]
        a.scatter(m["x"], m["y"], c=gp, cmap="turbo", s=110, ec="k", lw=.5)
        a.add_patch(Ellipse((cx, cy), xext + 250, yext + 250, fill=False, ec="#b00", lw=1, ls="--"))
        a.set_xlim(cx - half, cx + half); a.set_ylim(cy - half, cy + half)
        sb = 500
        a.plot([cx - half + 70, cx - half + 70 + sb], [cy - half + 60] * 2, "k-", lw=2.5)
        a.text(cx - half + 70, cy - half + 110, f"{sb}px", fontsize=6.5)
        a.set_xticks([]); a.set_yticks([]); a.set_ylabel(f"#{k+1}", fontsize=9)
        a.set_title(f"flow cell: {m['n']} nanowells, lane {m['lane']} tile {m['tile']}", fontsize=8)
        b = ax[k, 1]
        b.hlines(0, gp.min(), gp.max(), color="#b00", lw=2, zorder=1)
        b.scatter(gp, np.zeros(m["n"]), c=gp, cmap="turbo", s=110, ec="k", lw=.5, zorder=3)
        b.vlines(gp, -0.05, 0.05, color="#999", lw=.6, zorder=2)
        b.annotate("", (gp.min(), 0.28), (gp.max(), 0.28), arrowprops=dict(arrowstyle="<->", color="#333"))
        b.text((gp.min() + gp.max()) / 2, 0.36,
               f"{m['span']:.0f} kb · {m['n']} read-pairs · gap med {m['gap_med']:.0f} kb",
               ha="center", fontsize=8.5)
        b.set_ylim(-0.7, 0.7); b.set_yticks([])
        b.set_title(f"reconstructed molecule  [{chrom}:{m['pos'].min()/1e6:.3f} Mb]", fontsize=8)
        if k == n - 1:
            b.set_xlabel("genomic position within molecule (kb)")
    ranked = {"span": "by span (longest)", "reads": "by read support (most read-pairs)",
              "density": "by sparsity (lowest read-pairs/100 kb)"}[by]
    fig.suptitle(f"Top {n} molecules {ranked} in {region}  ({cram.split('/')[-1]})", fontsize=12, y=1.005)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")


def main():
    args = parse_args()
    n_top = max(1, args.top)
    tag = (f"top{n_top}_{args.by}_" if n_top > 1 else f"molecule_{args.by}_")
    out = args.out or tag + args.region.replace(":", "_").replace("-", "_") + ".png"
    reads = load_reads(args.cram, args.reference, args.region, args.mapq)
    print(f"loaded {len(reads[4])} mate-collapsed reads in {args.region}")
    mols = find_molecules(reads, args.eps, args.genomic_eps, args.min_reads, args.max_span_kb)
    if not mols:
        sys.exit(f"error: no molecule with >= {args.min_reads} read-pairs found. Try a larger region, "
                 "--min-reads 3, or confirm this is proximity (TruPath) data.")
    keyfn = {"span": lambda m: m["span"], "reads": lambda m: m["n"],
             "density": lambda m: m["reads_per_100kb"]}[args.by]
    mols.sort(key=keyfn, reverse=(args.by != "density"))   # density: ascending => sparsest first
    chosen = mols[:n_top]
    print(f"found {len(mols)} molecules; plotting top {len(chosen)} by {args.by}:")
    for i, m in enumerate(chosen, 1):
        print(f"  #{i}: span {m['span']:6.1f} kb | {m['n']:2d} read-pairs | "
              f"gap med {m['gap_med']:5.1f} kb (max {m['gap_max']:5.1f}) | lane {m['lane']} tile {m['tile']}")
    if len(chosen) == 1:
        plot_molecule(chosen[0], args.region, args.cram, out, by=args.by)
    else:
        plot_gallery(chosen, args.region, args.cram, out, by=args.by)


if __name__ == "__main__":
    main()
