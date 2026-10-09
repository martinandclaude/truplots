#!/usr/bin/env python3
"""
plot_vaf.py - Plot the variant allele frequency (VAF) distribution from a VCF/BCF
as a histogram, green for HET and blue for HOM.

You expect a clean germline sample to be bimodal: a HET peak near 0.5 and a HOM
peak near 1.0. Dashed guides are drawn at both.

VAF source (per record, in this order unless overridden):
  1. FORMAT/AD  -> sum(alt depths) / sum(all depths)   [recommended; caller-agnostic]
  2. FORMAT/AF or FORMAT/VAF                            [fallback if AD is absent]

Zygosity comes from FORMAT/GT:
  HET = alleles differ (0/1, 1/2, 0|1, ...)
  HOM = all alleles identical and non-ref (1/1, 2/2, or hemizygous 1)
  hom-ref (0/0) and no-call (./.) are skipped.

Examples:
  plot_vaf.py sample.vcf.gz
  plot_vaf.py sample.vcf.gz -s TUMOR -o vaf.pdf --min-dp 20 --pass-only --bins 80
  plot_vaf.py sample.vcf.gz --vaf-source af --title "NA12878 allele balance"
"""

import argparse
import os
import sys

import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HET_COLOR = "#2ca02c"   # green
HOM_COLOR = "#1f77b4"   # blue


def classify_gt(gt):
    """Return 'HET', 'HOM' or None (hom-ref / no-call) from a pysam GT tuple."""
    alleles = [a for a in gt if a is not None]
    if not alleles:
        return None                       # ./. no-call
    if all(a == 0 for a in alleles):
        return None                       # 0/0 hom-ref -> not a variant
    if len(set(alleles)) == 1:
        return "HOM"                      # 1/1, 2/2, hemizygous 1
    return "HET"                          # 0/1, 1/2, ...


def get_vaf(call, prefer):
    """Return (vaf, depth) for one sample call, or (None, None) if unavailable.

    prefer='ad' -> AD then AF/VAF;  prefer='af' -> AF/VAF then AD.
    """
    def from_ad():
        ad = call.get("AD")
        if ad is None:
            return None
        ad = [a for a in ad if a is not None]
        total = sum(ad)
        if len(ad) < 2 or total <= 0:
            return None
        alt = sum(ad[1:])                 # combined alt depth (covers multiallelic)
        return alt / total, total

    def from_af():
        for key in ("AF", "VAF"):
            af = call.get(key)
            if af is None:
                continue
            val = af[0] if isinstance(af, (tuple, list)) else af
            if val is None:
                continue
            dp = call.get("DP")
            return float(val), (int(dp) if dp is not None else None)
        return None

    order = (from_af, from_ad) if prefer == "af" else (from_ad, from_af)
    for fn in order:
        res = fn()
        if res is not None:
            return res
    return None, None


def collect(vcf_path, sample, prefer, min_dp, pass_only):
    """Iterate the VCF and return (records, sample, stats).

    records: list of (vaf, zyg).
    """
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    if not samples:
        sys.exit("ERROR: VCF has no genotype columns - cannot compute per-variant VAF.")
    if sample is None:
        sample = samples[0]
    elif sample not in samples:
        sys.exit(f"ERROR: sample '{sample}' not found. Available: {', '.join(samples)}")

    records = []
    stats = dict(filtered=0, nocall=0, novaf=0, lowdp=0)
    for rec in vcf:
        if pass_only:
            f = list(rec.filter.keys())
            if f and f != ["PASS"]:
                stats["filtered"] += 1
                continue
        call = rec.samples[sample]
        gt = call.get("GT")
        if gt is None:
            stats["nocall"] += 1
            continue
        zyg = classify_gt(gt)
        if zyg is None:
            stats["nocall"] += 1
            continue
        vaf, depth = get_vaf(call, prefer)
        if vaf is None:
            stats["novaf"] += 1
            continue
        if min_dp > 0 and (depth is None or depth < min_dp):
            stats["lowdp"] += 1
            continue
        records.append((float(vaf), zyg))
    vcf.close()
    return records, sample, stats


def plot(records, sample, out, title, dpi, nbins):
    vaf = np.array([r[0] for r in records], dtype=float)
    zyg = np.array([r[1] for r in records])
    edges = np.linspace(0, 1, nbins + 1)

    fig, ax = plt.subplots(figsize=(8, 5))

    counts = {}
    for label, color in (("HET", HET_COLOR), ("HOM", HOM_COLOR)):
        m = zyg == label
        counts[label] = int(m.sum())
        ax.hist(vaf[m], bins=edges, color=color, alpha=0.6,
                histtype="stepfilled", label=f"{label} (n={int(m.sum()):,})")

    for xline in (0.5, 1.0):               # expected HET / HOM centres
        ax.axvline(xline, color="0.65", lw=0.8, ls="--", zorder=0)

    ax.set_xlim(0, 1)
    ax.set_xlabel("Variant allele frequency")
    ax.set_ylabel("Number of variants")
    ax.set_title(title or f"VAF distribution - {sample}", fontweight="bold")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return counts


def default_output(vcf_path):
    b = os.path.basename(vcf_path)
    for ext in (".vcf.gz", ".vcf.bgz", ".bcf", ".vcf"):
        if b.endswith(ext):
            b = b[: -len(ext)]
            break
    return b + ".vaf_hist.png"


def main():
    ap = argparse.ArgumentParser(
        description="Plot the VAF distribution from a VCF as a histogram, green=HET / blue=HOM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("vcf", help="Input VCF / VCF.gz / BCF")
    ap.add_argument("-o", "--output", help="Output image (.png/.pdf/.svg). "
                    "Default: <vcf>.vaf_hist.png")
    ap.add_argument("-s", "--sample", help="Sample name (default: first sample in VCF)")
    ap.add_argument("--vaf-source", choices=["ad", "af"], default="ad",
                    help="'ad' = compute from FORMAT/AD (fallback to AF/VAF); "
                         "'af' = use FORMAT/AF or VAF (fallback to AD)")
    ap.add_argument("--bins", type=int, default=50, help="Number of histogram bins over 0-1")
    ap.add_argument("--min-dp", type=int, default=0,
                    help="Skip variants below this depth (0 = no filter; 20 is common)")
    ap.add_argument("--pass-only", action="store_true",
                    help="Keep only PASS (or unfiltered) variants")
    ap.add_argument("--title", help="Custom plot title")
    ap.add_argument("--dpi", type=int, default=150, help="Output resolution")
    args = ap.parse_args()

    if not os.path.exists(args.vcf):
        sys.exit(f"ERROR: file not found: {args.vcf}")
    if args.bins < 1:
        sys.exit("ERROR: --bins must be >= 1")
    out = args.output or default_output(args.vcf)

    records, sample, stats = collect(
        args.vcf, args.sample, args.vaf_source, args.min_dp, args.pass_only)

    if not records:
        sys.exit("ERROR: no plottable variants found "
                 f"(skipped: {stats}). Check sample, --vaf-source, or filters.")

    counts = plot(records, sample, out, args.title, args.dpi, args.bins)

    print(f"Sample plotted : {sample}", file=sys.stderr)
    print(f"Variants plotted: {len(records):,} "
          f"(HET={counts['HET']:,}, HOM={counts['HOM']:,})", file=sys.stderr)
    skipped = (f"filtered={stats['filtered']}, nocall/hom-ref={stats['nocall']}, "
               f"no-VAF={stats['novaf']}, low-DP={stats['lowdp']}")
    print(f"Skipped        : {skipped}", file=sys.stderr)
    print(f"Wrote          : {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
