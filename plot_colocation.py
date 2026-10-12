#!/usr/bin/env python3
"""
plot_colocation.py - plot DRAGEN TruPath colocation maps (<sample>.colocation.cooler) as heatmaps.

DRAGEN's TruPath colocation module splits the genome into fixed bins (default ~2 kb,
--colocation-bin-size) and counts, for every pair of bins, how many reads from the two bins
sat close together on the flow cell. The counts are stored as a cooler file. Most signal
sits on the diagonal (fragments of the same long template molecule); off-diagonal structure
can point to structural variants: deletions, inversions, and intra- or inter-chromosomal
translocations.

Views:
  genome-wide       (no --region)                     chr1-22, X, Y; auto-coarsened
  one region        --region chr5:60,000,000-80,000,000   square heatmap, or --triangle
  several regions   --region chr9 --region chr22       regions concatenated on both axes;
                                                       the off-diagonal blocks show signal
                                                       between them (e.g. a translocation)

--sv-vcf marks DRAGEN SV calls (<sample>.sv.vcf.gz) at their breakpoint pairs: DEL, DUP and
INV at (POS, END), BND at (POS, mate position). Calls go in the upper triangle only, so the
mirrored lower triangle shows the same signal unobstructed; a well-supported call sits on
off-diagonal colocation signal.

Reads single-resolution .cool/.cooler and multi-resolution .mcool/.mcooler files (cooler
schema v2/v3, read directly with h5py), e.g. after `cooler zoomify`. File bins are summed
into larger plot bins (cooler's coarsen convention), so very large views stay cheap. Raw
counts by default; --balance applies bins/weight if present (from `cooler balance`).

The colour scale adapts to each view, as in Illumina's HiGlass colocation figures: the top sits at
the near-diagonal band (99.9th percentile of off-diagonal values, so the diagonal saturates) and the
log scale spans --dynamic-range (default 30x) below it, so sparse background counts stay white.

--compare puts a second map on the same grid: the main file above the diagonal, the second below
it, each with its own colour scale. Use it to set a long-read map (ont_colocation.py) against TruPath;
a junction seen by both technologies shows as mirror-image spots across the diagonal.

Examples:
  plot_colocation.py sample.colocation.cooler
  plot_colocation.py sample.colocation.cooler --region chr5:60,000,000-80,000,000
  plot_colocation.py sample.colocation.cooler --region chrX:150,000,000-156,000,000 --triangle --depth 2e6
  plot_colocation.py sample.colocation.mcool --region chr9 --region chr22 -o chr9_chr22.png
  plot_colocation.py sample.colocation.mcool --sv-vcf sample.sv.vcf.gz
  plot_colocation.py sample.colocation.cooler --compare ont.colocation.cool --names TruPath,ONT --region chr2:15,500,000-16,300,000
  plot_colocation.py sample.colocation.cooler --info

Requires: Python 3 with numpy, h5py, matplotlib.
"""
import argparse
import collections
import contextlib
import dataclasses
import gzip
import math
import re
import sys
import textwrap
from pathlib import Path

import numpy as np

# HiGlass's default "fall" colormap, so plots match the TruPath HiGlass view.
FALL = ["#ffffff", "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c",
        "#fc4e2a", "#e31a1c", "#bd0026", "#800026", "#000000"]
PRIMARY = re.compile(r"(chr)?(\d+|X|Y)", re.IGNORECASE)
NICE_STEPS = (1, 2, 2.5, 5)
MAX_PLOT_BINS = 5000          # dense n x n float64 matrix: 5000 bins = 200 MB
# SV overlay colour and marker per type. The colours avoid the heatmap's yellow-red range and
# stay distinct under colour-vision deficiency; the marker shape repeats the type.
SV_STYLE = {"DEL": ("#2a78d6", "o"), "DUP": ("#008300", "s"),
            "INV": ("#4a3aa7", "D"), "BND": ("#d55181", "^")}


class InputError(Exception):
    pass


@dataclasses.dataclass
class Segment:
    """A genomic interval on the plot axes (0-based half-open), plus its plot-bin layout."""
    label: str
    chrom: str
    start: int
    end: int
    offset: int = 0           # first plot bin of this segment
    n: int = 0                # number of plot bins
    lo: int = 0               # file-bin range [lo, hi) covering the interval
    hi: int = 0
    x0: int = 0               # genomic start of file bin lo, where plot bins begin


@dataclasses.dataclass
class SV:
    """One SV junction from the VCF: two breakpoints, 0-based."""
    id: str
    kind: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    filter: str

    @property
    def passed(self):
        return self.filter in ("PASS", ".")


def text(value):
    return value.decode() if isinstance(value, bytes) else value


class Cooler:
    """Minimal reader for one cooler group (cooler schema v2/v3, fixed-size bins)."""

    def __init__(self, h5, group):
        g = h5[group]
        missing = [k for k in ("chroms", "bins", "pixels", "indexes") if k not in g]
        if missing:
            raise InputError(f"{group!r} is not a cooler group (missing {', '.join(missing)}).")
        self.group = group
        self.attrs = {k: text(v) for k, v in g.attrs.items()}
        version = int(self.attrs.get("format-version", 0))
        if version < 2:
            raise InputError(f"Unsupported cooler format version {version}; need 2 or 3.")
        if self.attrs.get("bin-type", "fixed") != "fixed":
            raise InputError("Variable-size bins are not supported; DRAGEN writes fixed bins.")
        self.binsize = int(self.attrs["bin-size"])
        self.chroms = [text(n) for n in g["chroms/name"][:]]
        self.lengths = g["chroms/length"][:].astype(np.int64)
        self.chrom_offset = g["indexes/chrom_offset"][:].astype(np.int64)
        self.bin1_offset = g["indexes/bin1_offset"]
        self.bin1, self.bin2, self.count = (g["pixels/bin1_id"], g["pixels/bin2_id"],
                                            g["pixels/count"])
        self.weight = g["bins/weight"] if "weight" in g["bins"] else None
        # v2 has no storage-mode attribute and is always upper-triangular.
        self.symmetric = self.attrs.get("storage-mode", "symmetric-upper") == "symmetric-upper"
        self.nbins = int(self.chrom_offset[-1])
        self.nnz = int(self.bin1.shape[0])


def split_uri(uri):
    path, _, group = uri.partition("::")
    return path, group or None


def find_resolutions(h5, group):
    """Return {bin size: group path}: one entry for a .cool, all zoom levels for an .mcool."""
    if group is None and "resolutions" in h5:
        return {int(k): f"/resolutions/{k}" for k in h5["resolutions"]}
    group = group or "/"
    if group not in h5:
        raise InputError(f"Group {group!r} not found in file.")
    try:
        return {int(text(h5[group].attrs["bin-size"])): group}
    except (KeyError, ValueError):
        raise InputError(f"{group!r} has no fixed bin-size; is it a cooler group?") from None


def match_contig(name, chroms):
    if name in chroms:
        return name
    alias = name[3:] if name.startswith("chr") else "chr" + name
    if alias in chroms:
        return alias
    raise InputError(f"Contig {name!r} not found in the cooler.")


def parse_region(value, clr):
    m = re.fullmatch(r"(.+):([\d,]+)-([\d,]+)", value.strip())
    name = m.group(1) if m else value.strip()
    chrom = match_contig(name, clr.chroms)
    length = int(clr.lengths[clr.chroms.index(chrom)])
    if not m:
        return Segment(chrom, chrom, 0, length)
    start, end = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if not 1 <= start <= end <= length:
        raise InputError(f"Invalid region {value!r}: need 1 <= start <= end <= {length:,}.")
    return Segment(f"{chrom}:{start:,}-{end:,}", chrom, start - 1, end)


def genome_segments(clr, all_contigs):
    names = [c for c in clr.chroms if all_contigs or PRIMARY.fullmatch(c)] or clr.chroms
    return [Segment(c, c, 0, int(clr.lengths[clr.chroms.index(c)])) for c in names]


def choose_binsize(span, base, max_bins, requested):
    """Plot bin size: a multiple of the file's base bins; a round number when that is close."""
    if requested:
        if requested % base:
            raise InputError(f"--binsize must be a multiple of the file's {base:,} bp bins.")
        return requested
    target = span / max_bins
    plain = math.ceil(target / base) * base
    if plain <= base:
        return base
    scale = 10 ** math.floor(math.log10(target))
    for size in (int(step * s) for s in (scale, scale * 10) for step in NICE_STEPS):
        if size >= target and size % base == 0:
            return size if size <= 2 * plain else plain
    return plain


def layout(clr, segments, factor):
    """Place segments on the plot axes; return the file-bin -> plot-bin map and the axis size."""
    fmap = np.full(clr.nbins, -1, dtype=np.int64)
    n = 0
    for seg in segments:
        c = clr.chroms.index(seg.chrom)
        off = int(clr.chrom_offset[c])
        # Plot bins sit on a fixed genomic grid (multiples of the plot bin size), so two files with
        # different file bin sizes line up bin for bin.
        seg.lo = off + seg.start // (factor * clr.binsize) * factor
        seg.hi = off + -(-seg.end // clr.binsize)
        seg.x0 = (seg.lo - off) * clr.binsize
        if (fmap[seg.lo:seg.hi] >= 0).any():
            raise InputError(f"Region {seg.label} overlaps another region.")
        seg.offset, seg.n = n, -(-(seg.hi - seg.lo) // factor)
        fmap[seg.lo:seg.hi] = n + np.arange(seg.hi - seg.lo) // factor
        n += seg.n
    if n > MAX_PLOT_BINS:
        raise InputError(f"The view would be {n:,} bins across (limit {MAX_PLOT_BINS:,}); "
                         "use a larger --binsize or a smaller region.")
    return fmap, n


def accumulate(clr, segments, factor, chunksize, balance):
    """Sum stored pixels into a dense, symmetric plot-bin matrix over the concatenated segments.

    Returns (matrix, total, within): total sums each stored bin pair once; within is the part
    whose two bins fall in the same segment (intra-chromosomal in the genome-wide view).
    """
    fmap, n = layout(clr, segments, factor)
    seg_of = np.repeat(np.arange(len(segments)), [s.n for s in segments])
    weight = None
    if balance:
        if clr.weight is None:
            raise InputError("No bins/weight in this cooler; run `cooler balance` first, or drop --balance.")
        weight = clr.weight[:]
    upper = np.zeros(n * n)
    total = within = 0.0
    # Pixels are sorted by bin1, then bin2. For a symmetric-upper file, any stored pixel with
    # both bins in view has its bin1 in view, so reading each segment's rows covers it.
    for seg in segments:
        p0, p1 = int(clr.bin1_offset[seg.lo]), int(clr.bin1_offset[seg.hi])
        for s in range(p0, p1, chunksize):
            e = min(s + chunksize, p1)
            b1, b2 = clr.bin1[s:e], clr.bin2[s:e]
            v = clr.count[s:e].astype(np.float64)
            if weight is not None:
                v *= weight[b1] * weight[b2]
            r, c = fmap[b1], fmap[b2]
            keep = (c >= 0) & np.isfinite(v)
            r, c, v = r[keep], c[keep], v[keep]
            if not v.size:
                continue
            total += v.sum()
            within += v[seg_of[r] == seg_of[c]].sum()
            # Rows within a chunk are contiguous, so bincount only spans those rows.
            base = int(r.min()) * n
            acc = np.bincount(r * n + c - base, weights=v)
            upper[base:base + acc.size] += acc
    m = upper.reshape(n, n)
    if clr.symmetric:
        # Mirror; as in `cooler coarsen`, a plot-bin diagonal holds each in-bin pair once.
        m = m + m.T - np.diag(np.diag(m))
    return m, total, within


def fmt_bp(bp):
    for unit, size in (("Mb", 1e6), ("kb", 1e3)):
        if bp >= size:
            return f"{bp / size:,.4g} {unit}"
    return f"{bp:,.0f} bp"


def bare(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def read_sv_vcf(path):
    """Read SV junctions: DEL/DUP/INV as (POS, END), BND as (POS, mate) with each mate pair once.

    Returns (svs, skipped); skipped counts records without a second breakpoint by type.
    """
    with open(path, "rb") as handle:
        opener = gzip.open if handle.read(2) == b"\x1f\x8b" else open   # also reads bgzip
    svs, seen, skipped = [], set(), collections.Counter()
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                skipped["malformed"] += 1
                continue
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in f[7].split(";"))
            kind = info.get("SVTYPE") or f[4].strip("<>").split(":")[0]
            try:
                pos = int(f[1]) - 1
                if kind == "BND":
                    mate = re.search(r"[\[\]](.+):(\d+)[\[\]]", f[4])
                    if not mate:
                        skipped["single breakend"] += 1
                        continue
                    chrom2, pos2 = mate.group(1), int(mate.group(2)) - 1
                elif kind in ("DEL", "DUP", "INV"):
                    chrom2 = f[0]
                    pos2 = (int(info["END"]) - 1 if "END" in info
                            else pos + abs(int(info["SVLEN"].split(",")[0])))
                else:
                    skipped[kind] += 1
                    continue
            except (KeyError, ValueError):
                skipped[f"{kind} without breakpoints"] += 1
                continue
            key = (kind, *sorted([(f[0], pos), (chrom2, pos2)]))
            if kind == "BND" and key in seen:      # the mate record of a pair already read
                continue
            seen.add(key)
            svs.append(SV(f[2], kind, f[0], pos, chrom2, pos2, f[6]))
    return svs, skipped


def place_svs(svs, segments, binsize, min_length, max_distance):
    """Upper-triangle plot coordinates (x >= y) for each call in view; counts of hidden calls.

    One segment: genomic coordinates. Several: plot-bin coordinates, as drawn by imshow.
    """
    def where(chrom, pos):
        for s in segments:
            if bare(s.chrom) == bare(chrom) and s.start <= pos < s.end:
                return pos if len(segments) == 1 else s.offset + (pos - s.x0) / binsize - 0.5
        return None

    placed, hidden = [], collections.Counter()
    for sv in svs:
        a, b = where(sv.chrom1, sv.pos1), where(sv.chrom2, sv.pos2)
        if a is None or b is None:
            hidden["outside the view"] += 1
        elif bare(sv.chrom1) == bare(sv.chrom2) and abs(sv.pos2 - sv.pos1) < min_length:
            hidden[f"shorter than {fmt_bp(min_length)}"] += 1
        elif max_distance and abs(sv.pos2 - sv.pos1) > max_distance:
            hidden["beyond --depth"] += 1
        else:
            placed.append((max(a, b), min(a, b), sv))
    return placed, hidden


def draw_svs(ax, placed, triangle):
    """Hollow markers with a white halo, so the signal under each call stays visible."""
    from matplotlib.lines import Line2D

    limits = ax.get_xlim(), ax.get_ylim()
    handles = []
    for kind, (color, marker) in SV_STYLE.items():
        calls = [(x, y, sv) for x, y, sv in placed if sv.kind == kind]
        for passed in (True, False):
            points = [((x + y) / 2, (x - y) / 2) if triangle else (x, y)
                      for x, y, sv in calls if sv.passed == passed]
            if not points:
                continue
            xs, ys = zip(*points)
            # Large, thin rings: the hollow centre stays clear even over a one-pixel spot.
            for edge, width, z in (("white", 2.6, 4), (color, 1.4, 5)):
                ax.scatter(xs, ys, s=150, marker=marker, facecolors="none", edgecolors=edge,
                           linewidths=width, alpha=1 if passed else 0.45, zorder=z)
        if calls:
            handles.append(Line2D([], [], ls="", marker=marker, markersize=7, markerfacecolor="none",
                                  markeredgecolor=color, markeredgewidth=1.6, label=f"{kind} ({len(calls)})"))
    if any(not sv.passed for *_, sv in placed):
        handles.append(Line2D([], [], ls="", marker="o", markersize=7, markerfacecolor="none",
                              markeredgecolor="#6b7280", markeredgewidth=1.6, alpha=0.45,
                              label="non-PASS (faded)"))
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    return handles


def print_info(path, h5, res):
    print(path)
    if len(res) > 1:
        print(f"  multi-resolution: {', '.join(fmt_bp(r) for r in sorted(res))}")
    clr = Cooler(h5, res[min(res)])
    a = clr.attrs
    print(f"  group {clr.group}: cooler format v{a.get('format-version')}, "
          f"{a.get('storage-mode', 'symmetric-upper')}, {fmt_bp(clr.binsize)} bins")
    print(f"  {len(clr.chroms):,} contigs, {clr.nbins:,} bins, {clr.nnz:,} non-zero pixels, "
          f"count dtype {clr.count.dtype}, balancing weights: {'yes' if clr.weight is not None else 'no'}")
    primary = [c for c in clr.chroms if PRIMARY.fullmatch(c)]
    print(f"  primary contigs ({len(primary)}): {' '.join(primary) or '-'}")
    for key in ("generated-by", "genome-assembly", "assembly", "creation-date"):
        if a.get(key):
            print(f"  {key}: {a[key]}")
    if a.get("metadata") and a["metadata"] not in ("{}", "null"):
        print(f"  metadata: {a['metadata']}")


def color_range(m, args, vmin=None, vmax=None, what="the selected view"):
    """Default colour limits. Top: the 99.9th percentile of non-zero off-diagonal values, i.e. the
    near-diagonal band, so the diagonal itself saturates and a few extreme bins do not set the scale.
    Bottom: on a log scale --dynamic-range below the top (values under it are drawn white); on a
    linear scale zero. vmin/vmax override either end."""
    positive = m[m > 0]
    if not positive.size:
        raise InputError(f"No colocation counts in {what}.")
    off = m[~np.eye(len(m), dtype=bool)]
    off = off[off > 0]
    if vmax is None:
        vmax = float(np.quantile(off if off.size else positive, 0.999))
    if vmin is None:
        vmin = 0.0 if args.linear else vmax / args.dynamic_range
    if not vmin < vmax:
        raise InputError(f"Colour scale minimum {vmin:g} is not below the maximum {vmax:g}; adjust --vmin/--vmax.")
    return vmin, vmax


def plot(m, segments, binsize, file_res, args, total, within, out, genome, placed=None, other=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as effects
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    cmap = (LinearSegmentedColormap.from_list("fall", FALL) if args.cmap == "fall"
            else matplotlib.colormaps[args.cmap]).with_extremes(bad="white", under="white")
    vmin, vmax = color_range(m, args, args.vmin, args.vmax)
    norm = Normalize(vmin, vmax) if args.linear else LogNorm(vmin, vmax)
    shown = np.ma.masked_less_equal(m, 0)        # empty bin pairs stay white
    image, draw = shown, dict(cmap=cmap, norm=norm)
    scale = "linear" if args.linear else "log"
    if other is not None:
        # Split view: main file on and above the diagonal, the second file below it, each on its
        # own adaptive scale (the two technologies count different things).
        names = other["names"]
        vmin_b, vmax_b = color_range(other["m"], args, what=f"{names[1]} for this view")
        norm_b = Normalize(vmin_b, vmax_b) if args.linear else LogNorm(vmin_b, vmax_b)
        upper = np.triu(np.ones(m.shape, dtype=bool))[..., None]
        image = np.where(upper, cmap(norm(shown)), cmap(norm_b(np.ma.masked_less_equal(other["m"], 0))))
        draw = {}

    name = Path(split_uri(args.cooler)[0]).name
    view = ("genome-wide" if genome else " + ".join(s.label for s in segments))
    kind = "balanced" if args.balance else "raw counts"
    if other is None:
        title = args.title or f"Colocation map · {name}"
        caption = (f"{view} · {fmt_bp(binsize)} bins (file {fmt_bp(file_res)}) · {kind}, "
                   f"{scale} scale {vmin:.3g}–{vmax:.3g}")
    else:
        title = args.title or f"Colocation maps · {names[0]} vs {names[1]}"
        caption = f"{view} · {fmt_bp(binsize)} bins · {kind}"
    if len(segments) > 1:
        caption += f" · {within / total:.1%} {'intra-chromosomal' if genome else 'within regions'}"
    if other is not None:
        caption += (f"\nupper: {names[0]}, {scale} {vmin:.3g}–{vmax:.3g} · "
                    f"lower: {names[1]}, {scale} {vmin_b:.3g}–{vmax_b:.3g}")
    # Wrap so a long caption cannot widen the figure past the plot.
    caption = "\n".join(textwrap.fill(line, 105) for line in caption.split("\n"))

    seg = segments[0]
    span = seg.end - seg.start
    unit, unit_name = (1e6, "Mb") if span >= 2e6 else (1e3, "kb")
    bp_axis = FuncFormatter(lambda x, _: f"{x / unit:,.6g}")
    x0 = seg.x0
    x1 = min(x0 + seg.n * binsize, seg.end) if len(segments) == 1 else None

    if args.triangle:
        depth = min(args.depth or span, x1 - x0)
        edges = np.minimum(x0 + np.arange(seg.n + 1) * binsize, x1).astype(float)
        X, Y = np.meshgrid(edges, edges)                 # X: column (bin2), Y: row (bin1)
        i, j = np.indices(m.shape)
        cells = np.ma.masked_where((j < i) | (m <= 0), m)
        fig, ax = plt.subplots(figsize=(10, max(2.4, min(6.0, 10 * depth / 2 / (x1 - x0)) + 1.0)))
        mappable = ax.pcolormesh((X + Y) / 2, (X - Y) / 2, cells, cmap=cmap, norm=norm,
                                 shading="flat", rasterized=True)
        ax.set_xlim(x0, x1)
        ax.set_ylim(0, depth / 2)
        ax.set_aspect("equal")
        ax.xaxis.set_major_formatter(bp_axis)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{2 * y / unit:,.6g}"))
        ax.set_xlabel(f"{seg.chrom} position ({unit_name})")
        ax.set_ylabel(f"distance ({unit_name})")
    else:
        fig, ax = plt.subplots(figsize=(8.6, 7.6))
        if len(segments) == 1:
            mappable = ax.imshow(image, interpolation="nearest", **draw,
                                 extent=(x0, x0 + seg.n * binsize, x0 + seg.n * binsize, x0))
            ax.set_xlim(x0, x1)
            ax.set_ylim(x1, x0)
            ax.xaxis.set_major_formatter(bp_axis)
            ax.yaxis.set_major_formatter(bp_axis)
            ax.set_xlabel(f"{seg.chrom} position ({unit_name})")
            ax.set_ylabel(f"{seg.chrom} position ({unit_name})")
        else:
            mappable = ax.imshow(image, interpolation="nearest", **draw)
            n = m.shape[0]
            for s in segments[1:]:
                for line in (ax.axhline, ax.axvline):
                    line(s.offset - 0.5, color="#6b7280", lw=0.5, alpha=0.7)
            # Label segments wide enough to hold text; the genome view drops "chr" to fit.
            labelled = [s for s in segments if s.n >= 0.008 * n]
            ticks = [s.offset + s.n / 2 - 0.5 for s in labelled]
            labels = [re.sub(r"^chr", "", s.label) if genome else s.label for s in labelled]
            ax.set_xticks(ticks, labels)
            ax.set_yticks(ticks, labels)
            ax.tick_params(length=0)
            if not genome:
                plt.setp(ax.get_yticklabels(), rotation=90, va="center")
        ax.spines[["top", "right"]].set_visible(True)
    # Inset colorbars follow the axes box, which aspect="equal" shrinks in the triangle view. In the
    # split view the top half belongs to the upper triangle and the bottom half to the lower one.
    left = 1.02 if args.triangle else 1.04
    label = "balanced colocation" if args.balance else "colocation count"
    bars = [(mappable, [left, 0, 0.018, 1], label)] if other is None else [
        (ScalarMappable(norm, cmap), [left, 0.53, 0.018, 0.47], f"{names[0]} (upper)"),
        (ScalarMappable(norm_b, cmap), [left, 0, 0.018, 0.47], f"{names[1]} (lower)")]
    for source, rect, text in bars:
        bar = fig.colorbar(source, cax=ax.inset_axes(rect))
        bar.set_label(text, fontsize=9)
        if not args.linear:                          # label 20, 30, 50, 70, 100, ... as Illumina does
            bar.ax.yaxis.set_major_locator(LogLocator(subs=(1, 2, 3, 5, 7)))
            bar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            bar.ax.yaxis.set_minor_formatter(NullFormatter())
        bar.ax.tick_params(labelsize=8)
        bar.outline.set_visible(False)
    if other is not None:
        for (tx, ty), ha, va, text in (((0.985, 0.985), "right", "top", names[0]),
                                       ((0.015, 0.015), "left", "bottom", names[1])):
            corner = ax.text(tx, ty, text, transform=ax.transAxes, ha=ha, va=va, fontsize=10,
                             fontweight="bold", color="#111827")
            corner.set_path_effects([effects.withStroke(linewidth=3, foreground="white")])
    ax.tick_params(labelsize=7 if genome else 8)
    ax.xaxis.label.set_fontsize(9)
    ax.yaxis.label.set_fontsize(9)
    ax.set_title(title + "\n" * (caption.count("\n") + 1), fontsize=12, fontweight="bold", loc="left")
    ax.text(0, 1.02, caption, transform=ax.transAxes, fontsize=8.5, color="#4b5563",
            ha="left", va="bottom")
    if placed is not None:
        handles = draw_svs(ax, placed, args.triangle)
        # Legend in a row under the axes, below the tick and axis labels.
        fig.canvas.draw()
        bottom = ax.get_tightbbox(fig.canvas.get_renderer()).y0
        y = ax.transAxes.inverted().transform((0, bottom))[1]
        legend_title = f"SV calls ({'all filters' if args.sv_all else 'PASS'}) · {Path(args.sv_vcf).name}"
        if not handles:
            legend_title += ": none in view"
        ax.legend(handles=handles, title=legend_title, loc="upper left", bbox_to_anchor=(0, y - 0.01),
                  ncol=max(1, len(handles)), frameon=False, fontsize=8, title_fontsize=8,
                  handletextpad=0.2, columnspacing=1.4, borderaxespad=0, alignment="left")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def file_stem(path):
    stem = Path(split_uri(path)[0]).name
    for ext in (".mcooler", ".mcool", ".cooler", ".cool"):
        if stem.endswith(ext):
            return stem[: -len(ext)]
    return stem


def default_output(path, segments, genome, triangle, compare=None):
    stem = file_stem(path) + (f".vs.{file_stem(compare).removesuffix('.colocation')}" if compare else "")
    view = "genome" if genome else "__".join(
        re.sub(r"[^\w.-]+", "_", s.label.replace(",", "")).strip("_") for s in segments)
    return f"{stem}.{view}{'.triangle' if triangle else ''}.png"


def read_view(h5, res, segments, binsize, args):
    """Sum one file into plot bins, from its coarsest zoom level that divides the plot bin size."""
    file_res = max(r for r in res if binsize % r == 0)
    clr = Cooler(h5, res[file_res])
    m, total, within = accumulate(clr, segments, binsize // file_res, args.chunksize, args.balance)
    return file_res, clr, m, total, within


def same_grid(segments, a, b, names):
    """The main file's segments on the second file's contig names. Both maps must share the
    reference build, so every contig in view must have the same length in both files."""
    out = []
    for seg in segments:
        try:
            chrom = match_contig(seg.chrom, b.chroms)
        except InputError:
            raise InputError(f"Contig {seg.chrom} is not in {names[1]}.") from None
        la, lb = int(a.lengths[a.chroms.index(seg.chrom)]), int(b.lengths[b.chroms.index(chrom)])
        if la != lb:
            raise InputError(f"{seg.chrom} is {la:,} bp in {names[0]} but {lb:,} bp in {names[1]}; "
                             "both maps need the same reference build.")
        out.append(dataclasses.replace(seg, chrom=chrom))
    return out


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cooler", help="Colocation cooler (.cool/.cooler), multi-resolution .mcool/.mcooler, "
                                  "or a URI such as sample.mcool::/resolutions/2000.")
    p.add_argument("-r", "--region", action="append", default=[],
                   help="chrom or chrom:start-end (1-based, inclusive). Repeat to place several regions "
                        "side by side. Default: genome-wide.")
    p.add_argument("--all-contigs", action="store_true",
                   help="Genome-wide view: also include alt/decoy/unplaced contigs and chrM "
                        "(default: chr1-22, X, Y).")
    p.add_argument("--binsize", type=int,
                   help="Plot bin size in bp, a multiple of the file's bin size. Default: chosen from --max-bins.")
    p.add_argument("--max-bins", type=int, default=1000,
                   help="Automatic bin size keeps the view at most this many bins across (default 1000, so every "
                        "bin gets at least one pixel at the default --dpi and a single-bin SV spot is never dropped).")
    p.add_argument("--triangle", action="store_true",
                   help="One region only: draw the upper triangle rotated 45 degrees, like a HiGlass horizontal heatmap.")
    p.add_argument("--depth", type=float,
                   help="With --triangle: largest genomic distance shown, in bp (default: the whole region).")
    p.add_argument("--sv-vcf", help="DRAGEN SV VCF (<sample>.sv.vcf.gz) whose calls are marked on the map.")
    p.add_argument("--sv-all", action="store_true", help="Also mark non-PASS calls (drawn faded).")
    p.add_argument("--sv-types", default="DEL,DUP,INV,BND",
                   help="Comma-separated SV types to mark (default DEL,DUP,INV,BND).")
    p.add_argument("--sv-min-length", type=int,
                   help="Hide intra-chromosomal calls shorter than this, in bp; they sit on the diagonal "
                        "(default: two plot bins).")
    p.add_argument("--compare", help="Second colocation cooler (e.g. from ont_colocation.py), drawn below the "
                                     "diagonal on the same grid with its own colour scale.")
    p.add_argument("--names", help="Names of the main and --compare maps, e.g. TruPath,ONT (default: file names).")
    p.add_argument("--balance", action="store_true", help="Use bins/weight (`cooler balance`) instead of raw counts.")
    p.add_argument("--linear", action="store_true", help="Linear colour scale (default: log).")
    p.add_argument("--vmax", type=float,
                   help="Colour scale maximum (default: 99.9th percentile of off-diagonal values; the diagonal "
                        "saturates). With --compare it applies to the main map; the second stays adaptive.")
    p.add_argument("--vmin", type=float,
                   help="Colour scale minimum; values below it are white (default: --vmax / --dynamic-range).")
    p.add_argument("--dynamic-range", type=float, default=30,
                   help="Log scale: how many times smaller the minimum is than the maximum (default 30).")
    p.add_argument("--cmap", default="fall", help="'fall' (HiGlass default) or any matplotlib colormap name.")
    p.add_argument("-o", "--out", help="Output image; format from the extension (.png/.pdf/.svg). "
                                       "Default: <cooler name>.<view>.png")
    p.add_argument("--title", help="Plot title (default: 'Colocation map · <file name>').")
    p.add_argument("--dpi", type=int, default=200, help="Output resolution (default 200).")
    p.add_argument("--chunksize", type=int, default=5_000_000,
                   help="Pixels read per chunk; lower it to save memory (default 5,000,000).")
    p.add_argument("--info", action="store_true", help="Print the cooler's metadata and exit.")
    args = p.parse_args(argv)
    if args.triangle and len(args.region) != 1:
        p.error("--triangle needs exactly one --region.")
    if args.compare and args.triangle:
        p.error("--compare draws a square split view; it cannot be combined with --triangle.")
    if args.names and not args.compare:
        p.error("--names is used with --compare.")
    if args.compare:
        names = args.names.split(",") if args.names else [file_stem(args.cooler), file_stem(args.compare)]
        if len(names) != 2 or not all(n.strip() for n in names):
            p.error("--names needs two comma-separated names, e.g. --names TruPath,ONT.")
        args.names = [n.strip() for n in names]
    if args.depth is not None and (not args.triangle or not args.depth > 0):
        p.error("--depth must be positive and is used with --triangle.")
    if any(v is not None and v < 1 for v in (args.binsize, args.max_bins, args.dpi, args.chunksize)):
        p.error("--binsize, --max-bins, --dpi and --chunksize must be positive.")
    args.sv_types = {t.strip().upper() for t in args.sv_types.split(",") if t.strip()}
    if not args.sv_types <= set(SV_STYLE):
        p.error(f"--sv-types must be a subset of {','.join(SV_STYLE)}.")
    if not args.sv_vcf and (args.sv_all or args.sv_min_length is not None):
        p.error("--sv-all and --sv-min-length are used with --sv-vcf.")
    if args.sv_min_length is not None and args.sv_min_length < 0:
        p.error("--sv-min-length must be zero or positive.")
    if not args.dynamic_range > 1:
        p.error("--dynamic-range must be greater than 1.")
    if not args.linear and args.vmin is not None and args.vmin <= 0:
        p.error("--vmin must be positive on a log scale (or add --linear).")
    if args.cmap != "fall":
        import matplotlib
        if args.cmap not in matplotlib.colormaps:
            p.error(f"Unknown colormap {args.cmap!r}.")
    return args


def svs_in_view(args, segments, binsize):
    """Read, filter and place the SV calls, and report what is shown and what is hidden."""
    svs, skipped = read_sv_vcf(args.sv_vcf)
    hidden = collections.Counter()
    kept = []
    for sv in svs:
        if sv.kind not in args.sv_types:
            hidden["type not in --sv-types"] += 1
        elif not (sv.passed or args.sv_all):
            hidden["non-PASS"] += 1
        else:
            kept.append(sv)
    min_length = 2 * binsize if args.sv_min_length is None else args.sv_min_length
    placed, out_of_view = place_svs(kept, segments, binsize, min_length,
                                    args.depth if args.triangle else None)
    hidden.update(out_of_view)
    kinds = collections.Counter(sv.kind for *_, sv in placed)
    print(f"SV calls from {args.sv_vcf}: {len(placed)} marked"
          + (f" ({', '.join(f'{k} {kinds[k]}' for k in SV_STYLE if kinds[k])})" if placed else ""))
    for reason, count in hidden.items():
        print(f"  hidden: {count} {reason}")
    if skipped:
        print(f"  skipped (no second breakpoint): {', '.join(f'{n} {k}' for k, n in skipped.items())}")
    for x, y, sv in sorted(placed, key=lambda p: (p[1], p[0]))[:25]:
        print(f"  {sv.kind:<4} {sv.filter:<12} {sv.chrom1}:{sv.pos1 + 1:,} - {sv.chrom2}:{sv.pos2 + 1:,}  {sv.id}")
    if len(placed) > 25:
        print(f"  ... and {len(placed) - 25} more")
    return placed


def main(argv=None):
    args = parse_args(argv)
    try:
        import h5py
    except ImportError:
        sys.exit("error: h5py is required (conda install h5py, or pip install h5py).")
    path, group = split_uri(args.cooler)
    cpath, cgroup = split_uri(args.compare) if args.compare else (None, None)
    try:
        for name in (path, args.sv_vcf, cpath):
            if name and not Path(name).is_file():
                raise InputError(f"File not found: {name}")
        with contextlib.ExitStack() as files:
            h5 = files.enter_context(h5py.File(path, "r"))
            res = find_resolutions(h5, group)
            if args.info:
                print_info(path, h5, res)
                return 0
            finest = Cooler(h5, res[min(res)])
            genome = not args.region
            segments = (genome_segments(finest, args.all_contigs) if genome
                        else [parse_region(r, finest) for r in args.region])
            span = sum(s.end - s.start for s in segments)
            base = finest.binsize
            if args.compare:
                h5b = files.enter_context(h5py.File(cpath, "r"))
                res_b = find_resolutions(h5b, cgroup)
                finest_b = Cooler(h5b, res_b[min(res_b)])
                base = math.lcm(base, finest_b.binsize)          # plot bins both files divide into
            binsize = choose_binsize(span, base, args.max_bins, args.binsize)
            file_res, clr, m, total, within = read_view(h5, res, segments, binsize, args)
            other = None
            if args.compare:
                segments_b = same_grid(segments, finest, finest_b, args.names)
                file_res_b, clr_b, m_b, total_b, _ = read_view(h5b, res_b, segments_b, binsize, args)
                other = dict(m=m_b, names=args.names)
        out = args.out or default_output(path, segments, genome, args.triangle, args.compare)
        print(f"{path}{'::' + clr.group if clr.group != '/' else ''}: "
              f"{fmt_bp(file_res)} bins, {clr.nnz:,} non-zero pixels")
        if other:
            print(f"compare {cpath}{'::' + clr_b.group if clr_b.group != '/' else ''}: "
                  f"{fmt_bp(file_res_b)} bins, {clr_b.nnz:,} non-zero pixels; "
                  f"{total_b:,.6g} {'weight' if args.balance else 'counts'} in view")
        print(f"view: {'genome-wide' if genome else ' + '.join(s.label for s in segments)} "
              f"({len(segments)} segment{'s' * (len(segments) > 1)}, {fmt_bp(span)}) "
              f"at {fmt_bp(binsize)} bins -> {m.shape[0]:,} x {m.shape[0]:,}")
        if total:
            share = f" ({within / total:.1%} within one {'chromosome' if genome else 'region'})" \
                if len(segments) > 1 else ""
            print(f"colocation {'weight' if args.balance else 'counts'} in view: {total:,.6g}{share}")
        placed = svs_in_view(args, segments, binsize) if args.sv_vcf else None
        plot(m, segments, binsize, file_res, args, total, within, out, genome, placed, other)
        print(f"wrote {out}")
        return 0
    except (InputError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
