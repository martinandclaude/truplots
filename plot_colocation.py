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

Reads single-resolution .cool/.cooler and multi-resolution .mcool/.mcooler files (cooler
schema v2/v3, read directly with h5py), e.g. after `cooler zoomify`. File bins are summed
into larger plot bins (cooler's coarsen convention), so very large views stay cheap. Raw
counts by default; --balance applies bins/weight if present (from `cooler balance`).

Examples:
  plot_colocation.py sample.colocation.cooler
  plot_colocation.py sample.colocation.cooler --region chr5:60,000,000-80,000,000
  plot_colocation.py sample.colocation.cooler --region chrX:150,000,000-156,000,000 --triangle --depth 2e6
  plot_colocation.py sample.colocation.mcool --region chr9 --region chr22 -o chr9_chr22.png
  plot_colocation.py sample.colocation.cooler --info

Requires: Python 3 with numpy, h5py, matplotlib.
"""
import argparse
import dataclasses
import math
import re
import sys
from pathlib import Path

import numpy as np

# HiGlass's default "fall" colormap, so plots match the TruPath HiGlass view.
FALL = ["#ffffff", "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c",
        "#fc4e2a", "#e31a1c", "#bd0026", "#800026", "#000000"]
PRIMARY = re.compile(r"(chr)?(\d+|X|Y)", re.IGNORECASE)
NICE_STEPS = (1, 2, 2.5, 5)
MAX_PLOT_BINS = 5000          # dense n x n float64 matrix: 5000 bins = 200 MB


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
        seg.lo = off + seg.start // clr.binsize
        seg.hi = off + -(-seg.end // clr.binsize)
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


def plot(m, segments, binsize, file_res, args, total, within, out, genome):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
    from matplotlib.ticker import FuncFormatter

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    cmap = (LinearSegmentedColormap.from_list("fall", FALL) if args.cmap == "fall"
            else matplotlib.colormaps[args.cmap]).with_extremes(bad="white")
    positive = m[m > 0]
    if not positive.size:
        raise InputError("No colocation counts in the selected view.")
    vmin = args.vmin if args.vmin is not None else positive.min()
    vmax = args.vmax if args.vmax is not None else positive.max()
    norm = Normalize(vmin, vmax) if args.linear else LogNorm(vmin, vmax)
    shown = np.ma.masked_less_equal(m, 0)        # empty bin pairs stay white

    name = Path(split_uri(args.cooler)[0]).name
    view = ("genome-wide" if genome else " + ".join(s.label for s in segments))
    title = args.title or f"Colocation map · {name}"
    caption = (f"{view} · {fmt_bp(binsize)} bins (file {fmt_bp(file_res)}) · "
               f"{'balanced' if args.balance else 'raw counts'}, "
               f"{'linear' if args.linear else 'log'} scale")
    if len(segments) > 1:
        caption += f" · {within / total:.1%} {'intra-chromosomal' if genome else 'within regions'}"

    seg = segments[0]
    span = seg.end - seg.start
    unit, unit_name = (1e6, "Mb") if span >= 2e6 else (1e3, "kb")
    bp_axis = FuncFormatter(lambda x, _: f"{x / unit:,.6g}")
    x0 = seg.start - seg.start % file_res        # plot bins start on a file-bin boundary
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
            mappable = ax.imshow(shown, cmap=cmap, norm=norm, interpolation="nearest",
                                 extent=(x0, x0 + seg.n * binsize, x0 + seg.n * binsize, x0))
            ax.set_xlim(x0, x1)
            ax.set_ylim(x1, x0)
            ax.xaxis.set_major_formatter(bp_axis)
            ax.yaxis.set_major_formatter(bp_axis)
            ax.set_xlabel(f"{seg.chrom} position ({unit_name})")
            ax.set_ylabel(f"{seg.chrom} position ({unit_name})")
        else:
            mappable = ax.imshow(shown, cmap=cmap, norm=norm, interpolation="nearest")
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
    # An inset colorbar follows the axes box, which aspect="equal" shrinks in the triangle view.
    bar = fig.colorbar(mappable, cax=ax.inset_axes([1.02 if args.triangle else 1.04, 0, 0.018, 1]))
    bar.set_label("balanced colocation" if args.balance else "colocation count", fontsize=9)
    bar.ax.tick_params(labelsize=8)
    bar.outline.set_visible(False)
    ax.tick_params(labelsize=7 if genome else 8)
    ax.xaxis.label.set_fontsize(9)
    ax.yaxis.label.set_fontsize(9)
    ax.set_title(f"{title}\n", fontsize=12, fontweight="bold", loc="left")
    ax.text(0, 1.02, caption, transform=ax.transAxes, fontsize=8.5, color="#4b5563",
            ha="left", va="bottom")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def default_output(path, segments, genome, triangle):
    stem = Path(path).name
    for ext in (".mcooler", ".mcool", ".cooler", ".cool"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    view = "genome" if genome else "__".join(
        re.sub(r"[^\w.-]+", "_", s.label.replace(",", "")).strip("_") for s in segments)
    return f"{stem}.{view}{'.triangle' if triangle else ''}.png"


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
    p.add_argument("--max-bins", type=int, default=1500,
                   help="Automatic bin size keeps the view at most this many bins across (default 1500).")
    p.add_argument("--triangle", action="store_true",
                   help="One region only: draw the upper triangle rotated 45 degrees, like a HiGlass horizontal heatmap.")
    p.add_argument("--depth", type=float,
                   help="With --triangle: largest genomic distance shown, in bp (default: the whole region).")
    p.add_argument("--balance", action="store_true", help="Use bins/weight (`cooler balance`) instead of raw counts.")
    p.add_argument("--linear", action="store_true", help="Linear colour scale (default: log).")
    p.add_argument("--vmin", type=float, help="Colour scale minimum (default: smallest non-zero value).")
    p.add_argument("--vmax", type=float, help="Colour scale maximum (default: largest value).")
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
    if args.depth is not None and (not args.triangle or not args.depth > 0):
        p.error("--depth must be positive and is used with --triangle.")
    if any(v is not None and v < 1 for v in (args.binsize, args.max_bins, args.dpi, args.chunksize)):
        p.error("--binsize, --max-bins, --dpi and --chunksize must be positive.")
    if not args.linear and args.vmin is not None and args.vmin <= 0:
        p.error("--vmin must be positive on a log scale (or add --linear).")
    if args.cmap != "fall":
        import matplotlib
        if args.cmap not in matplotlib.colormaps:
            p.error(f"Unknown colormap {args.cmap!r}.")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        import h5py
    except ImportError:
        sys.exit("error: h5py is required (conda install h5py, or pip install h5py).")
    path, group = split_uri(args.cooler)
    try:
        if not Path(path).is_file():
            raise InputError(f"File not found: {path}")
        with h5py.File(path, "r") as h5:
            res = find_resolutions(h5, group)
            if args.info:
                print_info(path, h5, res)
                return 0
            finest = Cooler(h5, res[min(res)])
            genome = not args.region
            segments = (genome_segments(finest, args.all_contigs) if genome
                        else [parse_region(r, finest) for r in args.region])
            span = sum(s.end - s.start for s in segments)
            binsize = choose_binsize(span, finest.binsize, args.max_bins, args.binsize)
            # From an .mcool, read the coarsest zoom level that divides the plot bin size.
            file_res = max(r for r in res if binsize % r == 0)
            clr = finest if file_res == finest.binsize else Cooler(h5, res[file_res])
            m, total, within = accumulate(clr, segments, binsize // file_res, args.chunksize, args.balance)
        out = args.out or default_output(path, segments, genome, args.triangle)
        print(f"{path}{'::' + clr.group if clr.group != '/' else ''}: "
              f"{fmt_bp(file_res)} bins, {clr.nnz:,} non-zero pixels")
        print(f"view: {'genome-wide' if genome else ' + '.join(s.label for s in segments)} "
              f"({len(segments)} segment{'s' * (len(segments) > 1)}, {fmt_bp(span)}) "
              f"at {fmt_bp(binsize)} bins -> {m.shape[0]:,} x {m.shape[0]:,}")
        if total:
            share = f" ({within / total:.1%} within one {'chromosome' if genome else 'region'})" \
                if len(segments) > 1 else ""
            print(f"colocation {'weight' if args.balance else 'counts'} in view: {total:,.6g}{share}")
        plot(m, segments, binsize, file_res, args, total, within, out, genome)
        print(f"wrote {out}")
        return 0
    except (InputError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
