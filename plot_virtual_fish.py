#!/usr/bin/env python3
"""
plot_virtual_fish.py - "virtual FISH": paint chromosomes or genomic regions onto a TruPath flow-cell lane.

Every Illumina read name records the nanowell the read was sequenced in
(INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y). This script colours each read pair by the chromosome or
region it maps to and draws it where it sat on the flow cell, like a chromosome-paint FISH image of
the lane. On TruPath data the reads of one long DNA molecule land in neighbouring nanowells, so a
paint shows up as single-coloured spots: the molecules' constellations.

One figure, three zoom levels:
  lane     every tile of one lane, placed by tile number (surface, swath, tile); colour is the paint
           mix and intensity the read density, from a uniform sample when there are many reads
  tile     one tile at full density
  window   a few thousand read-name units of that tile, one dot per read pair; reads with the same
           paint that are close on the flow cell and in the genome are joined into constellations

With two or more paints, constellations of different paints that sit together (one DNA molecule
carrying both, e.g. a BCR-ABL1 fusion) are counted lane-wide, like a dual-fusion FISH probe, next to
the count expected by chance from a tile-shifted control.

Paints (--paint, repeatable, up to 8): a chromosome (chr7), a region (chr9:130,700,000-130,900,000) or
a labelled region (ABL1=chr9:130,713,000-130,887,000). Without --paint every primary chromosome
(1-22, X, Y) gets its own colour; that reads the whole file.

Examples:
  plot_virtual_fish.py --input sample.cram --reference genome.fa --paint chr7
  plot_virtual_fish.py --input sample.cram --reference genome.fa \\
      --paint ABL1=chr9:130,713,000-130,887,000 --paint BCR=chr22:23,180,000-23,320,000
  plot_virtual_fish.py --input sample.bam --lane 2 --tile 1205        # all chromosomes

Read pairs are counted once (read 1, or unpaired reads); secondary, supplementary, QC-fail and
duplicate records are skipped. Requires samtools on PATH and Python 3 with numpy and matplotlib.
"""
import argparse
import array
import collections
import dataclasses
import math
import re
import shutil
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path

import numpy as np

# Paint colours in fixed order (validated categorical palette; the first three are the most
# distinct, so put the paints that matter most first).
PAINT_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
PRIMARY = re.compile(r"(chr)?(\d+|X|Y)", re.IGNORECASE)
FLUSH = 1 << 20


class InputError(Exception):
    pass


@dataclasses.dataclass
class Paint:
    label: str
    chrom: str
    start: int                # 0-based half-open
    end: int
    whole: bool

    @property
    def region(self):
        name = f"{{{self.chrom}}}" if ":" in self.chrom else self.chrom
        return name if self.whole else f"{name}:{self.start + 1}-{self.end}"


def run_capture(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode:
        raise InputError(f"samtools failed: {result.stderr.strip()}")
    return result.stdout


def read_contigs(args):
    cmd = [args.samtools, "view", "-H"] + (["-T", args.reference] if args.reference else []) + [args.input]
    contigs = {}
    for line in run_capture(cmd).splitlines():
        if line.startswith("@SQ"):
            tags = dict(t.split(":", 1) for t in line.split("\t")[1:] if ":" in t)
            contigs[tags["SN"]] = int(tags["LN"])
    if not contigs:
        raise InputError("No @SQ lines in the alignment header; is the file aligned?")
    return contigs


def match_contig(name, contigs):
    if name in contigs:
        return name
    alias = name[3:] if name.startswith("chr") else "chr" + name
    if alias in contigs:
        return alias
    raise InputError(f"Contig {name!r} not found in the alignment header.")


def parse_paint(value, contigs):
    label, sep, spec = value.partition("=")
    if not sep:
        label, spec = "", value
    m = re.fullmatch(r"(.+):([\d,]+)-([\d,]+)", spec.strip())
    chrom = match_contig(m.group(1) if m else spec.strip(), contigs)
    if not m:
        return Paint(label or chrom, chrom, 0, contigs[chrom], True)
    start, end = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if not 1 <= start <= end <= contigs[chrom]:
        raise InputError(f"Invalid paint {value!r}: need 1 <= start <= end <= {contigs[chrom]:,}.")
    return Paint(label or f"{chrom}:{start:,}-{end:,}", chrom, start - 1, end, False)


def check_overlaps(paints):
    for a in paints:
        for b in paints:
            if a is not b and a.chrom == b.chrom and a.start < b.end and b.start < a.end:
                raise InputError(f"Paints {a.label} and {b.label} overlap; a read can carry only one paint.")


def tone(hex_color, toward, amount):
    rgb = np.array([int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    return tuple(rgb + (toward - rgb) * amount)


def paint_colors(n):
    """RGB per paint. Up to 8: the palette. 24 chromosomes: the 8 hues in base, dark and light tones;
    the window panel's constellation labels name them, so identity never rests on colour alone."""
    if n <= len(PAINT_COLORS):
        return np.array([tone(c, 0, 0) for c in PAINT_COLORS[:n]])
    tones = [(0, 0), (0, 0.38), (1, 0.32)]
    return np.array([tone(PAINT_COLORS[i % 8], *tones[i // 8 % 3]) for i in range(n)])


class Sample:
    """Reads kept for plotting: every read on the protected (zoom) tile, plus a uniform sample of at
    most `budget` other reads. Reads are kept when a hash of their name falls below a threshold that
    halves each time the budget is exceeded, so the sample stays uniform while it is thinned."""

    def __init__(self, budget):
        self.budget, self.threshold, self.protected = budget, 1 << 32, None
        self.buffer = [array.array(code) for code in "iiihiI"]     # tile, x, y, paint, pos, hash
        self.chunks, self.sampled = [], 0

    def add(self, tile, x, y, paint, pos, digest):
        protected = tile == self.protected
        if digest >= self.threshold and not protected:
            return
        b = self.buffer
        b[0].append(tile)
        b[1].append(x)
        b[2].append(y)
        b[3].append(paint)
        b[4].append(pos)
        b[5].append(digest)
        if not protected:
            self.sampled += 1
            if self.sampled > self.budget:
                self.thin()
        if len(self.buffer[0]) >= FLUSH:
            self.flush()

    def flush(self):
        if len(self.buffer[0]):
            self.chunks.append([np.frombuffer(c, dtype=c.typecode).copy() for c in self.buffer])
            self.buffer = [array.array(c.typecode) for c in self.buffer]

    def thin(self):
        self.flush()
        self.threshold >>= 1
        kept, self.sampled = [], 0
        for chunk in self.chunks:
            outside = chunk[0] != self.protected
            keep = (chunk[5] < self.threshold) | ~outside
            kept.append([c[keep] for c in chunk])
            self.sampled += int((keep & outside).sum())
        self.chunks = kept

    def arrays(self):
        self.flush()
        if not self.chunks:
            return [np.zeros(0, dtype=c.typecode) for c in self.buffer]
        return [np.concatenate(parts) for parts in zip(*self.chunks)]


def collect(args, paints):
    """Stream read 1 records from samtools and keep those of the chosen lane that carry a paint.

    Lines are handled as bytes and samtools drops the aux tags; for CRAM, htslib decodes only the
    fields used here (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR), skipping sequence and qualities.
    """
    flags = 4 | 128 | 256 | 512 | 2048 | (0 if args.keep_duplicates else 1024)
    cmd = [args.samtools, "view", "-M", "-F", str(flags), "-q", str(args.mapq), "-@", str(args.threads),
           "--keep-tag", "RG"]
    if args.reference:
        cmd += ["-T", args.reference]
    if Path(args.input).suffix.lower() == ".cram":
        cmd += ["--input-fmt-option", "required_fields=0x3f"]
    if args.subsample:
        cmd += ["--subsample", str(args.subsample), "--subsample-seed", "1"]
    cmd += [args.input] + [p.region for p in paints]
    targets = collections.defaultdict(list)
    for i, p in enumerate(paints):
        targets[p.chrom.encode()].append((p.start, p.end, i))
    want_fc, _, want_lane = (v.encode() for v in args.lane.rpartition(":")) if args.lane else (b"", b"", None)
    want_tile = args.tile.encode() if args.tile else None
    lane, lanes, tiles = None, collections.Counter(), {}
    totals, bad = [0] * len(paints), 0
    sample = Sample(args.max_points)
    with tempfile.TemporaryFile(mode="w+t") as errors:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errors, bufsize=1 << 20)
        try:
            for line in proc.stdout:
                qname, _, rname, pos, _ = line.split(b"\t", 4)
                f = qname.split(b":")
                if len(f) < 7:
                    bad += 1
                    continue
                lanes[f[2], f[3]] += 1
                if lane is None:
                    if (want_lane and f[3] != want_lane) or (want_fc and f[2] != want_fc):
                        continue
                    lane = (f[2], f[3])
                if f[2] != lane[0] or f[3] != lane[1]:
                    continue
                p = int(pos) - 1
                for start, end, paint in targets.get(rname, ()):
                    if start <= p < end:
                        break
                else:
                    continue
                try:
                    x, y = int(f[5]), int(f[6])
                except ValueError:
                    bad += 1
                    continue
                tile = tiles.get(f[4])
                if tile is None:
                    tile = tiles[f[4]] = len(tiles)
                    if sample.protected is None and f[4] == (want_tile or f[4]):
                        sample.protected = tile
                totals[paint] += 1
                sample.add(tile, x, y, paint, p, zlib.crc32(qname))
            proc.stdout.close()
            code = proc.wait()
        except BaseException:
            proc.terminate()
            proc.wait()
            raise
        errors.seek(0)
        diagnostic = errors.read().strip()
    if code:
        raise InputError(f"samtools view failed: {diagnostic}")
    if lane is None:
        if bad and not lanes:
            raise InputError("Read names lack INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y; were the original "
                             "Illumina names kept?")
        seen = ", ".join(f"{fc.decode()}:{ln.decode()}" for fc, ln in sorted(lanes)) or "none"
        raise InputError(f"No painted reads in lane {args.lane or '(any)'}; lanes seen: {seen}.")
    lane = (lane[0].decode(), lane[1].decode())
    if want_tile and want_tile not in tiles:
        raise InputError(f"Tile {args.tile} has no painted reads in lane {lane[0]}:{lane[1]}.")
    names = [n.decode() for n in sorted(tiles, key=tiles.get)]
    lanes = collections.Counter({(fc.decode(), ln.decode()): n for (fc, ln), n in lanes.items()})
    return dict(lane=lane, lanes=lanes, tiles=names, totals=totals, bad=bad, sample=sample,
                columns=sample.arrays())


def tile_layout(names):
    """(surface, swath, number) per tile from SSTT names such as 1205; else name order, 20 per row."""
    if all(len(n) == 4 and n.isdigit() for n in names):
        return [(int(n[0]), int(n[1]), int(n[2:])) for n in names], True
    rank = {n: i for i, n in enumerate(sorted(names, key=lambda n: (len(n), n)))}
    return [(1, rank[n] // 20 + 1, rank[n] % 20 + 1) for n in names], False


def densest(x, y, size, bounds):
    """Centre of the size x size square (on a size/3 grid) holding the most reads."""
    (x0, x1), (y0, y1) = bounds
    step = size / 3
    hx = np.floor((x - x0) / step).astype(int)
    hy = np.floor((y - y0) / step).astype(int)
    grid = np.zeros((hx.max() + 3, hy.max() + 3))
    np.add.at(grid, (hx, hy), 1)
    window = sum(np.roll(np.roll(grid, -i, 0), -j, 1) for i in range(3) for j in range(3))
    i, j = np.unravel_index(window.argmax(), window.shape)
    cx = min(max(x0 + (i + 1.5) * step, x0 + size / 2), x1 - size / 2)
    cy = min(max(y0 + (j + 1.5) * step, y0 + size / 2), y1 - size / 2)
    return cx, cy


def constellations(x, y, paint, pos, radius, max_gap):
    """Single-linkage groups of reads with the same paint, within `radius` on the flow cell and
    `max_gap` bp in the genome. Returns a group label per read."""
    n = len(x)
    labels = np.arange(n)
    if n < 2:
        return labels
    cx, cy = np.floor(x / radius).astype(np.int64), np.floor(y / radius).astype(np.int64)
    key = (cx - cx.min()) * (cy.max() - cy.min() + 3) + (cy - cy.min())
    stride = int(cy.max() - cy.min() + 3)
    order = np.argsort(key, kind="stable")
    cells, starts = np.unique(key[order], return_index=True)
    bounds = dict(zip(cells.tolist(), zip(starts.tolist(), np.append(starts[1:], n).tolist())))
    left, right = [], []
    for cell, (s, e) in bounds.items():
        a = order[s:e]
        for offset in (0, 1, stride - 1, stride, stride + 1):        # each neighbouring cell pair once
            if cell + offset not in bounds:
                continue
            s2, e2 = bounds[cell + offset]
            b = order[s2:e2]
            ok = (((x[a, None] - x[b]) ** 2 + (y[a, None] - y[b]) ** 2 <= radius ** 2)
                  & (paint[a, None] == paint[b]) & (np.abs(pos[a, None] - pos[b]) <= max_gap))
            if offset == 0:
                ok &= a[:, None] < b
            i, j = np.nonzero(ok)
            left.append(a[i])
            right.append(b[j])
    i, j = np.concatenate(left), np.concatenate(right)
    while True:                                   # label propagation with pointer jumping
        low = np.minimum(labels[i], labels[j])
        new = labels.copy()
        np.minimum.at(new, i, low)
        np.minimum.at(new, j, low)
        new = new[new]
        if np.array_equal(new, labels):
            return labels
        labels = new


def centres(x, y, paint, groups, min_reads):
    """(x, y, paint, reads) of each constellation with at least min_reads reads."""
    ids, inverse, sizes = np.unique(groups, return_inverse=True, return_counts=True)
    keep = sizes >= min_reads
    cx = np.bincount(inverse, weights=x)[keep] / sizes[keep]
    cy = np.bincount(inverse, weights=y)[keep] / sizes[keep]
    cp = np.zeros(len(ids), dtype=int)
    cp[inverse] = paint                               # one paint per constellation by construction
    return cx, cy, cp[keep], sizes[keep]


def pair_up(cx, cy, cp, reach):
    """Greedy pairs (i, j) of constellations with different paints whose centres lie within reach."""
    used, pairs = set(), []
    for i in range(len(cx)):
        if i in used:
            continue
        d = np.hypot(cx - cx[i], cy - cy[i])
        for j in np.argsort(d):
            if d[j] > reach:
                break
            if j != i and j not in used and cp[j] != cp[i]:
                used |= {i, int(j)}
                pairs.append((i, int(j)))
                break
    return pairs


def colocalization(data, paints, args):
    """Lane-wide counts of constellation pairs with different paints on one footprint (centres within
    2 x --link-radius), per paint pair, with a control that shifts one paint's constellations to the
    next tiles. Needs every read, so it is skipped when the lane was sampled."""
    tile, x, y, paint, pos, _ = data["columns"]
    if not 2 <= len(paints) <= len(PAINT_COLORS) or data["sample"].threshold < 1 << 32 or len(x) > 1_000_000:
        return None
    spacing = float(x.max() - x.min()) + 10 * args.link_radius        # keeps tiles apart
    groups = constellations(x + tile * spacing, y, paint, pos, args.link_radius, args.link_gap)
    cx, cy, cp, _ = centres(x + tile * spacing, y, paint, groups, args.min_constellation)
    ct = np.floor(cx / spacing).astype(int)
    lx = cx - ct * spacing
    ntiles, reach = len(data["tiles"]), 2 * args.link_radius
    results = []
    for a in range(len(paints)):
        for b in range(a + 1, len(paints)):
            def count(shift):
                n = 0
                for t in np.unique(ct[cp == a]):
                    pa = (ct == t) & (cp == a)
                    pb = ((ct + shift) % ntiles == t) & (cp == b)
                    sel = pa | pb
                    n += len(pair_up(lx[sel], cy[sel], cp[sel], reach))
                return n
            shifts = [s for s in range(1, min(ntiles, 6))]
            control = np.mean([count(s) for s in shifts]) if shifts else float("nan")
            results.append(dict(a=paints[a].label, b=paints[b].label, observed=count(0), expected=control,
                                n_a=int((cp == a).sum()), n_b=int((cp == b).sum())))
    return results


def paint_raster(rows, cols, paint, shape, colors):
    """Blend paint colours per pixel; intensity follows read density, with a floor so lone reads show."""
    npx = shape[0] * shape[1]
    lin = rows * shape[1] + cols
    count = np.bincount(lin, minlength=npx).astype(float)
    image = np.ones((npx, 3))
    filled = count > 0
    if filled.any():
        mean = np.stack([np.bincount(lin, weights=colors[paint, k], minlength=npx)[filled]
                         for k in range(3)], axis=1) / count[filled, None]
        alpha = 0.35 + 0.65 * np.clip(count[filled] / np.quantile(count[filled], 0.9), 0, 1)
        image[filled] = 1 - alpha[:, None] * (1 - mean)
    return image.reshape(shape[0], shape[1], 3)


def plot(data, paints, args, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as effects
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    tile, x, y, paint, pos, digest = data["columns"]
    colors = paint_colors(len(paints))
    names, sample = data["tiles"], data["sample"]
    zoom = data["zoom_tile"]
    x0, x1, y0, y1 = int(x.min()), int(x.max()) + 1, int(y.min()), int(y.max()) + 1
    xr, yr = x1 - x0, y1 - y0
    ink, muted = "#111827", "#4b5563"

    # --- lane overview raster: Y runs along the lane, X across the swath -------------------------
    layout, physical = tile_layout(names)
    surfaces = sorted({s for s, _, _ in layout})
    swaths = sorted({w for _, w, _ in layout})
    ncols = max(num for _, _, num in layout)
    cw = int(min(40, max(4, 2400 // ncols)))
    ch = int(max(3, round(cw * xr / yr)))
    gap = max(2, ch // 2)
    rows = [(s, w) for s in surfaces for w in swaths]
    top = {sw: r * ch + surfaces.index(sw[0]) * gap for r, sw in enumerate(rows)}
    shape = (len(rows) * ch + (len(surfaces) - 1) * gap, ncols * cw)
    tile_top = np.array([top[(s, w)] for s, w, _ in layout])
    tile_left = np.array([(num - 1) * cw for _, _, num in layout])
    overview = digest.astype(np.int64) < sample.threshold
    t = tile[overview]
    cols = tile_left[t] + np.clip(((y[overview] - y0) / yr * cw).astype(int), 0, cw - 1)
    rws = tile_top[t] + np.clip(((x1 - 1 - x[overview]) / xr * ch).astype(int), 0, ch - 1)
    lane_image = paint_raster(rws, cols, paint[overview], shape, colors)

    # --- zoom tile and window ------------------------------------------------------------------
    on_tile = tile == zoom
    tx, ty, tp, tpos = x[on_tile], y[on_tile], paint[on_tile], pos[on_tile]
    dots = len(tx) <= 300_000          # sparse tiles as dots; full real tiles as a density raster
    if not dots:
        tw = 900
        th = max(50, int(round(tw * xr / yr)))
        tile_image = paint_raster(np.clip(((x1 - 1 - tx) / xr * th).astype(int), 0, th - 1),
                                  np.clip(((ty - y0) / yr * tw).astype(int), 0, tw - 1), tp, (th, tw), colors)
    size = args.window_size
    cx, cy = data["window_center"]
    wx0, wy0 = cx - size / 2, cy - size / 2
    inside = (tx >= wx0) & (tx < wx0 + size) & (ty >= wy0) & (ty < wy0 + size)
    wx, wy, wp, wpos = tx[inside], ty[inside], tp[inside], tpos[inside]
    groups = (constellations(wx, wy, wp, wpos, args.link_radius, args.link_gap)
              if len(wx) <= 200_000 else np.arange(len(wx)))
    group_ids, group_sizes = np.unique(groups, return_counts=True)
    big = group_ids[group_sizes >= args.min_constellation]
    in_big = np.isin(groups, big)
    data["window"] = dict(reads=len(wx), constellations=len(big), share=in_big.mean() if len(wx) else 0)

    # --- figure ---------------------------------------------------------------------------------
    lane_h = 13.2 * shape[0] / shape[1] + 0.5
    row_h = max(4.2, min(7.5, 8.2 * xr / yr))
    ncol = min(len(paints), 4 if len(paints) <= 8 else 8)
    header = 0.95 + 0.24 * math.ceil(len(paints) / ncol) + (0.3 if data.get("colocalization") else 0)
    fig_h = header + lane_h + row_h + 1.2
    fig = plt.figure(figsize=(14, fig_h))
    grid = fig.add_gridspec(2, 2, height_ratios=[lane_h, row_h], width_ratios=[1.65, 1],
                            hspace=0.9 / ((lane_h + row_h) / 2), wspace=0.12, left=0.06, right=0.98,
                            top=1 - header / fig_h, bottom=0.6 / fig_h)
    lane_ax = fig.add_subplot(grid[0, :])
    lane_ax.imshow(lane_image, interpolation="nearest", aspect="equal")
    for c in range(1, ncols):
        lane_ax.axvline(c * cw - 0.5, color="#e5e7eb", lw=0.4)
    for sw in rows:
        lane_ax.axhline(top[sw] - 0.5, color="#d1d5db", lw=0.5)
        lane_ax.axhline(top[sw] + ch - 0.5, color="#d1d5db", lw=0.5)
    zs, zw, zn = layout[zoom]
    lane_ax.add_patch(Rectangle(((zn - 1) * cw - 0.5, top[(zs, zw)] - 0.5), cw, ch,
                                fill=False, ec=ink, lw=1.4))
    step = 10 if ncols > 30 else 5 if ncols > 10 else 1
    ticks = [n for n in range(1, ncols + 1) if n == 1 or n % step == 0]
    lane_ax.set_xticks([(n - 0.5) * cw for n in ticks], [str(n) for n in ticks])
    lane_ax.set_yticks([top[sw] + ch / 2 for sw in rows], [f"{s}·{w}" for s, w in rows])
    lane_ax.set_xlabel("tile number" if physical else "tile (name order)", fontsize=9)
    lane_ax.set_ylabel("surface·swath" if physical else "row", fontsize=9)
    lane_ax.tick_params(labelsize=7, length=0)
    for side in lane_ax.spines.values():
        side.set_visible(False)
    share = sample.threshold / (1 << 32)
    lane_ax.set_title(f"Lane {data['lane'][0]}:{data['lane'][1]} · {len(names)} tiles · colour = paint mix, "
                      f"intensity = read density"
                      + ("" if share >= 1 else f" · uniform {share:.3g} sample"),
                      fontsize=9, color=muted, loc="left")
    if data.get("colocalization"):
        fig.text(0.06, 1 - 0.7 / fig_h - 0.24 * math.ceil(len(paints) / ncol) / fig_h,
                 "Co-localized constellations (one footprint, both paints): " + "; ".join(
                     f"{r['a']}+{r['b']} {r['observed']:,} (≈{r['expected']:.1f} by chance)"
                     for r in sorted(data["colocalization"], key=lambda r: r["expected"] - r["observed"])[:3])
                 + ("; more in the console output" if len(data["colocalization"]) > 3 else ""),
                 fontsize=9, color=ink, ha="left", va="top")

    tile_ax = fig.add_subplot(grid[1, 0])
    if dots:
        shuffle = np.random.default_rng(1).permutation(len(tx))
        tile_ax.scatter(ty[shuffle], tx[shuffle], c=colors[tp[shuffle]], linewidths=0, rasterized=True,
                        s=float(np.clip(4e4 / max(len(tx), 1), 0.3, 5)))
        tile_ax.set_xlim(y0, y1)
        tile_ax.set_ylim(x0, x1)
        tile_ax.set_aspect("equal")
    else:
        tile_ax.imshow(tile_image, interpolation="nearest", extent=(y0, y1, x0, x1), aspect="equal")
    tile_ax.add_patch(Rectangle((wy0, wx0), size, size, fill=False, ec=ink, lw=1.2))
    tile_ax.set_title(f"Tile {names[zoom]} · all {len(tx):,} painted read pairs", fontsize=9,
                      color=muted, loc="left")
    win_ax = fig.add_subplot(grid[1, 1])
    order = np.random.default_rng(0).permutation(len(wx))
    win_ax.scatter(wy[order], wx[order], c=colors[wp[order]], s=10, linewidths=0, rasterized=True)
    # Label each constellation (with two or more paints); two of different paints on one footprint
    # get one joint label.
    gx, gy, gp, _ = centres(wx, wy, wp, groups, args.min_constellation)
    if len(paints) == 1:
        gx = gx[:0]
    pairs = pair_up(gx, gy, gp, 2 * args.link_radius) if len(paints) <= 8 else []
    paired = {i for pair in pairs for i in pair}
    labels = [(gx[[i, j]].mean(), gy[[i, j]].mean(),
               "+".join(paints[k].label for k in sorted((gp[i], gp[j])))) for i, j in pairs]
    labels += [(gx[i], gy[i], paints[gp[i]].label) for i in range(len(gx)) if i not in paired]
    for lx, ly, label in labels:
        text = win_ax.text(ly, lx, re.sub(r"(^|\+)chr", r"\1", label), fontsize=7, color=ink,
                           ha="center", va="center")
        text.set_path_effects([effects.withStroke(linewidth=2.4, foreground="white")])
    win_ax.set_xlim(wy0, wy0 + size)
    win_ax.set_ylim(wx0, wx0 + size)
    win_ax.set_aspect("equal")
    win = data["window"]
    win_ax.set_title(f"Window · {win['constellations']} constellations (≥{args.min_constellation} reads) "
                     f"hold {win['share']:.0%} of {win['reads']:,}", fontsize=9, color=muted, loc="left")
    for ax in (tile_ax, win_ax):
        ax.ticklabel_format(style="plain", useOffset=False)
        ax.tick_params(labelsize=7)
        ax.set_xlabel("Y · read-name units (along the lane)", fontsize=8)
    tile_ax.set_ylabel("X · read-name units", fontsize=8)

    title = args.title or f"Virtual FISH · {Path(args.input).name}"
    fig.text(0.06, 1 - 0.1 / fig_h, title, fontsize=14, fontweight="bold", ha="left", va="top")
    total = sum(data["totals"])
    fig.text(0.06, 1 - 0.45 / fig_h, f"{total:,} painted read pairs in lane {data['lane'][0]}:{data['lane'][1]} · "
                                     f"MAPQ ≥ {args.mapq} · each dot is one read pair at its nanowell",
             fontsize=9, color=muted, ha="left", va="top")
    handles = [Line2D([], [], ls="", marker="o", markersize=6, markerfacecolor=colors[i],
                      markeredgewidth=0, label=f"{p.label} ({data['totals'][i]:,})")
               for i, p in enumerate(paints)]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.055, 1 - 0.7 / fig_h), frameon=False,
               ncol=ncol, fontsize=8, handletextpad=0.2, columnspacing=1.2, borderaxespad=0)
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def choose_zoom(data, args):
    """Zoom tile: --tile, else the protected first tile if the sample was thinned, else the tile with
    the densest spot. Window: --window, else the densest spot on that tile."""
    tile, x, y = data["columns"][:3]
    sample = data["sample"]
    bounds = ((x.min(), x.max() + 1), (y.min(), y.max() + 1))
    size = args.window_size
    if args.tile or sample.threshold < 1 << 32:
        zoom = sample.protected
    else:
        best = -1
        for t in np.unique(tile):
            on = tile == t
            cx, cy = densest(x[on], y[on], size, bounds)
            n = int(((np.abs(x[on] - cx) < size / 2) & (np.abs(y[on] - cy) < size / 2)).sum())
            if n > best:
                best, zoom = n, int(t)
    on = tile == zoom
    if args.window:
        cx, cy = args.window
        (xa, xb), (ya, yb) = bounds
        if not (xa <= cx <= xb and ya <= cy <= yb):
            raise InputError(f"--window {cx:g},{cy:g} lies outside the tiles (X {xa:,}-{xb:,}, Y {ya:,}-{yb:,}).")
    else:
        cx, cy = densest(x[on], y[on], size, bounds)
    data["zoom_tile"], data["window_center"] = zoom, (cx, cy)


def default_output(args, paints, genome):
    stem = Path(args.input).name
    for ext in (".cram", ".bam", ".sam"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    tag = "all" if genome else "_".join(re.sub(r"[^\w.-]+", "_", p.label.replace(",", "")) for p in paints)
    return f"{stem}.fish.{tag}.png"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "--cram", "--bam", dest="input", required=True,
                   help="Indexed, coordinate-sorted BAM/CRAM with original Illumina read names.")
    p.add_argument("--reference", "-T", help="Reference FASTA; required for CRAM.")
    p.add_argument("--paint", action="append", default=[],
                   help="chrom, chrom:start-end (1-based, inclusive) or LABEL=chrom:start-end; repeat for up "
                        "to 8 paints. Default: every primary chromosome.")
    p.add_argument("--lane", help="Lane to draw: a lane number or FLOWCELL:LANE. Default: the first lane in the file.")
    p.add_argument("--tile", help="Tile to zoom into, e.g. 1205. Default: the densest spot when every read "
                                  "fits in --max-points, otherwise the first tile with a painted read.")
    p.add_argument("--window", help="Window centre on the zoom tile as X,Y (read-name units). Default: the densest spot.")
    p.add_argument("--window-size", type=float, default=3000,
                   help="Window width in read-name units (default 3000).")
    p.add_argument("--link-radius", type=float, default=350,
                   help="Constellation link distance on the flow cell, read-name units (default 350).")
    p.add_argument("--link-gap", type=int, default=500_000,
                   help="Constellation link distance in the genome, bp (default 500,000).")
    p.add_argument("--min-constellation", type=int, default=4,
                   help="Reads needed to label a constellation in the window (default 4).")
    p.add_argument("--mapq", type=int, default=20, help="Minimum MAPQ (default 20).")
    p.add_argument("--keep-duplicates", action="store_true", help="Keep reads flagged as duplicates.")
    p.add_argument("--max-points", type=int, default=4_000_000,
                   help="Reads kept for the lane overview; beyond this a uniform sample is drawn. The zoom "
                        "tile always keeps every read (default 4,000,000).")
    p.add_argument("--subsample", type=float,
                   help="Pass only this fraction of reads through samtools, for a quick look (thins the zoom too).")
    p.add_argument("--threads", type=int, default=2, help="samtools decompression threads (default 2).")
    p.add_argument("-o", "--out", help="Output image (.png/.pdf/.svg). Default: <input>.fish.<paints>.png")
    p.add_argument("--title", help="Plot title (default: 'Virtual FISH · <input name>').")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--samtools", default="samtools", help="Path/name of the samtools executable.")
    args = p.parse_args(argv)
    if Path(args.input).suffix.lower() == ".cram" and not args.reference:
        p.error("CRAM input requires --reference.")
    if len(args.paint) > len(PAINT_COLORS):
        p.error(f"At most {len(PAINT_COLORS)} paints; omit --paint to paint every chromosome.")
    if args.window:
        try:
            wx, wy = args.window.split(",")
            args.window = (float(wx), float(wy))
        except ValueError:
            p.error("--window must be X,Y without thousands separators, e.g. --window 12000,30500.")
    if not 0 <= args.mapq <= 255:
        p.error("--mapq must be 0..255.")
    if args.subsample is not None and not 0 < args.subsample < 1:
        p.error("--subsample must be between 0 and 1.")
    if min(args.window_size, args.link_radius) <= 0 or min(args.max_points, args.min_constellation,
                                                          args.threads, args.dpi) < 1 or args.link_gap < 0:
        p.error("Sizes, distances and counts must be positive.")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        if not shutil.which(args.samtools):
            raise InputError(f"samtools executable not found: {args.samtools}")
        for name in (args.input, args.reference):
            if name and not Path(name).is_file():
                raise InputError(f"File not found: {name}")
        contigs = read_contigs(args)
        genome = not args.paint
        if genome:
            paints = [Paint(c, c, 0, n, True) for c, n in contigs.items() if PRIMARY.fullmatch(c)]
            if not paints:
                raise InputError("No primary chromosomes (1-22, X, Y) in the header; use --paint.")
        else:
            paints = [parse_paint(v, contigs) for v in args.paint]
            check_overlaps(paints)
            if len({p.label for p in paints}) < len(paints):
                raise InputError("Paint labels must be distinct.")
        print("paints: " + ", ".join(p.label if p.whole and p.label == p.chrom
                                     else f"{p.label}={p.chrom}:{p.start + 1:,}-{p.end:,}" for p in paints))
        data = collect(args, paints)
        choose_zoom(data, args)
        data["colocalization"] = colocalization(data, paints, args)
        out = args.out or default_output(args, paints, genome)
        plot(data, paints, args, out)
        lane = f"{data['lane'][0]}:{data['lane'][1]}"
        others = [f"{fc}:{ln} ({n:,})" for (fc, ln), n in sorted(data["lanes"].items()) if f"{fc}:{ln}" != lane]
        print(f"lane {lane}: {sum(data['totals']):,} painted read pairs on {len(data['tiles'])} tiles"
              + (f"; other lanes in the file: {', '.join(others)}" if others else ""))
        share = data["sample"].threshold / (1 << 32)
        if share < 1:
            print(f"lane overview drawn from a uniform {share:.3g} sample (--max-points {args.max_points:,})")
        if data["bad"]:
            print(f"skipped {data['bad']:,} reads whose names lack flow-cell coordinates", file=sys.stderr)
        for r in data["colocalization"] or []:
            print(f"co-localized {r['a']}+{r['b']}: {r['observed']:,} constellation pairs on one footprint "
                  f"(≈{r['expected']:.1f} expected by chance; {r['n_a']:,} {r['a']} and {r['n_b']:,} {r['b']} "
                  f"constellations of ≥{args.min_constellation} reads)")
        if 2 <= len(paints) <= len(PAINT_COLORS) and data["colocalization"] is None:
            print("co-localization not counted: it needs every read of the lane (raise --max-points, "
                  "or paint smaller regions)")
        w = data["window"]
        cx, cy = data["window_center"]
        print(f"zoom: tile {data['tiles'][data['zoom_tile']]}, window centre X={cx:,.0f} Y={cy:,.0f}: "
              f"{w['reads']:,} read pairs, {w['constellations']} constellations of ≥{args.min_constellation} "
              f"reads holding {w['share']:.0%} of them")
        print(f"wrote {out}")
        return 0
    except (InputError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
