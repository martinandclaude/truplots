#!/usr/bin/env python3
"""Paint one or two genomic probes on TruPath flowcell coordinates.

Requires Python >=3.10, samtools and matplotlib. Regions are 1-based inclusive;
internal intervals are 0-based half-open. No clustering uses genomic distance.
DRAGEN 4.6 HP/pp and optional HZ/pz are retained separately, as are PS and ps.
BX groups are inferred templates; spatial groups are exploratory candidates.
"""
from __future__ import annotations

import argparse
import collections
import csv
import dataclasses
import gzip
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import unquote

VERSION = "0.1.0"
TAGS = ("BX", "HP", "PS", "pp", "HZ", "pz", "PC", "li", "xq", "ps", "tr",
        "sd", "hl", "gs", "hs", "js", "lv")
COLORS = {0: "#9ca3af", 1: "#dc3545", 2: "#159b66", 3: "#9b59b6"}


class InputError(Exception):
    pass


@dataclasses.dataclass(frozen=True)
class Probe:
    label: str
    chrom: str
    start: int
    end: int
    origin: str = "region"

    @property
    def region(self):
        return f"{self.chrom}:{self.start + 1}-{self.end}"


@dataclasses.dataclass(frozen=True, order=True)
class Location:
    instrument: str
    run: str
    flowcell: str
    lane: int
    tile: str
    x: int
    y: int

    @property
    def tile_key(self):
        return (self.instrument, self.run, self.flowcell, self.lane, self.tile)


@dataclasses.dataclass
class Alignment:
    chrom: str
    start: int
    end: int
    cigar: str
    flag: int
    mapq: int
    blocks: list[tuple[int, int]]
    tags: dict[str, str]
    tag_types: dict[str, str]
    targets: int = 0


@dataclasses.dataclass
class Point:
    qname: str
    rg: str
    location: Location
    alignments: list[Alignment] = dataclasses.field(default_factory=list)
    targets: int = 0
    primary_targets: int = 0
    supplementary_targets: int = 0
    bx: str | None = None
    bx_conflict: bool = False
    group: str | None = None


def parse_region(value):
    label, sep, region = value.partition("=")
    if not sep:
        region, label = value, value
    try:
        chrom, coords = region.rsplit(":", 1)
        lo, hi = coords.replace(",", "").split("-", 1)
        start, end = int(lo), int(hi)
    except (ValueError, TypeError):
        raise InputError(f"Invalid region {value!r}; use [LABEL=]chr:start-end.") from None
    if not label or not chrom or start < 1 or end < start:
        raise InputError(f"Invalid region {value!r}; coordinates must satisfy 1 <= start <= end.")
    return Probe(label, chrom, start - 1, end)


def parse_location(qname):
    fields = qname.split(":")
    if len(fields) < 7 or not all(fields[:5]):
        raise ValueError("QNAME lacks instrument:run:flowcell:lane:tile:x:y")
    # Validate every value before adding the record to any collection.
    lane, x, y = int(fields[3]), int(fields[5]), int(fields[6])
    if lane < 1:
        raise ValueError("Invalid lane")
    return Location(*fields[:3], lane, fields[4], x, y)


def cigar_blocks(start, cigar):
    """Return actual aligned-base blocks and reference-consuming end.

    D and N advance reference coordinates but do not create probe support.
    An insertion or soft clip alone also does not support a genomic probe.
    """
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not ops or "".join(n + op for n, op in ops) != cigar:
        raise ValueError("Invalid CIGAR")
    pos, blocks = start, []
    for length, op in ops:
        n = int(length)
        if n < 1:
            raise ValueError("Zero-length CIGAR operation")
        if op in "M=X":
            if blocks and blocks[-1][1] == pos:
                blocks[-1] = (blocks[-1][0], pos + n)
            else:
                blocks.append((pos, pos + n))
            pos += n
        elif op in "DN":
            pos += n
    return blocks, pos


def parse_sam_record(line):
    f = line.rstrip("\n").split("\t")
    if len(f) < 11:
        raise ValueError("SAM record has fewer than 11 fields")
    flag, start, mapq = int(f[1]), int(f[3]) - 1, int(f[4])
    if flag & 4 or f[2] == "*" or start < 0:
        raise ValueError("Unmapped or positionless record")
    tags, types = {}, {}
    for token in f[11:]:
        parts = token.split(":", 2)
        if len(parts) == 3:
            key, typ, value = parts
            if key in TAGS or key == "RG":
                tags[key], types[key] = value, typ
    blocks, end = cigar_blocks(start, f[5])
    return f[0], tags.get("RG", ""), parse_location(f[0]), Alignment(
        f[2], start, end, f[5], flag, mapq, blocks, tags, types)


def annotation_attributes(value):
    if '"' in value:
        return dict(re.findall(r'(\S+)\s+"([^"]*)"', value))
    return {k: unquote(v) for item in value.strip().split(";")
            if "=" in item for k, v in [item.split("=", 1)]}


def resolve_genes(path, queries):
    """Resolve exact gene symbols/IDs; SYMBOL@contig disambiguates alternatives."""
    requests = [q.rsplit("@", 1) if "@" in q else (q, None) for q in queries]
    candidates = [dict() for _ in queries]
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                raise InputError(f"Annotation line {number} has {len(f)} fields, expected 9.")
            attr = annotation_attributes(f[8])
            is_gene = f[2].lower() in {"gene", "pseudogene"}
            # GTF fallback can aggregate transcripts/exons by gene_id. GFF3 needs gene features.
            symbol = attr.get("gene_name") or attr.get("gene") or (attr.get("Name") if is_gene else None)
            gid = attr.get("gene_id") or (attr.get("ID") if is_gene else None)
            names = {n for n in (symbol, gid) if n}
            for i, (query, contig) in enumerate(requests):
                if query not in names or (contig and f[0] != contig):
                    continue
                try:
                    lo, hi = int(f[3]) - 1, int(f[4])
                except ValueError:
                    raise InputError(f"Invalid annotation coordinates at line {number}.") from None
                if lo < 0 or hi <= lo:
                    raise InputError(f"Invalid annotation coordinates at line {number}.")
                key = (gid or symbol, f[0])
                bucket = candidates[i].setdefault(key, {"gene": [], "fallback": []})
                bucket["gene" if is_gene else "fallback"].append((lo, hi))
    result = []
    for query, matches in zip(queries, candidates):
        if not matches:
            raise InputError(f"Gene {query!r} was not found in {path}; use an exact gene_name/gene_id/Name/ID.")
        if len(matches) != 1:
            choices = ", ".join(f"{gid}@{chrom}" for gid, chrom in sorted(matches))
            raise InputError(f"Gene {query!r} is ambiguous: {choices}. Use a gene ID or SYMBOL@contig.")
        (_, chrom), bounds = next(iter(matches.items()))
        intervals = bounds["gene"] or bounds["fallback"]
        result.append(Probe(query, chrom, min(a for a, _ in intervals),
                            max(b for _, b in intervals), "gene"))
    return result


def run_capture(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode:
        raise InputError(f"samtools failed: {result.stderr.strip()}")
    return result.stdout


def parse_header(header):
    contigs, read_groups, programs = {}, [], []
    for line in header.splitlines():
        f = line.split("\t")
        attributes = dict(x.split(":", 1) for x in f[1:] if ":" in x)
        if f[0] == "@SQ":
            contigs[attributes["SN"]] = int(attributes["LN"])
        elif f[0] == "@RG":
            read_groups.append(attributes)
        elif f[0] == "@PG":
            programs.append(attributes)
    versions = sorted({p.get("VN", "unknown") for p in programs
                       if "dragen" in (p.get("PN", "") + p.get("ID", "")).lower()})
    return contigs, {"read_groups": read_groups, "programs": programs, "dragen_versions": versions}


def match_contig(chrom, contigs):
    if chrom in contigs:
        return chrom
    aliases = {chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
    if chrom in {"M", "MT", "chrM", "chrMT"}:
        aliases |= {"M", "MT", "chrM", "chrMT"}
    available = aliases.intersection(contigs)
    if len(available) == 1:
        return available.pop()
    raise InputError(f"Contig {chrom!r} absent or ambiguous in alignment header. Use an exact contig name.")


def prepare_probes(probes, contigs):
    result = []
    for p in probes:
        chrom = match_contig(p.chrom, contigs)
        if p.end > contigs[chrom]:
            raise InputError(f"Probe {p.region} exceeds contig length {contigs[chrom]}.")
        result.append(dataclasses.replace(p, chrom=chrom))
    return result


def query_windows(probes, padding, contigs):
    by_chrom = collections.defaultdict(list)
    for p in probes:
        by_chrom[p.chrom].append((max(0, p.start - padding), min(contigs[p.chrom], p.end + padding)))
    windows = []
    for chrom in sorted(by_chrom):
        merged = []
        for start, end in sorted(by_chrom[chrom]):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        windows.extend((chrom, start, end) for start, end in merged)
    return windows


def target_mask(alignment, probes):
    return sum(1 << i for i, p in enumerate(probes)
               if alignment.chrom == p.chrom and
               any(a < p.end and b > p.start for a, b in alignment.blocks))


def collect_points(lines, probes, args):
    points, stats = {}, collections.Counter()
    for line in lines:
        if line.startswith("@") or not line.strip():
            continue
        stats["records_fetched"] += 1
        try:
            qname, rg, loc, alignment = parse_sam_record(line)
        except (ValueError, IndexError):
            stats["malformed_or_coordinate_free_records"] += 1
            continue
        if args.read_group and rg not in args.read_group:
            stats["read_group_filtered_records"] += 1
            continue
        # Also enforce filters here so the function has identical behaviour on SAM fixtures.
        if alignment.flag & (4 | 256 | 512 | 1024) or alignment.mapq < args.mapq:
            stats["filtered_records"] += 1
            continue
        if args.exclude_supplementary and alignment.flag & 2048:
            stats["supplementary_filtered_records"] += 1
            continue
        alignment.targets = target_mask(alignment, probes)
        key = (qname, rg)
        if key not in points:
            if len(points) >= args.max_fragments:
                raise InputError(f"More than {args.max_fragments:,} fragments fetched. Narrow probes/padding or raise --max-fragments.")
            points[key] = Point(qname, rg, loc)
        point = points[key]
        point.alignments.append(alignment)
        point.targets |= alignment.targets
        if alignment.flag & 2048:
            point.supplementary_targets |= alignment.targets
        else:
            point.primary_targets |= alignment.targets
        stats["records_retained"] += 1
        for tag in alignment.tags:
            if tag in TAGS:
                stats[f"tag_{tag}_records"] += 1
    kept = []
    for point in points.values():
        bx_values = {a.tags["BX"] for a in point.alignments if a.tags.get("BX")}
        point.bx_conflict = len(bx_values) > 1
        point.bx = next(iter(bx_values)) if len(bx_values) == 1 else None
        if point.bx_conflict:
            stats["conflicting_BX_fragments"] += 1
        if args.haplotype is not None and not any(
                a.tags.get("HZ" if args.phase_source == "raw" else "HP") == args.haplotype
                for a in point.alignments):
            stats["haplotype_filtered_fragments"] += 1
            continue
        kept.append(point)
    kept.sort(key=lambda p: (p.location.tile_key, p.rg, p.qname))
    stats["fragments_retained"] = len(kept)
    stats["target_fragments"] = sum(bool(p.targets) for p in kept)
    stats["context_fragments"] = sum(not p.targets for p in kept)
    return kept, dict(stats)


def load_points(args, probes, windows):
    with tempfile.TemporaryDirectory(prefix="virtual-fish-") as tmp:
        bed = Path(tmp) / "windows.bed"
        bed.write_text("".join(f"{c}\t{a}\t{b}\n" for c, a, b in windows))
        mask = 4 | 256 | 512 | 1024 | (2048 if args.exclude_supplementary else 0)
        cmd = [args.samtools, "view", "-M", "-L", str(bed), "-F", str(mask), "-q", str(args.mapq)]
        if args.reference:
            cmd += ["-T", args.reference]
        cmd.append(args.input)
        # File-backed stderr prevents deadlock if samtools emits many diagnostics.
        with tempfile.TemporaryFile(mode="w+t") as errors:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errors, text=True)
            try:
                assert proc.stdout is not None
                points, stats = collect_points(proc.stdout, probes, args)
                proc.stdout.close()
                code = proc.wait()
            except BaseException:
                proc.terminate()
                proc.wait()
                raise
            errors.seek(0)
            diagnostic = errors.read().strip()
            if code:
                raise InputError(f"samtools region query failed: {diagnostic}")
            if diagnostic:
                print(diagnostic, file=sys.stderr)
    return points, stats


def point_namespace(point):
    return (*point.location.tile_key, point.rg)


def near_edges(points, radius, max_edges):
    """Generate radius neighbours among target fragments in each tile/read group."""
    grid = collections.defaultdict(list)
    edges, checks = [], 0
    for i, p in enumerate(points):
        if not p.targets:
            continue
        loc = p.location
        cx, cy = math.floor(loc.x / radius), math.floor(loc.y / radius)
        namespace = point_namespace(p)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in grid.get((*namespace, cx + dx, cy + dy), ()):
                    checks += 1
                    if checks > max_edges * 100:
                        raise InputError("Spatial query is too dense. Reduce --radius or narrow probes.")
                    other = points[j].location
                    distance = math.hypot(loc.x - other.x, loc.y - other.y)
                    if distance <= radius:
                        edges.append((j, i, distance))
                        if len(edges) > max_edges:
                            raise InputError(f"More than {max_edges:,} spatial edges; reduce --radius or raise --max-edges.")
        grid[(*namespace, cx, cy)].append(i)
    return edges


def make_groups(points, edges, mode):
    has_bx = any(p.bx for p in points if p.targets)
    selected = ("bx" if has_bx else "spatial") if mode == "auto" else mode
    buckets = collections.defaultdict(list)
    if selected == "bx":
        for i, point in enumerate(points):
            if point.bx:
                buckets[(*point_namespace(point), point.bx)].append(i)
    elif selected == "spatial":
        parent = list(range(len(points)))

        def root(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for a, b, _ in edges:
            parent[root(b)] = root(a)
        for i, point in enumerate(points):
            if point.targets:
                buckets[(root(i),)].append(i)
    groups = []
    for members in buckets.values():
        targets = 0
        for i in members:
            targets |= points[i].targets
        if not targets:
            continue
        group_id = f"{selected}:{len(groups) + 1}"
        for i in members:
            points[i].group = group_id
        a_only = sum(points[i].targets == 1 for i in members)
        b_only = sum(points[i].targets == 2 for i in members)
        groups.append({"id": group_id, "kind": "inferred BX template" if selected == "bx" else "spatial candidate",
                       "bx": points[members[0]].bx if selected == "bx" else None,
                       "tile_key": list(points[members[0]].location.tile_key),
                       "read_group": points[members[0]].rg, "members": members, "targets": targets,
                       "target_fragments": sum(bool(points[i].targets) for i in members),
                       "A_only_fragments": a_only, "B_only_fragments": b_only,
                       "dual_probe_fragments": sum(points[i].targets == 3 for i in members),
                       "distinct_A_B_support": bool(a_only and b_only)})
    groups.sort(key=lambda g: (-g["distinct_A_B_support"], -g["target_fragments"], g["id"]))
    return selected, groups


def phase_marker(point, source, kind):
    tag = "HZ" if source == "raw" else "HP"
    values = {a.tags[tag] for a in point.alignments if a.tags.get(tag)}
    if not values:
        return "x"
    if (kind == "mrjd" and source == "reported") or len(values) > 1:
        return "D"
    value = next(iter(values))
    return {"1": "o", "2": "s"}.get(value, "D")


def detect_kind(points, path, requested):
    if requested != "auto":
        return requested
    if ".mrjd." in Path(path).name.lower() or any(
            a.tags.get("HP") not in {None, "1", "2"} for p in points for a in p.alignments):
        return "mrjd"
    return "germline"


def summarize(points, probes, groups, edges, stats, args, header_info, kind, selected):
    nearby = [(a, b, d) for a, b, d in edges if {points[a].targets, points[b].targets} == {1, 2}]
    versions = header_info["dragen_versions"]
    if args.phase_source == "raw":
        phase_description = "optional pre-VC HZ/pz; no fallback to HP/pp"
    elif versions and all(v.startswith("4.6") for v in versions):
        phase_description = "reported HP/pp; DRAGEN 4.6 default is post-VC haplotagging"
    elif versions and all(v.startswith("4.5") for v in versions):
        phase_description = "reported HP/pp; DRAGEN 4.5 TruPath uses pre-VC read phasing"
    else:
        phase_description = "reported HP/pp; phase stage not inferred from this header"
    if kind == "mrjd" and args.phase_source == "reported":
        phase_description = "MRJD HP copy assignments (PC confidence), not germline HP 1/2"
    warnings = []
    if stats.get("malformed_or_coordinate_free_records"):
        warnings.append(f"Skipped {stats['malformed_or_coordinate_free_records']} malformed or coordinate-free records.")
    if stats.get("conflicting_BX_fragments"):
        warnings.append("Conflicting BX tags were excluded from template grouping.")
    if selected == "spatial":
        warnings.append("Spatial candidates are radius-connected groups, not validated molecules; chains can merge unrelated reads.")
    if selected == "bx" and not any(p.bx for p in points if p.targets):
        warnings.append("No target-bearing fragments carry a usable BX tag; all target hits remain visible.")
    if args.phase_source == "raw" and not any(a.tags.get("HZ") for p in points for a in p.alignments):
        warnings.append("HZ tags are absent; raw phasing is unavailable. Enable --vc-include-raw-read-phase-scores in DRAGEN to retain them.")
    if args.phase_source == "reported" and not versions:
        warnings.append("No DRAGEN version found in @PG; HP stage is unknown.")
    if kind == "mrjd":
        warnings.append("MRJD alignments may project copies into several paralog regions; probe membership is not an independent locus-origin assignment.")
    if any(p.targets and not p.primary_targets for p in points):
        warnings.append("Some target support comes only from supplementary alignments; see per-alignment evidence in JSON/TSV.")
    if len(probes) == 2 and probes[0].chrom == probes[1].chrom and (
            probes[0].start < probes[1].end and probes[1].start < probes[0].end):
        warnings.append("Probe intervals overlap; dual-probe hits can follow directly from interval overlap.")
    for i, probe in enumerate(probes):
        if not any(p.targets & (1 << i) for p in points):
            warnings.append(f"No retained aligned-base support for probe {probe.label} ({probe.region}).")
    return {"schema_version": VERSION, "title": args.title or Path(args.input).name,
            "input": str(Path(args.input).resolve()),
            "reference": str(Path(args.reference).resolve()) if args.reference else None,
            "annotation": str(Path(args.annotation).resolve()) if args.annotation else None,
            "annotation_build": args.annotation_build, "probes": [dataclasses.asdict(p) for p in probes],
            "coordinate_convention": "0-based half-open genomic intervals; original QNAME X/Y units",
            "settings": {k: getattr(args, k) for k in ("mapq", "padding", "radius", "grouping", "phase_source", "haplotype", "exclude_supplementary")},
            "header": header_info, "input_kind": kind, "phase_description": phase_description,
            "grouping_used": selected, "counts": {**stats,
                "A_fragments": sum(bool(p.targets & 1) for p in points),
                "B_fragments": sum(bool(p.targets & 2) for p in points),
                "dual_probe_fragments": sum(p.targets == 3 for p in points),
                "target_tiles": len({p.location.tile_key for p in points if p.targets}),
                "groups": len(groups), "groups_with_distinct_A_B_support": sum(g["distinct_A_B_support"] for g in groups),
                "nearby_A_only_B_only_pairs": len(nearby)},
            "interpretation": "Exploratory flowcell visualization. Neighbour counts have no calibrated significance; grouping and phase labels are inferred. Query windows do not recover complete templates.",
            "warnings": warnings}, nearby


def write_tables(prefix, points, groups, nearby, summary):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(str(prefix) + ".points.tsv", "w", newline="") as handle:
        columns = ["point_id", "qname", "RG", "instrument", "run", "flowcell", "lane", "tile", "x", "y",
                   "targets", "primary_targets", "supplementary_targets", "BX", "BX_conflict", "group", "phase_records"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for i, p in enumerate(points):
            writer.writerow({"point_id": i, "qname": p.qname, "RG": p.rg, **dataclasses.asdict(p.location),
                             "targets": p.targets, "primary_targets": p.primary_targets,
                             "supplementary_targets": p.supplementary_targets, "BX": p.bx or "",
                             "BX_conflict": p.bx_conflict, "group": p.group or "",
                             "phase_records": json.dumps([{ "chrom": a.chrom, "start": a.start,
                                  "flag": a.flag, **a.tags} for a in p.alignments], separators=(",", ":"))})
    with open(str(prefix) + ".groups.tsv", "w", newline="") as handle:
        columns = ["id", "kind", "bx", "tile_key", "read_group", "target_fragments", "A_only_fragments", "B_only_fragments", "dual_probe_fragments", "distinct_A_B_support", "members"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for g in groups:
            writer.writerow({**g, "tile_key": json.dumps(g["tile_key"]), "members": ",".join(map(str, g["members"]))})
    with open(str(prefix) + ".nearby.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["A_point_id", "B_point_id", "distance_QNAME_units", "same_BX_group"])
        for a, b, distance in nearby:
            if points[a].targets == 2:
                a, b = b, a
            same = bool(points[a].group and points[a].group == points[b].group and points[a].group.startswith("bx:"))
            writer.writerow([a, b, f"{distance:.6f}", same])
    payload = {**summary, "groups": groups, "points": [dataclasses.asdict(p) for p in points]}
    Path(str(prefix) + ".json").write_text(json.dumps(payload, indent=2) + "\n")


def scatter_points(ax, points, source, kind, large=False):
    buckets = collections.defaultdict(list)
    for p in points:
        buckets[(p.targets, phase_marker(p, source, kind))].append(p)
    for (targets, marker), members in sorted(buckets.items()):
        size = (65 if large else 18) if targets else (14 if large else 5)
        ax.scatter([p.location.x for p in members], [p.location.y for p in members],
                   c=COLORS[targets], marker=marker, s=size, alpha=.9 if targets else .25,
                   linewidths=.7, zorder=3 if targets else 1)
    for p in points:
        if p.targets and p.supplementary_targets & ~p.primary_targets:
            ax.scatter([p.location.x], [p.location.y], s=95 if large else 34,
                       facecolors="none", edgecolors="#222", linewidths=.75, zorder=4)


def fit_xy(ax, points):
    if not points:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    else:
        xs, ys = [p.location.x for p in points], [p.location.y for p in points]
        extent = max(max(xs) - min(xs), max(ys) - min(ys), 100)
        half = extent * .58
        cx, cy = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
        ax.set_xlim(cx - half, cx + half)
        ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal", adjustable="box")
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.set_xlabel("X · read-name units", fontsize=8)
    ax.set_ylabel("Y · read-name units", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=.13)


def plot_outputs(prefix, points, probes, groups, summary, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    kind = summary["input_kind"]
    labels = {0: "queried context", 1: probes[0].label}
    if len(probes) == 2:
        labels.update({2: probes[1].label, 3: "both probes · one read pair"})
    legend = [Line2D([], [], color=COLORS[k], marker="o", ls="", label=v) for k, v in labels.items()]
    phase_tag = "HZ" if args.phase_source == "raw" else "HP"
    if kind == "germline" or args.phase_source == "raw":
        legend += [Line2D([], [], color="#555", marker=m, ls="", label=l)
                   for m, l in [("o", f"{phase_tag} 1"), ("s", f"{phase_tag} 2"), ("x", "unassigned"), ("D", "mixed labels")]]
    else:
        legend += [Line2D([], [], color="#555", marker="D", ls="", label="MRJD copy label"),
                   Line2D([], [], color="#555", marker="x", ls="", label="unassigned")]
    tiles = sorted({p.location.tile_key for p in points if p.targets})
    files = []
    for page in range(max(1, math.ceil(len(tiles) / args.tiles_per_page))):
        chosen = tiles[page * args.tiles_per_page:(page + 1) * args.tiles_per_page]
        columns = min(3, max(1, len(chosen)))
        rows = max(1, math.ceil(len(chosen) / columns))
        fig, axes = plt.subplots(rows, columns, figsize=(5 * columns, 4.5 * rows + 1.2), squeeze=False)
        for ax, tile in zip(axes.flat, chosen):
            members = [p for p in points if p.location.tile_key == tile]
            scatter_points(ax, members, args.phase_source, kind)
            fit_xy(ax, members)
            inst, run, fc, lane, number = tile
            counts = collections.Counter(p.targets for p in members)
            ax.set_title(f"{fc} · run {run} · lane {lane} · tile {number}\n"
                         f"A {counts[1] + counts[3]}  B {counts[2] + counts[3]}  context {counts[0]}", fontsize=9)
        for ax in list(axes.flat)[len(chosen):]:
            ax.set_axis_off()
        if not chosen:
            axes[0, 0].text(.5, .5, "No retained probe hits", ha="center", va="center", transform=axes[0, 0].transAxes)
        fig.suptitle(f"TruPath virtual probes · {summary['title']}", fontsize=15, y=.99)
        fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(.5, .95), ncol=min(4, len(legend)), fontsize=8)
        fig.text(.5, .018, "Tile panels show queried extents; each tile has its own scale. DNA locations after extraction.\n" + summary["phase_description"], ha="center", fontsize=8, color="#555")
        fig.tight_layout(rect=(0, .12, 1, .85 if rows == 1 else .91))
        suffix = ".atlas" + (f".{page + 1:03d}" if len(tiles) > args.tiles_per_page else "")
        for fmt in ("png", "svg"):
            out = str(prefix) + suffix + "." + fmt
            fig.savefig(out, dpi=args.dpi, facecolor="white")
            files.append(out)
        plt.close(fig)
    chosen_groups = [g for g in groups if g["target_fragments"] >= args.min_group_reads][:args.top]
    if chosen_groups:
        n = len(chosen_groups)
        fig = plt.figure(figsize=(15, 3.2 * n + 1.0))
        grid = fig.add_gridspec(n, 1 + len(probes), width_ratios=[1.25] + [1] * len(probes))
        for row, g in enumerate(chosen_groups):
            members = [points[i] for i in g["members"]]
            ax = fig.add_subplot(grid[row, 0])
            scatter_points(ax, members, args.phase_source, kind, large=True)
            fit_xy(ax, members)
            bx_label = (g["bx"][:42] + "…") if g["bx"] and len(g["bx"]) > 42 else g["bx"]
            _, run, flowcell, lane, tile = g["tile_key"]
            ax.set_title(f"{g['id']} · {g['kind']}\n"
                         f"{bx_label or 'radius connected'} · {g['target_fragments']} target pairs\n"
                         f"{flowcell} · run {run} · lane {lane} · tile {tile}", fontsize=8)
            for col, probe in enumerate(probes, 1):
                b = fig.add_subplot(grid[row, col])
                b.axvspan(probe.start + 1, probe.end, color=COLORS[1 << (col - 1)], alpha=.12)
                offsets = []
                for offset, point in enumerate(members):
                    for alignment in point.alignments:
                        if alignment.chrom != probe.chrom:
                            continue
                        for start, end in alignment.blocks:
                            if start >= probe.end + args.padding or end <= probe.start - args.padding:
                                continue
                            b.plot([start + 1, end], [offset, offset], color=COLORS[point.targets],
                                   lw=2, ls="--" if alignment.flag & 2048 else "-", alpha=.9)
                            offsets.append(offset)
                lo, hi = max(1, probe.start + 1 - args.padding), probe.end + args.padding
                b.set_xlim(lo - .5, hi + .5)
                b.set_ylim(-1, max(offsets, default=0) + 1)
                b.set_title(f"{probe.label}\n{probe.region}", fontsize=9)
                b.set_xlabel(f"{probe.chrom} · genomic position (bp)", fontsize=8)
                b.set_ylabel("read-pair row", fontsize=8)
                b.ticklabel_format(style="plain", useOffset=False)
                b.tick_params(labelsize=7)
                b.grid(alpha=.13)
        fig.suptitle(f"Target-bearing groups · {summary['title']}", fontsize=14, y=.997)
        fig.text(.5, .008, "Colours mark probe membership; dashed genomic segments are supplementary. Groups and haplotypes are inferred.", ha="center", fontsize=8, color="#555")
        fig.tight_layout(rect=(0, .035, 1, .97))
        for fmt in ("png", "svg"):
            out = str(prefix) + ".groups." + fmt
            fig.savefig(out, dpi=args.dpi, facecolor="white")
            files.append(out)
        plt.close(fig)
    return files


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "--cram", "--bam", dest="input", required=True, help="Indexed, coordinate-sorted BAM/CRAM.")
    p.add_argument("--reference", "-T", help="Matching reference FASTA; required for CRAM.")
    p.add_argument("--region", action="append", default=[], help="One or two [LABEL=]chr:start-end intervals; repeat this option.")
    p.add_argument("--gene", action="append", default=[], help="Exact gene symbol/ID; repeat for two genes. SYMBOL@contig disambiguates.")
    p.add_argument("--genes", nargs="+", default=[], help="Alternative to repeated --gene.")
    p.add_argument("--annotation", help="Reference-matched GTF/GFF3, optionally gzip-compressed.")
    p.add_argument("--annotation-build", help="Annotation build label recorded in provenance, e.g. GRCh38. No liftover is performed.")
    p.add_argument("--out-prefix", default="virtual_fish", help="Prefix for PNG/SVG/TSV/JSON outputs.")
    p.add_argument("--title", help="Plot title; defaults to the input filename.")
    p.add_argument("--mapq", type=int, default=20, help="Minimum ordinary SAM MAPQ (not the xq tag); default 20.")
    p.add_argument("--padding", type=int, default=10000, help="Genomic context around probes; does not expand probe membership. Default 10000 bp.")
    p.add_argument("--grouping", choices=("auto", "bx", "spatial", "none"), default="auto", help="Auto uses BX if present on targets, otherwise exploratory spatial groups.")
    p.add_argument("--radius", type=float, default=350, help="Spatial neighbour radius in QNAME coordinate units; default 350. No genomic distance cutoff.")
    p.add_argument("--phase-source", choices=("reported", "raw"), default="reported", help="Reported HP/pp, or optional pre-VC HZ/pz. Raw never substitutes HP.")
    p.add_argument("--haplotype", help="Keep a fragment if any fetched alignment has this exact HP (or HZ in raw mode), including MRJD copy strings.")
    p.add_argument("--input-kind", choices=("auto", "germline", "mrjd"), default="auto", help="MRJD HP is a copy label, not the germline numeric haplotype.")
    p.add_argument("--read-group", action="append", default=[], help="Restrict to exact RG IDs; repeat as needed.")
    p.add_argument("--exclude-supplementary", action="store_true", help="Exclude supplementary evidence; by default it is retained without counting another point.")
    p.add_argument("--max-fragments", type=int, default=500000, help="Memory guard: fail rather than silently downsample. Default 500000.")
    p.add_argument("--max-edges", type=int, default=1000000, help="Density guard for spatial radius edges; default 1000000.")
    p.add_argument("--top", type=int, default=8, help="Number of group zooms; every probe hit remains in the atlas. Default 8.")
    p.add_argument("--min-group-reads", type=int, default=2, help="Minimum target read pairs for a zoom; default 2.")
    p.add_argument("--tiles-per-page", type=int, default=12, help="Atlas pagination; default 12 tiles per PNG/SVG page.")
    p.add_argument("--dpi", type=int, default=160)
    p.add_argument("--tables-only", action="store_true", help="Write TSV/JSON without plotting or requiring matplotlib.")
    p.add_argument("--samtools", default="samtools", help="Path/name of samtools executable.")
    p.add_argument("--version", action="version", version=VERSION)
    args = p.parse_args(argv)
    if not 1 <= len(args.region) + len(args.gene) + len(args.genes) <= 2:
        p.error("Specify one or two probes using --region, --gene or --genes.")
    if (args.gene or args.genes) and not args.annotation:
        p.error("--gene/--genes requires --annotation.")
    if Path(args.input).suffix.lower() == ".cram" and not args.reference:
        p.error("CRAM input requires --reference.")
    if not 0 <= args.mapq <= 255 or args.padding < 0:
        p.error("--mapq must be 0..255 and --padding must be non-negative.")
    if not math.isfinite(args.radius) or args.radius <= 0:
        p.error("--radius must be finite and positive.")
    if any(getattr(args, k) < 1 for k in ("max_fragments", "max_edges", "top", "min_group_reads", "tiles_per_page", "dpi")):
        p.error("Limits, --top, --min-group-reads, --tiles-per-page and --dpi must be positive.")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        if not shutil.which(args.samtools):
            raise InputError(f"samtools executable not found: {args.samtools}")
        for name in ("input", "reference", "annotation"):
            value = getattr(args, name)
            if value and not Path(value).is_file():
                raise InputError(f"{name} file not found: {value}")
        if not args.tables_only:
            try:
                import matplotlib  # noqa: F401
            except ImportError:
                raise InputError("Install matplotlib, or use --tables-only.") from None
        probes = [parse_region(value) for value in args.region]
        probes += resolve_genes(args.annotation, args.gene + args.genes) if args.gene or args.genes else []
        if len({p.label for p in probes}) != len(probes):
            raise InputError("Probe labels must be distinct; use explicit LABEL=region names.")
        header_cmd = [args.samtools, "view", "-H"]
        if args.reference:
            header_cmd += ["-T", args.reference]
        header = run_capture(header_cmd + [args.input])
        contigs, header_info = parse_header(header)
        probes = prepare_probes(probes, contigs)
        windows = query_windows(probes, args.padding, contigs)
        print("Probes: " + "; ".join(f"{p.label}={p.region}" for p in probes))
        points, stats = load_points(args, probes, windows)
        if not stats.get("records_retained") and stats.get("malformed_or_coordinate_free_records"):
            raise InputError("No usable coordinate-bearing records. Check that original Illumina QNAMEs were retained in this file.")
        # Explicit haplotype filters must not silently produce an empty result when tags are unavailable.
        if args.haplotype is not None and not points:
            raise InputError("No fragments match --haplotype in the fetched windows; check phase source and tags.")
        kind = detect_kind(points, args.input, args.input_kind)
        edges = near_edges(points, args.radius, args.max_edges)
        selected, groups = make_groups(points, edges, args.grouping)
        summary, nearby = summarize(points, probes, groups, edges, stats, args, header_info, kind, selected)
        summary["query_windows"] = [{"chrom": c, "start": a, "end": b} for c, a, b in windows]
        prefix = Path(args.out_prefix)
        write_tables(prefix, points, groups, nearby, summary)
        plots = [] if args.tables_only else plot_outputs(prefix, points, probes, groups, summary, args)
        print(f"Retained {stats['target_fragments']:,} target fragments across {summary['counts']['target_tiles']} tiles; grouping={selected}.")
        print(summary["phase_description"])
        for warning in summary["warnings"]:
            print("Note: " + warning, file=sys.stderr)
        print(f"Wrote {prefix}.json, .points.tsv, .groups.tsv, .nearby.tsv" + (f" and {len(plots)} plots." if plots else "."))
        return 0
    except (InputError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
