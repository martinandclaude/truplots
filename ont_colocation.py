#!/usr/bin/env python3
"""
ont_colocation.py - build a TruPath-style colocation map from long-read alignments (ONT, PacBio).

A TruPath colocation map counts, for every pair of genomic bins, how often reads from the two bins
came from the same long DNA molecule. A long read is itself one molecule, so the same map can be
built from long-read alignments: each read adds one to every pair of bins it covers, including the
bins of its split (supplementary) alignments listed in the SA tag. Deletions, inversions,
translocations and amplicon or ecDNA junctions then show up as the same off-diagonal signals, and

  plot_colocation.py TRUPATH.cooler --compare ONT.cool

draws both technologies in one figure: TruPath above the diagonal, long reads below.

--region counts only reads with an alignment in those regions (e.g. an amplicon). All of a read's
segments are still counted wherever they map, so co-amplified segments elsewhere appear too; give
them their own --region to also count reads that lie entirely inside them. Without --region every
read in the file is counted; for a whole genome choose a coarser --binsize, since 2 kb bins
genome-wide need several GB of memory.

The output follows the cooler v3 schema (written with h5py), so cooler zoomify / balance, HiGlass and
plot_colocation.py all read it. --like copies the bin size of the TruPath cooler so the two maps share
one grid.

Examples:
  ont_colocation.py --input ont.bam --like sample.colocation.cooler \\
      --region chr2:13,000,000-18,000,000 -o ont.mycn.cool
  plot_colocation.py sample.colocation.cooler --compare ont.mycn.cool --names TruPath,ONT \\
      --region chr2:15,500,000-16,300,000

Requires: Python 3 with pysam, numpy, h5py.
"""
import argparse
import collections
import datetime
import functools
import json
import re
import sys
from pathlib import Path

import numpy as np

REF_OPS = {0, 2, 3, 7, 8}              # M D N = X consume the reference
GAP_OPS = {2, 3}                       # D N: a long one leaves bins uncovered
CIGAR_CODES = {op: i for i, op in enumerate("MIDNSHP=X")}
COMPACT = 20_000_000                   # pending bin pairs before they are merged


class InputError(Exception):
    pass


def parse_cigar(text):
    return [(CIGAR_CODES[op], int(n)) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", text)]


def blocks(start, cigar, min_gap):
    """Reference intervals covered by one alignment, split only at deletions or skips of at
    least min_gap bp (a shorter gap cannot leave a whole bin uncovered)."""
    out, first, pos = [], start, start
    for op, n in cigar:
        if op in GAP_OPS and n >= min_gap:
            if pos > first:
                out.append((first, pos))
            pos += n
            first = pos
        elif op in REF_OPS:
            pos += n
    if pos > first:
        out.append((first, pos))
    return out


def segments(read, mapq, min_gap):
    """(contig, start, end) blocks of every alignment of the read: this record plus its SA tag."""
    out = []
    if read.mapping_quality >= mapq:
        out += [(read.reference_name, s, e) for s, e in blocks(read.reference_start, read.cigartuples, min_gap)]
    if read.has_tag("SA"):
        for entry in read.get_tag("SA").rstrip(";").split(";"):
            rname, pos, _strand, cigar, q, _nm = entry.split(",")
            if int(q) >= mapq:
                out += [(rname, s, e) for s, e in blocks(int(pos) - 1, parse_cigar(cigar), min_gap)]
    return out


@functools.lru_cache(maxsize=4096)
def upper_pairs(k):
    return np.triu_indices(k)


class PairCounts:
    """Counts of bin pairs (bin1 <= bin2) as sorted int64 keys bin1 * nbins + bin2."""

    def __init__(self, nbins):
        self.nbins = nbins
        self.keys = np.zeros(0, dtype=np.int64)
        self.counts = np.zeros(0, dtype=np.int64)
        self.pending, self.npending = [], 0

    def add(self, bins):
        i, j = upper_pairs(len(bins))
        self.pending.append(bins[i] * self.nbins + bins[j])
        self.npending += len(i)
        if self.npending >= COMPACT:
            self.compact()

    def compact(self):
        if not self.pending:
            return
        keys = np.concatenate([self.keys] + self.pending)
        weights = np.concatenate([self.counts, np.ones(self.npending, dtype=np.int64)])
        self.keys, inverse = np.unique(keys, return_inverse=True)
        self.counts = np.bincount(inverse, weights=weights).astype(np.int64)
        self.pending, self.npending = [], 0


def match_contig(name, contigs):
    if name in contigs:
        return name
    alias = name[3:] if name.startswith("chr") else "chr" + name
    if alias in contigs:
        return alias
    raise InputError(f"Contig {name!r} not found in the alignment header.")


def parse_region(value, contigs):
    m = re.fullmatch(r"(.+):([\d,]+)-([\d,]+)", value.strip())
    chrom = match_contig(m.group(1) if m else value.strip(), contigs)
    if not m:
        return chrom, 0, contigs[chrom]
    start, end = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if not 1 <= start <= end <= contigs[chrom]:
        raise InputError(f"Invalid region {value!r}: need 1 <= start <= end <= {contigs[chrom]:,}.")
    return chrom, start - 1, end


def like_cooler(uri, contigs):
    """Bin size of the reference (TruPath) cooler, and contigs whose lengths differ from the BAM's."""
    import h5py
    from plot_colocation import Cooler, find_resolutions, split_uri
    path, group = split_uri(uri)
    if not Path(path).is_file():
        raise InputError(f"File not found: {path}")
    with h5py.File(path, "r") as h5:
        res = find_resolutions(h5, group)
        clr = Cooler(h5, res[min(res)])
        lengths = dict(zip(clr.chroms, clr.lengths.tolist()))
    differ = []
    for name, length in contigs.items():
        try:
            other = lengths[match_contig(name, lengths)]
        except InputError:
            continue
        if other != length:
            differ.append(f"{name} {length:,} vs {other:,}")
    return clr.binsize, differ


def count_reads(args, contigs, offsets, nbins):
    """Stream the alignments and add each read's covered bins to the pair counts."""
    import pysam
    mode = "rc" if Path(args.input).suffix.lower() == ".cram" else "rb"
    bam = pysam.AlignmentFile(args.input, mode, reference_filename=args.reference, threads=args.threads)
    pairs, stats, seen = PairCounts(nbins), collections.Counter(), set()
    regions = [parse_region(r, contigs) for r in args.region]
    records = (read for chrom, start, end in regions for read in bam.fetch(chrom, start, end)) \
        if regions else bam.fetch(until_eof=True)
    for read in records:
        if read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate:
            continue
        if regions:
            # A read can show up in several regions, or through a supplementary record; count it once.
            if read.query_name in seen:
                continue
            seen.add(read.query_name)
        elif read.is_supplementary:
            continue                   # counted through its primary record's SA tag
        if args.min_length and read.infer_read_length() < args.min_length:
            stats["shorter than --min-length"] += 1
            continue
        segs = [(offsets[c], s, e) for c, s, e in segments(read, args.mapq, args.binsize) if c in offsets]
        if not segs:
            stats["no alignment with MAPQ >= --mapq"] += 1
            continue
        bins = np.unique(np.concatenate([np.arange(off + s // args.binsize, off + (e - 1) // args.binsize + 1)
                                         for off, s, e in segs]))
        pairs.add(bins)
        stats["reads counted"] += 1
        stats["split reads (SA)"] += read.has_tag("SA")
    bam.close()
    pairs.compact()
    return pairs, stats


def write_cooler(path, contigs, binsize, pairs, metadata, assembly):
    """Write a cooler v3 file (symmetric-upper storage) with h5py."""
    import h5py
    names = list(contigs)
    lengths = np.array([contigs[n] for n in names], dtype=np.int64)
    per_chrom = -(-lengths // binsize)
    chrom_offset = np.concatenate([[0], np.cumsum(per_chrom)]).astype(np.int64)
    nbins = int(chrom_offset[-1])
    starts = np.concatenate([np.arange(n, dtype=np.int64) * binsize for n in per_chrom])
    ends = np.minimum(starts + binsize, np.repeat(lengths, per_chrom))
    bin1, bin2 = pairs.keys // nbins, pairs.keys % nbins
    bin1_offset = np.searchsorted(bin1, np.arange(nbins + 1), side="left").astype(np.int64)
    gz = dict(compression="gzip", compression_opts=6, shuffle=True)
    with h5py.File(path, "w") as f:
        f.attrs.update({"format": "HDF5::Cooler", "format-version": 3, "format-url": "https://github.com/open2c/cooler",
                        "bin-type": "fixed", "bin-size": binsize, "storage-mode": "symmetric-upper",
                        "nchroms": len(names), "nbins": nbins, "nnz": len(bin1), "sum": int(pairs.counts.sum()),
                        "generated-by": "truplots ont_colocation.py",
                        "creation-date": datetime.datetime.now().isoformat(timespec="seconds"),
                        "metadata": json.dumps(metadata)})
        if assembly:
            f.attrs["genome-assembly"] = assembly
        width = max(len(n) for n in names)
        chroms = f.create_group("chroms")
        chroms.create_dataset("name", data=np.array(names, dtype=f"S{width}"), **gz)
        chroms.create_dataset("length", data=lengths.astype(np.int32), **gz)
        bins = f.create_group("bins")
        enum = h5py.enum_dtype({n: i for i, n in enumerate(names)}, basetype="<i4")
        bins.create_dataset("chrom", data=np.repeat(np.arange(len(names), dtype=np.int32), per_chrom),
                            dtype=enum, **gz)
        bins.create_dataset("start", data=starts.astype(np.int32), **gz)
        bins.create_dataset("end", data=ends.astype(np.int32), **gz)
        pixels = f.create_group("pixels")
        for name, data, dtype in (("bin1_id", bin1, np.int64), ("bin2_id", bin2, np.int64),
                                  ("count", pairs.counts, np.int32)):
            pixels.create_dataset(name, data=data.astype(dtype), maxshape=(None,),
                                  chunks=(min(max(len(data), 1), 1 << 20),), **gz)
        indexes = f.create_group("indexes")
        indexes.create_dataset("chrom_offset", data=chrom_offset, **gz)
        indexes.create_dataset("bin1_offset", data=bin1_offset, **gz)
    return nbins


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "--bam", dest="input", required=True,
                   help="Indexed, coordinate-sorted long-read BAM/CRAM (e.g. minimap2), with SA tags on split reads.")
    p.add_argument("--reference", "-T", help="Reference FASTA; required for CRAM.")
    p.add_argument("-o", "--out", help="Output cooler (default: <input>.colocation.cool).")
    p.add_argument("--binsize", type=int, help="Bin size in bp (default 2000, or that of --like).")
    p.add_argument("--like", help="TruPath colocation cooler (or .mcool) whose bin size to copy; also warns "
                                  "when contig lengths differ (another reference build).")
    p.add_argument("-r", "--region", action="append", default=[],
                   help="chrom or chrom:start-end (1-based, inclusive); count reads with an alignment here. "
                        "Repeatable. Default: every read in the file.")
    p.add_argument("--mapq", type=int, default=20, help="Minimum MAPQ per alignment segment (default 20).")
    p.add_argument("--min-length", type=int, default=0, help="Skip reads shorter than this, in bp (default 0).")
    p.add_argument("--assembly", help="Genome assembly label stored in the file, e.g. hg38.")
    p.add_argument("--threads", type=int, default=2, help="Decompression threads (default 2).")
    args = p.parse_args(argv)
    if Path(args.input).suffix.lower() == ".cram" and not args.reference:
        p.error("CRAM input requires --reference.")
    if args.binsize is not None and args.binsize < 1 or args.mapq < 0 or args.min_length < 0 or args.threads < 1:
        p.error("--binsize and --threads must be positive; --mapq and --min-length zero or more.")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        import pysam
    except ImportError:
        sys.exit("error: pysam is required (conda install pysam, or pip install pysam).")
    try:
        if not Path(args.input).is_file():
            raise InputError(f"File not found: {args.input}")
        with pysam.AlignmentFile(args.input, "rc" if Path(args.input).suffix.lower() == ".cram" else "rb",
                                 reference_filename=args.reference) as bam:
            contigs = dict(zip(bam.references, bam.lengths))
        if not contigs:
            raise InputError("No contigs in the alignment header; is the file aligned?")
        if args.like:
            like_binsize, differ = like_cooler(args.like, contigs)
            if args.binsize and args.binsize != like_binsize:
                raise InputError(f"--binsize {args.binsize:,} differs from the {like_binsize:,} bp bins of --like.")
            args.binsize = like_binsize
            if differ:
                print(f"warning: contig lengths differ from {args.like} (another reference build?): "
                      + ", ".join(differ[:5]) + (" ..." if len(differ) > 5 else ""), file=sys.stderr)
        args.binsize = args.binsize or 2000
        per_chrom = [-(-n // args.binsize) for n in contigs.values()]
        offsets = dict(zip(contigs, np.concatenate([[0], np.cumsum(per_chrom)[:-1]]).tolist()))
        nbins = int(sum(per_chrom))
        pairs, stats = count_reads(args, contigs, offsets, nbins)
        if not stats["reads counted"]:
            raise InputError("No reads counted; check --region, --mapq and that the file is aligned.")
        out = args.out or Path(args.input).name.rsplit(".", 1)[0] + ".colocation.cool"
        metadata = {"input": str(Path(args.input).resolve()), "regions": args.region, "bin-size": args.binsize,
                    "min-mapq": args.mapq, "min-length": args.min_length,
                    "counting": "each read adds 1 to every pair of bins covered by its alignments (incl. SA)",
                    **{k.replace(" ", "-"): v for k, v in stats.items()}}
        write_cooler(out, contigs, args.binsize, pairs, metadata, args.assembly)
        print(f"{stats['reads counted']:,} reads counted ({stats['split reads (SA)']:,} split across "
              f"several alignments) into {len(pairs.keys):,} bin pairs at {args.binsize:,} bp")
        for reason in ("shorter than --min-length", "no alignment with MAPQ >= --mapq"):
            if stats[reason]:
                print(f"  skipped {stats[reason]:,} reads: {reason}")
        print(f"wrote {out}")
        return 0
    except (InputError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
