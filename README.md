# trupath-qc

Visualisation tooling for Illumina **TruPath** (proximity-mapped read) sequencing QC.

TruPath builds the library on a patterned flow cell so that one long input DNA molecule seeds
sequencing clusters across a small patch of **neighbouring nanowells** — a "constellation". DRAGEN
exploits the fact that *clusters close together on the flow cell come from positions close together in the
genome* to add long-range information to short reads. This repo reconstructs and plots individual
molecules straight from the aligned reads, using the physical flow-cell coordinates encoded in every
Illumina read name (`INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y`).

> Works on any Illumina CRAM/BAM whose read names carry `LANE:TILE:X:Y`, but is only *meaningful* on
> proximity (TruPath) data — standard WGS has no constellations.

## Setup

```bash
conda env create -f environment.yml
conda activate trupath-qc
```


## `plot_molecule.py`

Finds and plots a single molecule (or the top *N*) in a genomic region as a two-panel figure:

- **left** — the molecule's nanowells on the flow cell, coloured by position within the molecule;
- **right** — the same reads laid out on the genome (the reconstructed molecule), annotated with span and
  the genomic gap between read-pairs.

A "molecule" is recovered by clustering reads within each flow-cell tile (DBSCAN on `x`, `y`, and a scaled
genomic coordinate), so a cluster is tight *both* on the chip and in the genome.

### Usage

```bash
python3 plot_molecule.py \
  --cram   sample.cram \
  --reference genome.fa \
  --region chr2:178525989-178830800 \
  --out    molecule.png
```

### Options

| option | meaning |
|--------|---------|
| `--by span` | (default) plot the **longest** molecule |
| `--by reads` | plot the **best-supported** molecule (most read-pairs) |
| `--by density` | plot the **sparsest** molecule (fewest read-pairs per 100 kb = biggest gaps) |
| `--top N` | plot the top *N* molecules as a stacked gallery instead of one |
| `--min-reads` | minimum read-pairs per molecule (default 5) |
| `--max-span-kb` | reject over-long chained clusters above this span |
| `--eps` / `--genomic-eps` | clustering tightness: pixel radius, and the bp distance mapped to it (default 350 px / 60 kb) |
| `--mapq` | minimum mapping quality (default 1) |

The tool prints per-molecule stats (span, read-pairs, median/max read-pair gap, tile) to stdout and
annotates them on the figure.

### Requirements

- [`samtools`](https://www.htslib.org/) on `PATH`
- Python 3 with `numpy`, `matplotlib`, `scikit-learn`

## `plot_vaf.py`

Plots the **variant allele frequency (VAF) distribution** of a VCF/BCF as a histogram — green for
heterozygous calls, blue for homozygous. A clean germline sample is bimodal: a HET peak near 0.5 and a HOM
peak near 1.0. Handy as an allele-balance / sample-QC check: a low-VAF shoulder or a HET peak shifted off
0.5 points to contamination or a sample mixture, whereas simply *more* HET calls still centred on 0.5
indicate a genuinely high-heterozygosity (e.g. admixed) genome.

VAF is taken from `FORMAT/AD` (combined alt depth / total depth — caller-agnostic and multiallelic-aware),
falling back to `FORMAT/AF` or `VAF`; zygosity comes from `FORMAT/GT`.

### Usage

```bash
python3 plot_vaf.py sample.vcf.gz -o vaf.png
python3 plot_vaf.py sample.vcf.gz -s TUMOR --min-dp 20 --pass-only --bins 80
```

### Options

| option | meaning |
|--------|---------|
| `-s, --sample` | sample/column to plot (default: first in the VCF) |
| `--vaf-source` | `ad` (default, from `FORMAT/AD`) or `af` (from `FORMAT/AF`/`VAF`) |
| `--min-dp` | drop variants below this depth (e.g. 20) |
| `--pass-only` | keep only PASS / unfiltered variants |
| `--bins` | histogram bins over 0–1 (default 50) |
| `-o/--output`, `--title`, `--dpi` | output path / title / resolution |

### Requirements

- Python 3 with [`pysam`](https://pysam.readthedocs.io/), `numpy`, `matplotlib`
