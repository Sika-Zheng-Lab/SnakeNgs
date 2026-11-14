# preprocessing_RNAseq.smk

Snakemake workflow for preprocessing **paired-end** and **single-end** bulk RNA-seq data. The `layout` parameter in the config file controls the mode.

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/Sika-Zheng-Lab/SnakeNgs) repository.

## Workflow

### Paired-end (`layout: "paired"`)

1. FASTQ trimming/filtering and read-level QC using [fastp](https://github.com/OpenGene/fastp) with the default parameters.
2. Alignment using [STAR](https://github.com/alexdobin/STAR) with the parameter `--outFilterMultimapNmax 1`.
3. Convert the SAM file to BAM file and sort using [samtools](http://www.htslib.org/).
4. Mark duplicates in a QC-only BAM with `samtools` (retain duplicates), then run [RustQC](https://github.com/seqeralabs/RustQC) `rna` for alignment statistics, gene-body coverage, feature distribution, duplication, splice-junction and library-complexity metrics. Paired-end runs also include inner-distance and insert-size statistics.
5. Make bigWig files using [deepTools](https://deeptools.readthedocs.io/en/develop/) `bamCoverage` with the parameter `--binSize 1`.
6. Make summary statistics using [MultiQC](https://multiqc.info/).

### Single-end (`layout: "single"`)

1. FASTQ trimming/filtering and read-level QC using [fastp](https://github.com/OpenGene/fastp) with the default parameters.
2. Alignment using [STAR](https://github.com/alexdobin/STAR) with the parameter `--outFilterMultimapNmax 1`.
3. Convert the SAM file to BAM file and sort using [samtools](http://www.htslib.org/).
4. Mark duplicates in a QC-only BAM with `samtools` (retain duplicates), then run [RustQC](https://github.com/seqeralabs/RustQC) `rna` for alignment statistics, gene-body coverage, feature distribution, duplication, splice-junction and library-complexity metrics.
5. Make bigWig files using [deepTools](https://deeptools.readthedocs.io/en/develop/) `bamCoverage`.
6. Make summary statistics using [MultiQC](https://multiqc.info/).

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/preprocessing_RNAseq.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

### Paired-end

``` yaml
workdir: path/to/output
samples: ["SRRXXXXXX", "SRRYYYYYY", "SRRZZZZZZ"]
star_index: path/to/star_index
gtf: path/to/reference_transcriptome.gtf
layout: "paired"
```

- `path/to/output` should contain `fastq` directory with the following structure:

``` bash
output/
└── fastq
    ├── SRRXXXXXX_1.fastq.gz
    ├── SRRXXXXXX_2.fastq.gz
    ├── SRRYYYYYY_1.fastq.gz
    ├── SRRYYYYYY_2.fastq.gz
    ├── SRRZZZZZZ_1.fastq.gz
    └── SRRZZZZZZ_2.fastq.gz
```

### Single-end

``` yaml
workdir: path/to/output
samples: ["SRRXXXXXX", "SRRYYYYYY", "SRRZZZZZZ"]
star_index: path/to/star_index
gtf: path/to/reference_transcriptome.gtf
layout: "single"
```

- `path/to/output` should contain `fastq` directory with the following structure:

``` bash
output/
└── fastq
    ├── SRRXXXXXX.fastq.gz
    ├── SRRYYYYYY.fastq.gz
    └── SRRZZZZZZ.fastq.gz
```

### Common settings

- `/path/to/reference_transcriptome.gtf` is the reference transcriptome in GTF format (e.g. `Homo_sapiens.GRCh38.106.gtf` for human transcriptome).

Please refer to the [tutorial](../tutorial/RNAseq_preprocessing.md) for more information.

## RustQC configuration and outputs

RustQC is pinned to **v0.2.1**. Existing configs remain valid. An optional
`rustqc_stranded: unstranded` setting accepts `unstranded` (default), `forward`,
or `reverse`; choose the value appropriate for the library preparation.
`layout` determines whether RustQC receives `--paired`.

- `rustqc/<sample>/`: RustQC's native directory tree, including `samtools/`,
  `rseqc/`, `featurecounts/`, `dupradar/`, `qualimap/` and, when estimable, `preseq/`.
- `rustqc/<sample>/rustqc_summary.json`: machine-readable run summary.
- `log/rustqc_<sample>.log`: diagnostics, including metrics that could not be
  estimated (for example, library complexity with too few reads).
- `multiqc/multiqc_report.html`: fastp, STAR and RustQC outputs for the configured
  samples. RustQC outputs are parsed by MultiQC's existing tool modules.

MultiQC reads the gene-level featureCounts summary and displays biotypes in a
separate custom plot, avoiding same-sample overwrites. Zero-event junction
annotation reports are excluded from MultiQC's percentage plot to avoid its
zero-division error; the original RustQC reports remain on disk.

The original `star/<sample>/<sample>_Aligned.out.bam` remains the input for
bigWig generation and downstream analyses. `rustqc_bam/<sample>.markdup.bam`
is a temporary QC input; reads are marked, not removed. Sorting and marking
require additional temporary disk space under `rustqc_bam/`.

This replaces the RNA workflow's Picard `CollectRnaSeqMetrics` and
`CollectInsertSizeMetrics` reports. RustQC produces different metrics and file
formats: it does **not** reproduce Picard's exact metric definitions or PDF
histogram. Insert sizes are available in `samtools/<sample>.stats` (`IS` rows),
and RNA coverage/read-origin metrics in the Qualimap/RSeQC outputs. The old
`metrics/*.picard.analysis.*`, `refFlat.txt` and ribosomal-interval files are no
longer generated or included in MultiQC; existing files can remain on disk.

RustQC v0.2.1 only provides RNA alignment QC. **fastp** remains necessary for
FASTQ trimming/filtering, **MultiQC** for report aggregation, and the ChIP/ATAC
workflow retains Picard insert-size QC and deepTools fingerprinting. Long-read
RNA-seq retains Sequali; the scRNA-seq and de novo workflows are unchanged.
See the [RustQC CLI](https://github.com/seqeralabs/RustQC/blob/v0.2.1/src/cli.rs)
for the supported inputs and analyses.

## Docker image used in the workflow

- [quay.io/biocontainers/fastp:0.23.4--hadf994f_2](https://quay.io/repository/biocontainers/fastp)
- [quay.io/biocontainers/star:2.7.11a--h0033a41_0](https://quay.io/repository/biocontainers/star)
- [quay.io/biocontainers/samtools:1.18--h50ea8bc_1](https://quay.io/repository/biocontainers/samtools)
- [quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1](https://quay.io/repository/biocontainers/deeptools)
- [ghcr.io/seqeralabs/rustqc:v0.2.1](https://github.com/seqeralabs/RustQC/releases/tag/v0.2.1)
- [multiqc/multiqc:v1.35](https://hub.docker.com/r/multiqc/multiqc)

## Testing

With Snakemake, pytest, PyYAML and pysam installed, run `pytest tests/rustqc -q`
from the repository root. For the integration tests, also provide RustQC 0.2.1,
samtools and MultiQC 1.35 on `PATH`, then run:

``` bash
RUSTQC_INTEGRATION=1 pytest tests/rustqc -q
```

The tests check paired/single and grouped-sample RNA/ChIP DAGs, invalid settings,
and missing annotations. Integration tests execute duplicate marking, RustQC
and MultiQC on synthetic BAMs, check read counts and sample separation, and
verify that original BAMs are unchanged. They do not run FASTQ alignment or
validate biological results on full-sized datasets. The GitHub Actions
`Test RustQC migration` job runs both DAG and integration tests.
