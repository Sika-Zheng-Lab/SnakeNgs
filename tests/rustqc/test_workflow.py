"""DAG regressions and real BAM -> duplicate marking -> RustQC -> MultiQC tests.

Integration tests need rustqc 0.2.1, samtools, MultiQC 1.35 and pysam on PATH.
Set RUSTQC_INTEGRATION=1 to require these dependencies instead of skipping.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SNAKEMAKE = shutil.which("snakemake")


def run_workflow(config, *args, success=True):
    result = subprocess.run(
        [SNAKEMAKE, "-s", str(ROOT / "snakefile" / config["snakefile"]),
         "--configfile", str(config["path"]), "--cores", "2", *args],
        cwd=ROOT, capture_output=True, text=True, timeout=180,
    )
    output = result.stdout + result.stderr
    if success:
        assert result.returncode == 0, output
    else:
        assert result.returncode != 0, output
    return output


def fixture_config(tmp_path, assay="RNAseq", layout="paired", grouped=False):
    work = tmp_path / "work"
    (work / "fastq").mkdir(parents=True)
    index = tmp_path / "index"
    index.mkdir()
    (index / "chrLength.txt").write_text("10000\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\texon\t101\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_biotype "protein_coding";\n'
        'chr1\ttest\texon\t2001\t3000\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_biotype "protein_coding";\n'
        'chr1\ttest\texon\t4001\t6000\t.\t-\t.\tgene_id "g2"; transcript_id "t2"; gene_biotype "rRNA";\n'
    )
    # Similar names expose report/sample-name collisions.
    samples = ["sample.one", "sample.two"]
    runs = ["lane1", "lane2", "lane3"] if grouped else samples
    for run in runs:
        suffixes = ["_1", "_2"] if layout == "paired" else [""]
        for suffix in suffixes:
            (work / "fastq" / f"{run}{suffix}.fastq.gz").touch()
    cfg = dict(workdir=str(work), samples=samples, layout=layout,
               star_index=str(index), gtf=str(gtf),
               bowtie2_index=str(index / "genome"), bowtie2_args="")
    if grouped:
        cfg["samples"] = {samples[0]: runs[:2], samples[1]: runs[2:]}
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return dict(path=path, cfg=cfg, work=work, samples=samples,
                snakefile=f"preprocessing_{assay}.smk")


@pytest.mark.parametrize("assay", ["RNAseq", "ChIPseq"])
@pytest.mark.parametrize("layout", ["paired", "single"])
@pytest.mark.parametrize("grouped", [False, True])
def test_full_dag(tmp_path, assay, layout, grouped):
    config = fixture_config(tmp_path, assay, layout, grouped)
    out = run_workflow(config, "all", "--dry-run", "--printshellcmds")
    assert "multiqc" in out
    assert ("rule merge_fastq_" in out) == grouped
    if assay == "RNAseq":
        assert "rule markdup_rna_qc:" in out and "rule rustqc_rna:" in out
        assert "CollectRnaSeqMetrics" not in out and "CollectInsertSizeMetrics" not in out
        commands = [line for line in out.splitlines() if "rustqc rna " in line]
        assert len(commands) == 2
        assert all(("--paired" in command) == (layout == "paired") for command in commands)
        assert all("--threads 2" in command for command in commands)
        # No explicit target must also schedule the full pipeline.
        default = run_workflow(config, "--dry-run")
        assert "rule rustqc_rna:" in default and "rule multiqc:" in default
    else:
        # RustQC has no ChIP mode: preserve insert sizes, fingerprints and trimming.
        assert "rule plotFingerprint:" in out and "fastp -i" in out
        assert ("rule CollectInsertSizeMetrics:" in out) == (layout == "paired")
        assert "rustqc rna " not in out


@pytest.mark.parametrize("key,value", [("layout", "typo"), ("rustqc_stranded", "typo")])
def test_invalid_config(tmp_path, key, value):
    config = fixture_config(tmp_path)
    config["cfg"][key] = value
    config["path"].write_text(yaml.safe_dump(config["cfg"]))
    assert f"{key} must be" in run_workflow(config, "--dry-run", success=False)


def test_missing_gtf(tmp_path):
    config = fixture_config(tmp_path)
    Path(config["cfg"]["gtf"]).unlink()
    assert "genes.gtf" in run_workflow(config, "all", "--dry-run", success=False)


def make_bam(path, paired, spliced=False):
    import pysam
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 10000}]}
    reads = []
    for i in range(240):
        # Each locus has three copies; the fixture must acquire duplicate flags.
        start = 150 + (i // 3) * 4 if i < 120 else 4100 + ((i - 120) // 3) * 4
        for mate in range(2 if paired else 1):
            read = pysam.AlignedSegment()
            read.query_name = f"read{i}"
            read.query_sequence = "ACGT" * 25
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.flag = (99 if mate == 0 else 147) if paired else 0
            read.reference_id = 0
            read.reference_start = start + 150 * mate
            read.mapping_quality = 60
            read.cigarstring = "100M"
            if paired:
                read.next_reference_id = 0
                read.next_reference_start = start + (150 if mate == 0 else 0)
                read.template_length = 250 if mate == 0 else -250
            if spliced and i < 6:
                read.reference_start = 950 if mate == 0 else 2100
                read.cigarstring = "50M1000N50M" if mate == 0 else "100M"
                if paired:
                    read.next_reference_start = 2100 if mate == 0 else 950
                    read.template_length = 1250 if mate == 0 else -1250
            read.set_tag("NH", 1)
            read.set_tag("NM", 0)
            reads.append(read)
    with pysam.AlignmentFile(path, "wb", header=header) as handle:
        for read in sorted(reads, key=lambda r: r.reference_start):
            handle.write(read)
    return len(reads)


@pytest.mark.parametrize("layout", ["paired", "single"])
def test_real_qc_and_multiqc(tmp_path, layout):
    missing = [name for name in ("rustqc", "samtools", "multiqc") if not shutil.which(name)]
    if missing:
        if os.environ.get("RUSTQC_INTEGRATION") == "1":
            pytest.fail(f"Missing integration dependencies: {missing}")
        pytest.skip(f"Missing integration dependencies: {missing}")
    import pysam
    config = fixture_config(tmp_path, layout=layout)
    config["cfg"]["rustqc_stranded"] = "reverse"
    config["path"].write_text(yaml.safe_dump(config["cfg"]))
    work = config["work"]
    checksums = {}
    for sample in config["samples"]:
        bam = work / f"star/{sample}/{sample}_Aligned.out.bam"
        total = make_bam(bam, layout == "paired", spliced=sample == "sample.two")
        checksums[sample] = hashlib.sha256(bam.read_bytes()).hexdigest()
        (bam.parent / f"{sample}_Log.final.out").write_text(
            "Number of input reads | 240\nUniquely mapped reads number | 240\n"
            "Uniquely mapped reads % | 100.00%\n"
            "Number of reads mapped to multiple loci | 0\n% of reads mapped to multiple loci | 0.00%\n"
            "Number of reads mapped to too many loci | 0\n% of reads mapped to too many loci | 0.00%\n"
            "% of reads unmapped: too many mismatches | 0.00%\n"
            "% of reads unmapped: too short | 0.00%\n% of reads unmapped: other | 0.00%\n"
        )
        fastp = work / f"fastp/log/{sample}.json"
        fastp.parent.mkdir(parents=True, exist_ok=True)
        summary = dict(total_reads=total, total_bases=total * 100, q20_bases=total * 100,
                       q30_bases=total * 100, q20_rate=1.0, q30_rate=1.0,
                       read1_mean_length=100, read2_mean_length=100, gc_content=0.5)
        fastp.write_text(json.dumps(dict(
            summary=dict(fastp_version="0.23.4", before_filtering=summary, after_filtering=summary),
            command="fastp -i input.fastq.gz", filtering_result=dict(passed_filter_reads=total,
                low_quality_reads=0, too_many_N_reads=0, too_short_reads=0),
            duplication=dict(rate=0.0),
        ), indent=2))
    run_workflow(config, "multiqc/multiqc_report.html", "--notemp", "--printshellcmds",
                 "--allowed-rules", "markdup_rna_qc", "rustqc_rna", "multiqc")
    for sample in config["samples"]:
        bam = work / f"star/{sample}/{sample}_Aligned.out.bam"
        assert hashlib.sha256(bam.read_bytes()).hexdigest() == checksums[sample]
        with pysam.AlignmentFile(work / f"rustqc_bam/{sample}.markdup.bam") as handle:
            reads = list(handle)
        assert len(reads) == total
        assert sum(r.is_duplicate for r in reads) > 0
        report = work / "rustqc" / sample
        assert (report / "rustqc_summary.json").stat().st_size > 0
        stats = report / "samtools" / f"{sample}.stats"
        assert "SN\traw total sequences:\t" + str(total) in stats.read_text()
        flagstat = report / "samtools" / f"{sample}.flagstat"
        assert flagstat.read_text().startswith(f"{total} + 0 in total")
        assert (report / "featurecounts" / f"{sample}.featureCounts.tsv.summary").exists()
    assert (work / "multiqc/multiqc_report.html").stat().st_size > 1000
    sources = (work / "multiqc/multiqc_data/multiqc_sources.txt").read_text()
    for sample in config["samples"]:
        assert sample in sources
        assert f"{sample}.stats" in sources
        assert f"{sample}.featureCounts.tsv.summary" in sources
    assert "fastp" in sources.lower() and "star" in sources.lower()
    assert "sample.one.junction_annotation.log" not in sources
    assert "sample.two.junction_annotation.log" in sources
    assert ".featureCounts.biotype.tsv.summary" not in sources
    biotypes = (work / "multiqc/multiqc_data/multiqc_featurecounts_biotype_plot.txt").read_text()
    assert "sample.one" in biotypes and "sample.two" in biotypes
    # The skipped zero-junction report remains available for inspection.
    assert (work / "rustqc/sample.one/rseqc/junction_annotation/sample.one.junction_annotation.log").exists()
