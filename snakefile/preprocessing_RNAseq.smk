
'''
Snakefile for RNA-seq preprocessing (paired-end and single-end)

Usage:
    snakemake -s preprocessing_RNAseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity

config.yaml must contain:
    layout: "paired" or "single"
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
star_index = config["star_index"]
gtf = config["gtf"]
layout = config["layout"]
if layout not in ("paired", "single"):
    raise ValueError("layout must be 'paired' or 'single'")
rustqc_stranded = config.get("rustqc_stranded", "unstranded")
if rustqc_stranded not in ("unstranded", "forward", "reverse"):
    raise ValueError("rustqc_stranded must be 'unstranded', 'forward', or 'reverse'")

include: "common/containers.smk"
include: "common/functions.smk"

samples, samples_dict = normalize_samples(config["samples"])

include: "common/sort_bam.smk"
include: "common/bigwig.smk"
include: "common/rnaseq_metrics.smk"

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

# When samples is given as a dict {sample: [run1, run2, ...]}, fastq files for
# each run are concatenated into fastq_merged/ before QC. When samples is a
# plain list, the merge rule is not emitted and rule qc reads fastq/ directly.
def _qc_input_paired(wildcards):
    base = "fastq_merged" if samples_dict else "fastq"
    return {
        "R1": f"{base}/{wildcards.sample}_1.fastq.gz",
        "R2": f"{base}/{wildcards.sample}_2.fastq.gz",
    }

def _qc_input_single(wildcards):
    base = "fastq_merged" if samples_dict else "fastq"
    return {"R1": f"{base}/{wildcards.sample}.fastq.gz"}

if layout == "paired":

    rule all:
        default_target: True
        input:
            multiqc = "multiqc/multiqc_report.html",
            bam = expand("star/{sample}/{sample}_Aligned.out.bam", sample = samples),
            bigwig = expand("bigwig/{sample}.bw", sample = samples),
            rustqc = expand("rustqc/{sample}", sample = samples)

    if samples_dict:

        rule merge_fastq_paired:
            input:
                R1 = lambda wc: [f"fastq/{run}_1.fastq.gz" for run in samples_dict[wc.sample]],
                R2 = lambda wc: [f"fastq/{run}_2.fastq.gz" for run in samples_dict[wc.sample]]
            output:
                R1 = "fastq_merged/{sample}_1.fastq.gz",
                R2 = "fastq_merged/{sample}_2.fastq.gz"
            benchmark:
                "benchmark/merge_fastq_{sample}.txt"
            log:
                "log/merge_fastq_{sample}.log"
            shell:
                "cat {input.R1} > {output.R1} 2> {log} && "
                "cat {input.R2} > {output.R2} 2>> {log}"

    rule qc:
        container:
            CONTAINERS["fastp"]
        input:
            unpack(_qc_input_paired)
        output:
            R1 = "fastp/{sample}_1.fastq.gz",
            R2 = "fastp/{sample}_2.fastq.gz",
            json = "fastp/log/{sample}.json"
        params:
            html = "fastp/log/{sample}.html"
        threads:
            8
        benchmark:
            "benchmark/fastp_{sample}.txt"
        log:
            "log/fastp_{sample}.log"
        shell:
            "fastp -i {input.R1} -I {input.R2} "
            "-o {output.R1} -O {output.R2} -w {threads} "
            "-h {params.html} -j {output.json} >& {log}"

    rule mapping:
        container:
            CONTAINERS["star"]
        input:
            R1 = "fastp/{sample}_1.fastq.gz",
            R2 = "fastp/{sample}_2.fastq.gz",
        output:
            sam = "star/{sample}/{sample}_Aligned.out.sam",
            log = "star/{sample}/{sample}_Log.final.out"
        params:
            outdir = "star/{sample}/{sample}_"
        threads:
            workflow.cores
        benchmark:
            "benchmark/star_{sample}.txt"
        log:
            "log/star_{sample}.log"
        shell:
            "STAR --runThreadN {threads} --genomeDir {star_index} "
            "--outFilterMultimapNmax 1 "
            "--readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --outFileNamePrefix {params.outdir} >& {log}"


else:

    rule all:
        default_target: True
        input:
            multiqc = "multiqc/multiqc_report.html",
            bam = expand("star/{sample}/{sample}_Aligned.out.bam", sample = samples),
            bigwig = expand("bigwig/{sample}.bw", sample = samples),
            rustqc = expand("rustqc/{sample}", sample = samples)

    if samples_dict:

        rule merge_fastq_single:
            input:
                R1 = lambda wc: [f"fastq/{run}.fastq.gz" for run in samples_dict[wc.sample]]
            output:
                R1 = "fastq_merged/{sample}.fastq.gz"
            benchmark:
                "benchmark/merge_fastq_{sample}.txt"
            log:
                "log/merge_fastq_{sample}.log"
            shell:
                "cat {input.R1} > {output.R1} 2> {log}"

    rule qc:
        container:
            CONTAINERS["fastp"]
        input:
            unpack(_qc_input_single)
        output:
            R1 = "fastp/{sample}.fastq.gz",
            json = "fastp/log/{sample}.json"
        params:
            html = "fastp/log/{sample}.html"
        threads:
            8
        benchmark:
            "benchmark/fastp_{sample}.txt"
        log:
            "log/fastp_{sample}.log"
        shell:
            "fastp -i {input.R1} "
            "-o {output.R1} -w {threads} "
            "-h {params.html} -j {output.json} >& {log}"

    rule mapping:
        container:
            CONTAINERS["star"]
        input:
            R1 = "fastp/{sample}.fastq.gz"
        output:
            sam = "star/{sample}/{sample}_Aligned.out.sam",
            log = "star/{sample}/{sample}_Log.final.out"
        params:
            outdir = "star/{sample}/{sample}_"
        threads:
            workflow.cores
        benchmark:
            "benchmark/star_{sample}.txt"
        log:
            "log/star_{sample}.log"
        shell:
            "STAR --runThreadN {threads} --genomeDir {star_index} "
            "--outFilterMultimapNmax 1 "
            "--readFilesIn {input.R1} --readFilesCommand zcat --outFileNamePrefix {params.outdir} >& {log}"

rule multiqc:
    container:
        CONTAINERS["multiqc"]
    input:
        json = expand("fastp/log/{sample}.json", sample = samples),
        starlog = expand("star/{sample}/{sample}_Log.final.out", sample = samples),
        rustqc = expand("rustqc/{sample}", sample = samples),
        config = workflow.source_path("common/multiqc_rna.yaml")
    output:
        "multiqc/multiqc_report.html"
    benchmark:
        "benchmark/multiqc.txt"
    log:
        "log/multiqc.log"
    shell:
        "multiqc --force --config {input.config:q} -o multiqc/ "
        "{input.json:q} {input.starlog:q} {input.rustqc:q} > {log:q} 2>&1"

## DON'T CHANGE ABOVE THIS LINE ##
