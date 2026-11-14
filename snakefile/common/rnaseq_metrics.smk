# RustQC requires duplicate-marked reads. Keep all reads and use a separate
# temporary BAM so downstream analyses still receive the original STAR BAM.
rule markdup_rna_qc:
    container:
        CONTAINERS["samtools"]
    input:
        bam = "star/{sample}/{sample}_Aligned.out.bam"
    output:
        bam = temp("rustqc_bam/{sample}.markdup.bam")
    params:
        extra_threads = lambda wildcards, threads: max(0, threads - 1)
    threads:
        4
    benchmark:
        "benchmark/markdup_rna_qc_{sample}.txt"
    log:
        "log/markdup_rna_qc_{sample}.log"
    shell:
        r"""
        exec > {log:q} 2>&1
        qc_tmp=$(mktemp -d rustqc_bam/.markdup.XXXXXXXX)
        trap 'rm -rf "$qc_tmp"' EXIT
        samtools sort -n -@ {params.extra_threads} -T "$qc_tmp/name" \
            -o "$qc_tmp/name.bam" {input.bam:q}
        samtools fixmate -m -@ {params.extra_threads} "$qc_tmp/name.bam" "$qc_tmp/fixmate.bam"
        rm "$qc_tmp/name.bam"
        samtools sort -@ {params.extra_threads} -T "$qc_tmp/coord" \
            -o "$qc_tmp/coord.bam" "$qc_tmp/fixmate.bam"
        rm "$qc_tmp/fixmate.bam"
        samtools markdup -@ {params.extra_threads} "$qc_tmp/coord.bam" {output.bam:q}
        """


rule rustqc_rna:
    container:
        CONTAINERS["rustqc"]
    input:
        bam = "rustqc_bam/{sample}.markdup.bam",
        gtf = gtf
    output:
        report = directory("rustqc/{sample}")
    params:
        paired = "--paired" if layout == "paired" else "",
        stranded = rustqc_stranded
    threads:
        8
    benchmark:
        "benchmark/rustqc_{sample}.txt"
    log:
        "log/rustqc_{sample}.log"
    shell:
        "rustqc rna {input.bam:q} --gtf {input.gtf:q} "
        "{params.paired} --stranded {params.stranded:q} --threads {threads} "
        "--sample-name {wildcards.sample:q} --outdir {output.report:q} "
        "--json-summary {output.report:q}/rustqc_summary.json > {log:q} 2>&1"
