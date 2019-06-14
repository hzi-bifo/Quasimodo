rule rmdup:
    input:
        rules.bwa.output.sortedbam
    output:
        rmdupbam = seq_dir + "/bam/{sample}.{ref_name}.rmdup.bam"
    conda:
        "../config/conda_env.yaml"
    log:
        report_dir + "/picard/{sample}.{ref_name}.rmdup.metrics.txt"
    threads: threads
    shell:
        """
        picard MarkDuplicates REMOVE_DUPLICATES=true \
            INPUT={input} \
            OUTPUT={output} \
            METRICS_FILE={log}
        samtools index {output.rmdupbam}
        """
