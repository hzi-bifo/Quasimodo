rule bwa:
    input:
        reads = get_fastq,
        ref = ref_seq
#        ref_bwa_idx = ref_seq_idx
    output:
        sortedbam = seq_dir + "/bam/{sample}.{ref_name}.bam"
    conda:
        "config/conda_env.yaml"
    log:
        report_dir + "/bwa/{sample}.{ref_name}.log"
    threads: 16
    shell:
        """
        bwa mem -k 31 -t {threads} {input.ref[0]} {input[0]} {input[1]} |\
            samtools view -Shb - 2>> {log} |\
            samtools sort -@ {threads} -m 10G - -o {output.sortedbam} >> {log} 2>&1
        """
