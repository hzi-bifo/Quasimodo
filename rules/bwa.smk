rule bwa:
    input:
        r1 = rules.rm_phix.output.cl_r1,
        r2 = rules.rm_phix.output.cl_r2,
        ref = ref_seq
    output:
        sortedbam = seq_dir + "/bam/{sample}.{ref_name}.bam"
    conda:
        "../config/conda_env.yaml"
    log:
        report_dir + "/bwa/{sample}.{ref_name}.log"
    threads: threads
    shell:
        """
        bwa mem -k 31 -t {threads} {input.ref[0]} {input.r1} {input.r2} |\
            samtools view -Shb - 2>> {log} |\
            samtools sort -@ {threads} -m 10G - -o {output.sortedbam} >> {log} 2>&1
        """
