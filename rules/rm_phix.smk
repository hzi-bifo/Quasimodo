rule rm_phix:
    input:
        reads = get_fastq,
        phix = phix_ref,
        phix_ref_idx = phix_ref + ".bwt"
    output:
        cl_bam = temp(seq_dir + "/cl_bam/{sample}.qc.nophix.bam"),
        cl_r1 = seq_dir + "/cl_fq/{sample}.qc.nophix.r1.fq",
        cl_r2 = seq_dir + "/cl_fq/{sample}.qc.nophix.r2.fq"
    wildcard_constraints:
        sample = "[^\.]+"
    threads: 16
    shell:
        """
        bwa mem -k 31 -t {threads} {input.phix} {input.reads[0]} {input.reads[1]}|samtools view -Shb - |\
            samtools view -b -f 12 -F 256 -|samtools sort -n - -@ {threads} -m 10G -o {output.cl_bam}
        bedtools bamtofastq -i {output.cl_bam} -fq {output.cl_r1} -fq2 {output.cl_r2}
        """
