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


if rm_human_ecoli:
    human_ref = "/net/sgi/viral_genomics/MHH/human_genome/hg19.genome.bwa"
    ecoli_ref = "ncbi_ref/Ecoli.NC_000913.fa"

    rule rm_human_ecoli:
        input:
            r1 = rules.rm_phix.output.cl_r1, 
            r2 = rules.rm_phix.output.cl_r2,
            human = human_ref,
            human_ref_idx = human_ref + ".bwt",
            ecoli = ecoli_ref,
            ecoli_ref_idx = ecoli_ref + ".bwt"
        output:
            nohuman_r1 = seq_dir + "/cl_fq/{sample}.qc.nohuman.r1.fq",
            nohuman_r2 = seq_dir + "/cl_fq/{sample}.qc.nohuman_r1.r1.fq",
            cl_r1 = seq_dir + "/cl_fq/{sample}.qc.cl.r1.fq",
            cl_r2 = seq_dir + "/cl_fq/{sample}.qc.cl.r2.fq"
        wildcard_constraints:
            sample = "[^\.]+"
        threads: 16
        shell:
            """
            bwa mem -k 31 -t {threads} {input.human} {input.r1} {input.r2}|samtools view -Shb - |\
                samtools view -b -f 12 -F 256 -|samtools sort -n - -@ {threads} -m 10G |\
                bedtools bamtofastq -i - -fq {output.nohuman_r1} -fq2 {output.nohuman_r2}
            bwa mem -k 31 -t {threads} {input.ecoli} {output.nohuman_r1} {output.nohuman_r2}|samtools view -Shb - |\
                samtools view -b -f 12 -F 256 -|samtools sort -n - -@ {threads} -m 10G |\
                bedtools bamtofastq -i - -fq {output.cl_r1} -fq2 {output.cl_r2}
            """