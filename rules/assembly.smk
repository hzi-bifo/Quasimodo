rule spades:
    input:
        get_fastq
    output:
        scaffolds = assembly_dir + "/spades/{sample}/scaffolds.fasta",
        renamed_scaffolds = assembly_dir + "/spades/{sample}.spades.scaffolds.fa"
    params:
        outdir = assembly_dir + "/spades/{sample}"
    threads: 16
    benchmark:
        report_dir + "/benchmarks/{sample}.spades.benchmark.txt"
    shell:
        """
        spades.py -k 21,33,55,77,99,127 --careful -1 {input[0]} \
            -2 {input[1]} -o {params.outdir} -t {threads}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule tadpole:
    input:
        get_fastq
    output:
        cor_fq_r1 = assembly_dir + "/tadpole/{sample}/{sample}.tadpole.corr.r1.fq",
        cor_fq_r2 = assembly_dir + "/tadpole/{sample}/{sample}.tadpole.corr.r2.fq",
        scaffolds = assembly_dir + "/tadpole/{sample}/{sample}.tadpole.contigs.fa",
        renamed_scaffolds = assembly_dir + "/tadpole/{sample}.tadpole.scaffolds.fa"
    threads: 16
    benchmark:
        report_dir + "/benchmarks/{sample}.tadpole.benchmark.txt"
    shell:
        """
        tadpole.sh in={input[0]} in2={input[1]} out={output.cor_fq_r1} \
                      out2={output.cor_fq_r2} mode=correct threads={threads}
        tadpole.sh in={output.cor_fq_r1} in2={output.cor_fq_r2} \
                      out={output.scaffolds} threads={threads}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule megahit:
    input:
        get_fastq
    output:
        scaffolds = assembly_dir + "/megahit/{sample}/{sample}.megahit.contigs.fa",
        renamed_scaffolds = assembly_dir + "/megahit/{sample}.megahit.scaffolds.fa"
    benchmark:
        report_dir + "/benchmarks/{sample}.megahit.benchmark.txt"
    params:
        megahit_out = assembly_dir + "/megahit/{sample}",
        prefix = "{sample}.megahit"
    threads: 16
    shell:
        """
        rmdir {params.megahit_out}
        megahit -t {threads} --continue --k-min 21 --k-max 151 -1 {input[0]} \
            -2 {input[1]} -o {params.megahit_out}  --out-prefix {params.prefix}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule ray:
    input:
        get_fastq
    output:
        scaffolds = assembly_dir + "/ray/{sample}/Scaffolds.fasta",
        renamed_scaffolds = assembly_dir + "/ray/{sample}.ray.scaffolds.fa"
    benchmark:
        report_dir + "/benchmarks/{sample}.ray.benchmark.txt"
    params:
        ray_out = assembly_dir + "/ray/{sample}"
    threads: 16
    shell:
        """
        rmdir {params.ray_out}
        mpiexec -n {threads} Ray -k31 -p {input[0]} {input[1]} -o {params.ray_out}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule idba:
    input:
        get_fastq
    output:
        scaffolds = assembly_dir + "/idba/{sample}/scaffold.fa",
        renamed_scaffolds = assembly_dir + "/idba/{sample}.idba.scaffolds.fa"
    benchmark:
        report_dir + "/benchmarks/{sample}.idba.benchmark.txt"
    params:
        merged_pe = seq_dir + "/qc_fq/{sample}_12.fa",
        idba_out = assembly_dir + "/idba/{sample}"
    threads: 16
    shell:
        """
        rmdir {params.idba_out}
        fq2fa --merge {input[0]} {input[1]} {params.merged_pe}
        idba_ud -r {params.merged_pe} --num_threads {threads} -o {params.idba_out}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule abyss:
    input:
        get_fastq
    output:
        scaffolds = assembly_dir + "/abyss/{sample}/{sample}.abyss-scaffolds.fa",
        renamed_scaffolds = assembly_dir + "/abyss/{sample}.abyss.scaffolds.fa"
    benchmark:
        report_dir + "/benchmarks/{sample}.abyss.benchmark.txt"
    params:
        abyss_out = "{sample}.abyss",
        abyss_outdir = assembly_dir + "/abyss/{sample}"
    threads: 16
    shell:
        """
        abyss-pe np={threads} name={params.abyss_out} k=96 in='{input[0]} {input[1]}'
        mv {params.abyss_out}* {params.abyss_outdir}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

