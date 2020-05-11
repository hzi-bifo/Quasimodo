r1 = rules.rm_phix.output.cl_r1 if not rm_human_ecoli else rules.rm_human_ecoli.output.cl_r1
r2 = rules.rm_phix.output.cl_r2 if not rm_human_ecoli else rules.rm_human_ecoli.output.cl_r2

# cpus = 24

rule spades:
    input:
        r1 = r1,
        r2 = r2
    output:
        scaffolds = assembly_dir + "/spades/{sample}/scaffolds.fasta",
        renamed_scaffolds = assembly_dir + \
            "/spades/{sample}.spades.scaffolds.fa"
    conda:
        "../config/conda_env.yaml"
    params:
        outdir = assembly_dir + "/spades/{sample}"
    threads: threads
    priority: 83
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.spades.benchmark.txt"
    shell:
        """
        spades.py -k 21,33,55,77,99,127 --careful -1 {input.r1} \
            -2 {input.r2} -o {params.outdir} -t {threads}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule metaspades:
    input:
        r1 = r1,
        r2 = r2
    output:
        scaffolds = assembly_dir + "/metaspades/{sample}/scaffolds.fasta",
        renamed_scaffolds = assembly_dir + \
            "/metaspades/{sample}.metaspades.scaffolds.fa"
    conda:
        "../config/conda_env.yaml"
    params:
        outdir = assembly_dir + "/metaspades/{sample}"
    threads: threads
    priority: 85
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.metaspades.benchmark.txt"
    shell:
        """
        metaspades.py -k 21,33,55,77,99,127 -1 {input.r1} \
            -2 {input.r2} -o {params.outdir} -t {threads}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule tadpole:
    input:
        r1 = r1,
        r2 = r2
    output:
        cor_fq_r1 = assembly_dir + \
            "/tadpole/{sample}/{sample}.tadpole.corr.r1.fq",
        cor_fq_r2 = assembly_dir + \
            "/tadpole/{sample}/{sample}.tadpole.corr.r2.fq",
        scaffolds = assembly_dir + \
            "/tadpole/{sample}/{sample}.tadpole.contigs.fa",
        renamed_scaffolds = assembly_dir + \
            "/tadpole/{sample}.tadpole.scaffolds.fa"
    conda:
        "../config/conda_env.yaml"
    threads: threads
    priority: 95
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.tadpole.benchmark.txt"
    shell:
        """
        tadpole.sh in={input.r1} in2={input.r2} out={output.cor_fq_r1} \
                      out2={output.cor_fq_r2} mode=correct threads={threads}
        tadpole.sh in={output.cor_fq_r1} in2={output.cor_fq_r2} \
                      out={output.scaffolds} threads={threads}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

# TODO sometimes broken in snakemake conda env
rule megahit:
    input:
        r1 = r1,
        r2 = r2
    output:
        scaffolds = assembly_dir + \
            "/megahit/{sample}/{sample}.megahit.contigs.fa",
        renamed_scaffolds = assembly_dir + \
            "/megahit/{sample}.megahit.scaffolds.fa"
    conda:
        "../config/conda_env.yaml"
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.megahit.benchmark.txt"
    params:
        megahit_out = assembly_dir + "/megahit/{sample}",
        prefix = "{sample}.megahit"
    threads: threads
    priority: 100
    # rmdir --ignore-fail-on-non-empty {params.megahit_out}
    shell:
        """
        rm -rf {params.megahit_out}
        megahit -t {threads} --k-min 21 --k-max 151 -1 {input.r1} \
            -2 {input.r2} -o {params.megahit_out}  --out-prefix {params.prefix}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """
        # rm -rf {params.megahit_out}

rule ray:
    input:
        r1 = r1,
        r2 = r2
    output:
        scaffolds = assembly_dir + "/ray/{sample}/Scaffolds.fasta",
        renamed_scaffolds = assembly_dir + "/ray/{sample}.ray.scaffolds.fa"
    conda:
        "../config/conda_env.yaml"
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.ray.benchmark.txt"
    params:
        ray_out = assembly_dir + "/ray/{sample}"
    threads: threads
    priority: 20
    shell:
        """
        rm -rf {params.ray_out}
        mpiexec -n {threads} Ray -k31 -p {input.r1} {input.r2} -o {params.ray_out}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule idba:
    input:
        r1 = r1,
        r2 = r2
    output:
        scaffolds = assembly_dir + "/idba/{sample}/scaffold.fa",
        renamed_scaffolds = assembly_dir + "/idba/{sample}.idba.scaffolds.fa"
    conda:
        "../config/conda_env.yaml"
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.idba.benchmark.txt"
    params:
        merged_pe = seq_dir + "/cl_fq/{sample}_12.fa",
        idba_out = assembly_dir + "/idba/{sample}"
    threads: threads
    priority: 70
    shell:
        """
        rm -rf {params.idba_out}
        fq2fa --merge {input.r1} {input.r2} {params.merged_pe}
        idba_ud -r {params.merged_pe} --num_threads {threads} -o {params.idba_out}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule abyss:
    input:
        r1 = r1,
        r2 = r2
    output:
        scaffolds = assembly_dir + \
            "/abyss/{sample}/{sample}.abyss-scaffolds.fa",
        renamed_scaffolds = assembly_dir + "/abyss/{sample}.abyss.scaffolds.fa"
    conda:
        "../config/conda_env.yaml"
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.abyss.benchmark.txt"
    params:
        abyss_out = "{sample}.abyss",
        abyss_outdir = assembly_dir + "/abyss/{sample}"
    threads: threads
    priority: 90
    shell:
        """
        abyss-pe np={threads} name={params.abyss_out} k=96 in='{input.r1} {input.r2}'
        mv {params.abyss_out}* {params.abyss_outdir}
        cp {output.scaffolds} {output.renamed_scaffolds}
        """

rule iva:
    input:
        r1 = seq_dir + "/cl_fq/{sample}.qc.cl.r1.fq",
        r2 = seq_dir + "/cl_fq/{sample}.qc.cl.r2.fq"
    output:
        contigs = assembly_dir + "/iva/{sample}/contig.fasta",
        scaffolds = assembly_dir + "/iva/{sample}.iva.scaffolds.fa"
    benchmark:
        report_dir + "/benchmarks/assembler/assembly/{sample}.iva.benchmark.txt"
    priority: 82
    params:
        out_dir = assembly_dir + "/iva/{sample}"
    conda:
        "../config/conda_iva.yaml"
    threads: threads
    shell:
        """
        rm -rf {params.out_dir}
        iva -f {input.r1} -r {input.r2} -t {threads} {params.out_dir}
        cp {output.contigs} {output.scaffolds}
        """


rule mk_fq_dir:
    input:
        r1 = seq_dir + "/cl_fq/{sample}.qc.cl.r1.fq",
        r2 = seq_dir + "/cl_fq/{sample}.qc.cl.r2.fq"
    output:
        or1 = seq_dir + "/cl_fq/{sample}/{sample}.qc.cl.r1.fq",
        or2 = seq_dir + "/cl_fq/{sample}/{sample}.qc.cl.r2.fq"
    shell:
        """
        cp {input.r1} {output.or1}
        cp {input.r2} {output.or2}
        """


rule vicuna:
    input:
        r1 = seq_dir + "/cl_fq/{sample}/{sample}.qc.cl.r1.fq",
        r2 = seq_dir + "/cl_fq/{sample}/{sample}.qc.cl.r2.fq"
    output:
        config = assembly_dir + "/vicuna/{sample}/vicuna.conf",
        contigs = assembly_dir + "/vicuna/{sample}/contig.fasta",
        scaffolds = assembly_dir + "/vicuna/{sample}.vicuna.scaffolds.fa"
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.vicuna.benchmark.txt"
    priority: 82
    params:
        fq_dir = seq_dir + "/cl_fq/{sample}/",
        out_dir = assembly_dir + "/vicuna/{sample}/"
    threads: threads
    shell:
        """
        python program/vicuna_config_generator.py {params.fq_dir} {params.out_dir} > {output.config}
        export OMP_NUM_THREADS={threads}
        ./libs/VICUNA_v1.3/executable/vicuna-omp.static.linux64 {output.config}
        cp {output.contigs} {output.scaffolds}
        """


# Merge read paired ends
rule pear:
    input:
        r1 = rules.rm_phix.output.cl_r1,
        r2 = rules.rm_phix.output.cl_r2,
    
    output:
        merged = os.path.abspath(
            seq_dir + "/pear_merge/{sample}.pear.assembled.fastq"),
        unmerged_r1 = os.path.abspath(
            seq_dir + "/pear_merge/{sample}.pear.unassembled.forward.fastq"),
        unmerged_r2 = os.path.abspath(
            seq_dir + "/pear_merge/{sample}.pear.unassembled.reverse.fastq")
    conda:
        "config/conda_savage_env.yaml"
    params:
        out = seq_dir + "/pear_merge/{sample}.pear"
    threads: threads
    shell:
        """
        pear -j {threads} -f {input.r1} -r {input.r2} -o {params.out}
        """

# The original script random_split_fastq.py provided by Savage is not suitable for PE fq files converted
# by BAMToFastq
rule modify_savage:
    output:
        "modify_savage.done"

    conda:
        "config/conda_savage_env.yaml"
    shell:
        """
        path_savage=${{PATH%%/bin:*}}
        mv $path_savage/opt/savage-0.4.0/scripts/random_split_fastq.{{py,bak}}
        cp program/random_split_fastq.py $path_savage/opt/savage-0.4.0/scripts/random_split_fastq.py
        touch {output}
        """

# Run Savage for haplotype reconstruction
rule savage_full_ref:
    input:
        merged = rules.pear.output.merged,
        unmerged_r1 = rules.pear.output.unmerged_r1,
        unmerged_r2 = rules.pear.output.unmerged_r2,
        full_ref = lambda w: merlin_ref if w.sample.startswith(
            "TM") else ad_ref,
        full_ref_fai = lambda w: merlin_ref + \
            ".fai" if w.sample.startswith("TM") else ad_ref + ".fai",
        modified = rules.modify_savage.output
    output:
        assembly_dir + "/savage/{sample}/contigs_stage_c.fasta"
    conda:
        "config/conda_savage_env.yaml"
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.savage.benchmark.txt"
    params:
        cd = cd,
        savage_sample_dir = assembly_dir + "/savage/{sample}"
    threads: 24
    shell:
        """
        mkdir -p {params.savage_sample_dir}
        cd {params.savage_sample_dir}
        savage -t {threads} --ref {input.full_ref} -m 150 --split 4 --merge_contigs 0 \
            --contig_len_stage_c 50 -s {input.merged} -p1 {input.unmerged_r1} -p2 {input.unmerged_r2}
        cd {params.cd}
        """

rule rename_savage:
    input:
        rules.savage_full_ref.output
    output:
        assembly_dir + "/savage/{sample}.savage.scaffolds.fa"
    shell:
        """
        cp {input} {output}
        """


# rule pehaplo:
#     input:
#         r1 = os.path.abspath(r1),
#         r2 = os.path.abspath(r2)
#     output:
#         contigs = assembly_dir + "/pehaplo/{sample}/Contigs.fa",
#         scaffolds = assembly_dir + "/pehaplo/{sample}.pehaplo.scaffolds.fa"
#     benchmark:
#         report_dir + "/benchmarks/assembly/{sample}.pehaplo.benchmark.txt"
#     params:
#         out_dir = assembly_dir + "/pehaplo/{sample}",
#         cd = cd
#     threads: cpus
#     priority: 80
#     # mkdir -p {params.out_dir}
#     shell:
#         """
#         export PATH=/usr/bin:{params.cd}/libs/PEHaplo/bin:/home/zldeng/miniconda3/envs/assembly/bin:/home/zldeng/miniconda3/bin:$PATH
        
#         cd {params.out_dir}
#         python2 {params.cd}/libs/PEHaplo/pehaplo.py -f1 {input.r1} -f2 {input.r2} \
#             -l 180 -l1 210 -r 300 -F 600 -std 150 -n 3 -correct yes -t {threads}
        
#         cd {params.cd}
#         cp {output.contigs} {output.scaffolds}
#         """

# Reference based assembly
rule virgena:
    input:
        r1 = os.path.abspath(r1),
        r2 = os.path.abspath(r2),
        ref = ref_seq
    output:
        contigs = assembly_dir + "/virgena/{sample}/contigs.fasta",
        # scaffolds = assembly_dir + "/virgena/{sample}.virgena.scaffolds.fa"
    benchmark:
        report_dir + "/benchmarks/assembler/{sample}.virgena.benchmark.txt"
    priority: 72
    params:
        out_dir = assembly_dir + "/virgena/{sample}"
    threads: 24
    shell:
        """
        python program/virgena_config_generator.py {input.r1} {input.r2} {input.ref[0]} {params.out_dir} \
            -t {threads} > {params.out_dir}/config.xml
        java -jar libs/virgena/VirGenA.jar assemble -c {params.out_dir}/config.xml
        """