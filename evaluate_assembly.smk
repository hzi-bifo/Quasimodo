include: "rules/load_config.smk"

assembly_dir = "/".join([project_dir, "results/asssembly"])
metaquast_dir = "/".join([project_dir, "results/metaquast"])
assemblers = ["spades", "tadpole", "megahit", "ray", "idba", "abyss", "savage"]
metaquast_criteria = ["num_contigs", "Largest_contig", "Genome_fraction", 
                   "Duplication_ratio", "Largest_alignment", "LGA50", 
                   "NGA50", "num_misassemblies", "num_mismatches_per_100_kbp"]

## Get current working directory
cwd = os.getcwd()

## Samples to corresponding TM or TM mixture folders
def make_mix():
    return ["{}/{}".format(sample.split("-")[0], sample) for sample in sample_list]

onsuccess:
    print("The assembly evaluation is done!")
#    shell("mail -s 'The The assembly evaluation is done' youremail@provider.com")

rule all:
    input:
        expand("{assemblyDir}/{assembler}/{sample}.{assembler}.scaffolds.fa", 
            assemblyDir=assembly_dir, sample=list(samples["sample"]), assembler=assemblers),
        expand("{metaquastDir}/{strain_sample}/report.html", metaquastDir=metaquast_dir, 
            assembler=assemblers, strain_sample=make_mix()),
        expand(metaquast_dir + "/summary_for_figure/{mix}.{criteria}_merged.tsv",
            mix=["TM", "TA"], criteria=metaquast_criteria),
        expand(results_dir + "/final_figures/assembly_metaquast_evaluation.pdf"),
        expand(results_dir + "/final_tables/assembly_metaquast_ranking.tsv")
        

## Build index for reference
include: "rules/index.smk"

## Run all consensus assembly tools
include: "rules/assembly.smk"

## Merge read paired ends
rule pear:
    input:
        reads = get_fastq,
        assembly_done = expand("{assemblyDir}/{assembler}/{sample}.{assembler}.scaffolds.fa", 
            assemblyDir=assembly_dir, sample=list(samples["sample"]), 
            assembler=[i for i in assemblers if i!="savage"])
    output:
        merged = os.path.abspath(seq_dir + "/pear_merge/{sample}.pear.assembled.fastq"),
        unmerged_r1 = os.path.abspath(seq_dir + "/pear_merge/{sample}.pear.unassembled.forward.fastq"),
        unmerged_r2 = os.path.abspath(seq_dir + "/pear_merge/{sample}.pear.unassembled.reverse.fastq")
    conda:
        "config/conda_savage_env.yaml"
    params:
        out = seq_dir + "/pear_merge/{sample}.pear"
    threads: 16
    shell:
        """
        pear -j {threads} -f {input.reads[0]} -r {input.reads[1]} -o {params.out}
        """

## The original script random_split_fastq.py provided by Savage is not suitable for PE fq files converted
### by BAMToFastq
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


## Run Savage for haplotype reconstruction
rule savage_full_ref:
    input:
        merged = rules.pear.output.merged,
        unmerged_r1 = rules.pear.output.unmerged_r1,
        unmerged_r2 = rules.pear.output.unmerged_r2,
        full_ref = lambda w: merlin_ref if w.sample.startswith("TM") else ad_ref,
        full_ref_fai = lambda w: merlin_ref + ".fai" if w.sample.startswith("TM") else ad_ref + ".fai",
        modified = rules.modify_savage.output
    output:
        assembly_dir + "/savage/{sample}/contigs_stage_c.fasta"
    conda:
        "config/conda_savage_env.yaml"
    benchmark:
        report_dir + "/benchmarks/{sample}.savage.benchmark.txt"
    params:
        cwd = cwd,
        savage_sample_dir = assembly_dir + "/savage/{sample}"
    threads: 16
    shell:
        """
        mkdir -p {params.savage_sample_dir}
        cd {params.savage_sample_dir}
        savage -t {threads} --ref {input.full_ref} -m 150 --split 4 --merge_contigs 0 \
            --contig_len_stage_c 50 -s {input.merged} -p1 {input.unmerged_r1} -p2 {input.unmerged_r2}
        cd {params.cwd}
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

## Evaluate assemblies using metaquast
rule metaquast:
    input:
        scaffolds = lambda wc: expand(assembly_dir + "/{assembler}/{{sample}}.{assembler}.scaffolds.fa", 
            assembler=assemblers),
        ref_fai = lambda wc: [tb_ref + ".fai", ad_ref + ".fai"] if 
                wc.mix == "TA" else [tb_ref + ".fai", merlin_ref + ".fai"]
    output:
        report = metaquast_dir + "/{mix}/{sample, [A-Z]+-[0-9\-]+}/report.html"
    conda:
        "config/conda_env.yaml"
    params:
        metaquast_outdir = metaquast_dir + "/{mix}/{sample}",
        ref = lambda wc: ",".join([tb_ref, ad_ref]) if wc.mix == "TA" else ",".join([tb_ref, merlin_ref])
    shell:
        """
        metaquast.py --unique-mapping -o {params.metaquast_outdir} -R {params.ref} {input.scaffolds}
        """

## Summarize all evaluations
rule summarize:
    input:
        expand("{metaquastDir}/{strain_sample}/report.html", metaquastDir=metaquast_dir, 
            assembler=assemblers, strain_sample=make_mix())
    output:
        metaquast_dir + "/summary_for_figure/{mix}.{criteria}_merged.tsv"
    conda:
        "config/conda_env.yaml"
    params:
        input_files = metaquast_dir + "/{mix}/*/summary/TSV/{criteria}.tsv"
    shell:
        """
        paste -d"\t" {params.input_files}|sed '1s/\.scaffolds//g' |csvtk transpose -Tt -|\
            awk 'NR==1{{print}}$1!="Assemblies"{{print}}'|sed '1s/\.GFP\|\.BAC//g' > {output}
        """

## Visualize the evaluation
rule visualize:
    input:
        expand(metaquast_dir + "/summary_for_figure/{mix}.{criteria}_merged.tsv", 
            mix=["TA", "TM"], criteria=metaquast_criteria)
    output:
        figure = results_dir + "/final_figures/assembly_metaquast_evaluation.pdf",
        table = results_dir + "/final_tables/assembly_metaquast_ranking.tsv"
    conda:
        "config/conda_env.yaml"
    params:
        input_dir = metaquast_dir + "/summary_for_figure"
    script:
        "scripts/metaquast_visualize.R"

