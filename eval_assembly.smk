# Load parameters from config file
include: "rules/load_config.smk"

assembly_dir = "/".join([project_dir, "results/assembly"])
metaquast_dir = "/".join([project_dir, "results/metaquast"])
assemblers = ["spades", "metaspades", "tadpole", "abyss", 
              "megahit", "ray", "idba", "vicuna", "iva", "savage"]  # "haploflow", "pehaplo", "quasirecomb",

metaquast_criteria = ["num_contigs", "Largest_contig", "Genome_fraction",
                      "Duplication_ratio", "Largest_alignment", "LGA50",
                      "NGA50", "num_mismatches_per_100_kbp"]

# Get current working directory
# cwd = os.getcwd()

# Samples to corresponding TM or TM mixture folders

# ruleorder: metaspades > megahit > tadpole > abyss > spades > ray > idba > savage_full_ref

def make_mix():
    return ["{}/{}".format(sample.split("-")[0], sample) for sample in sample_list]

def get_assembly_ref(wc):
    if wc.sample.startswith('TA'):
        if wc.sample.endswith('-1-0'):
            ref_list = [tb_ref]
        elif wc.sample.endswith('-0-1'):
            ref_list = [ad_ref]
        else:
            ref_list = [tb_ref, ad_ref]
    else:
        if wc.sample.endswith('-1-0'):
            ref_list = [tb_ref]
        elif wc.sample.endswith('-0-1'):
            ref_list = [merlin_ref]
        else:
            ref_list = [tb_ref, merlin_ref]
    return ','.join(ref_list)


onsuccess:
    print("The assembly evaluation is done!")
#    shell("mail -s 'The assembly evaluation is done' youremail@provider.com")

rule all:
    input:
        metaquast_report = expand(metaquast_dir + "/{strain_sample}/report.html",
                strain_sample=make_mix()),
        all_sample_metaquast_table = results_dir + "/final_tables/all_sample_metaquast.tsv",
        figure = results_dir + "/final_figures/assembly_metaquast_evaluation.pdf"
        # expand("{assemblyDir}/{assembler}/{sample}.{assembler}.scaffolds.fa",
        #        assemblyDir=assembly_dir, sample=sample_list, assembler=assemblers),
        
        # expand(metaquast_dir + "/summary_for_figure/{mix}.{criteria}.merged.tsv",
        #        mix=["TM", "TA"], criteria=metaquast_criteria),
        # results_dir + "/final_figures/assembly_metaquast_evaluation.pdf",
        # results_dir + "/final_tables/assembly_metaquast_ranking.tsv",
        # results_dir + "/final_tables/assembly_metaquast_score.tsv",
        # results_dir + "/final_tables/assembly_metaquast_scaled.tsv"


# Build index for reference
include: "rules/index.smk"

# If not run on reads, copy the resulting scaffolds provided within the software for benchmarking
if not run_on_reads:
    # Extract assembly
    if not os.path.exists(cd + "/data/assembly"):
        shell("tar -xzvf {} -C {}".format(cd + "/data/assembly.tar.gz",
                                          cd + "/data/"))

    rule cp_assembly:
        input: cd + \
            "/data/assembly/{assembler}/{sample}.{assembler}.scaffolds.fa"
        output: assembly_dir + "/{assembler}/{sample}.{assembler}.scaffolds.fa"
        shell:
            """
            cp {input} {output}
            """

else:
    # Remove remaining Phix reads
    include: "rules/decontamination.smk"

    # Run all consensus assembly tools
    include: "rules/assembly.smk"


# Evaluate assemblies using metaquast
rule metaquast:
    input:
        scaffolds = lambda wc: expand(assembly_dir + "/{assembler}/{{sample}}.{assembler}.scaffolds.fa",
                                      assembler=assemblers),
        ref_fai = lambda wc: [tb_ref + ".fai", ad_ref + ".fai"] if
        wc.mix == "TA" else [tb_ref + ".fai", merlin_ref + ".fai"]
    output:
        report = metaquast_dir + "/{mix}/{sample, [A-Z]+-[0-9\-]+}/report.html",
        tsv_report = metaquast_dir + "/{mix}/{sample, [A-Z]+-[0-9\-]+}/combined_reference/report.tsv"

    conda:
        "config/conda_env.yaml"
    threads: threads
    params:
        metaquast_outdir = metaquast_dir + "/{mix}/{sample}",
        ref = get_assembly_ref
        # ref = lambda wc: ",".join(
        #     [tb_ref, ad_ref]) if wc.mix == "TA" else ",".join([tb_ref, merlin_ref])
    shell:
        """
        metaquast.py --unique-mapping -o {params.metaquast_outdir} -R {params.ref} {input.scaffolds} -t {threads}
        """

# Summarize all evaluations
rule summarize:
    input:
        expand(metaquast_dir + "/{strain_sample}/report.html",
               strain_sample=make_mix())
    output:
        metaquast_dir + "/summary_for_figure/{mix}.{criteria}.merged.tsv"
    # conda:
    #     "config/conda_env.yaml"
    params:
        input_files = metaquast_dir + "/{mix}/*/summary/TSV/{criteria}.tsv",
        joiner = cd + '/program/join_tsv.py'
    shell:
        """
        python {params.joiner} {params.input_files}|sed '1s/\.scaffolds//g' |csvtk transpose -Tt -|\
            awk 'NR==1{{print}}$1!="Assemblies"{{print}}'|sed '1s/\.GFP\|\.BAC//g' > {output}
        """

#    paste -d"\t" {params.input_files}|sed '1s/\.scaffolds//g' |csvtk transpose -Tt -|\
#             awk 'NR==1{{print}}$1!="Assemblies"{{print}}'|sed '1s/\.GFP\|\.BAC//g' > {output}

# Visualize the evaluation
rule visualize:
    input:
        individual_ref_reports = expand(metaquast_dir + "/summary_for_figure/{mix}.{criteria}.merged.tsv",
                mix=["TA", "TM"], criteria=metaquast_criteria),
        combined_ref_reports = expand("{metaquastDir}/{strain_sample}/combined_reference/report.tsv", metaquastDir=metaquast_dir,
               assembler=assemblers, strain_sample=make_mix())
    output:
        figure = results_dir + "/final_figures/assembly_metaquast_evaluation.pdf",
        all_sample_table = results_dir + "/final_tables/all_sample_metaquast.tsv",
        table = results_dir + "/final_tables/assembly_metaquast_ranking.tsv",
        radarplot_table = results_dir + "/final_tables/assembly_metaquast_scaled.tsv",
        score = results_dir + "/final_tables/assembly_metaquast_score.tsv"
    conda:
        "config/conda_env.yaml"
    params:
        input_dir = metaquast_dir + "/summary_for_figure"
    script:
        "scripts/metaquast_visualize.R"
