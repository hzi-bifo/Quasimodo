include: "rules/load_config_custom.smk"
try:
    scaffolds = list(map(str.strip, config["scaffolds"].split(",")))
except AttributeError as e:
    raise PathNotGiven(
        "The scaffold files from assembly are not specified.")

metaquast_dir = results_dir + "/metaquast"
metaquast_criteria = ["num_contigs", "Largest_contig", "Genome_fraction",
                      "Duplication_ratio", "Largest_alignment", "LGA50",
                      "NGA50", "num_misassemblies", "num_mismatches_per_100_kbp"]

rule all:
    input:
        assembly_benchmark_figure = results_dir + \
            "/final_figures/assembly_benchmark.pdf",
      #  assembly_rank_table = results_dir + "/final_tables/assembly_rank.txt"

rule build_idx:
    input:
        refs
    output:
        fa_idx = [ref + ".fai" for ref in refs]
    conda:
        "config/conda_env.yaml"
    shell:
        """
        for fa in {input}
        do
            samtools faidx $fa
        done
        """

rule metaquast:
    input:
        refs = refs,
        scaffolds = scaffolds,
        refs_fai = [ref + ".fai" for ref in refs]
    output:
        report = metaquast_dir + "/report.html",
        table = expand(metaquast_dir + "/summary/TSV/{criteria}.tsv",
                       criteria=metaquast_criteria)
    conda:
        "config/conda_env.yaml"
    threads: threads
    params:
        metaquast_outdir = metaquast_dir,
        refs_comb = config["refs"]
    shell:
        """
        metaquast.py --unique-mapping -o {params.metaquast_outdir} -R {params.refs_comb} \
            {input.scaffolds} -t {threads}
        """

# Visualize the evaluation
rule visualize:
    input:
        expand(metaquast_dir + "/summary/TSV/{criteria}.tsv",
               criteria=metaquast_criteria)
    output:
        assembly_benchmark_figure = results_dir + \
            "/final_figures/assembly_benchmark.pdf",
       # assembly_rank_table = results_dir + "/final_tables/assembly_rank.txt"
    params: number_ref = len(refs)
    conda:
        "config/conda_env.yaml"
    script:
        "scripts/custom_assembly_benchmark.R"
