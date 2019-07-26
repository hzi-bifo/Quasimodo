import os

configfile: "config/customize_data.yaml"

refs = list(map(str.strip, config["refs"].split(",")))
project_dir = config["projectDir"].rstrip("/")
threads = config["threads"]
results_dir = project_dir + "/results"

vcfs = list(map(str.strip, config["vcfs"].split(",")))
snp_dir = results_dir + "/snp"
snpcall_dir = snp_dir + "/callers"
g1_name, g2_name = (os.path.splitext(os.path.basename(ref))[0] for ref in refs)

gdiff_name = g1_name + "_" + g2_name

snpcallers = [os.path.splitext(os.path.basename(vcf))[0] for vcf in vcfs]


def get_vcf(wc):
    for vcf_file in vcfs:
        if vcf_file.endswith(wc.snpcaller + ".vcf"):
            return vcf_file


rule all:
    input:
        snp_benchmark_figure = results_dir + "/final_figures/snpcall_benchmark.pdf",
        snp_venn_figure = results_dir + "/final_figures/snpcall_venn.pdf",
        snp_benchmark_table = results_dir + "/final_tables/snpcall_benchmark.txt"

# The first given ref should be the ref used to generate the VCFs
rule gdiff:
    input:
        refs
    output:
        delta = snp_dir + "/nucmer/" + gdiff_name + ".delta",
        snps = snp_dir + "/nucmer/" + gdiff_name + ".maskrepeat.snps"
    conda:
        "config/conda_env.yaml"
    params:
        genome_diff_prefix = snp_dir + "/nucmer/" + gdiff_name
    shell:
        """
        nucmer --prefix={params.genome_diff_prefix} {input}
        show-snps -CTIHlr <(delta-filter -r -q {output.delta}) > {output.snps}
        """

rule extract_TP:
    input:
        vcf = get_vcf,
        genome_diff = rules.gdiff.output.snps
    output:
        filtered = snpcall_dir + "/{snpcaller}.filtered.vcf",
        fp = snpcall_dir + "/fp/{snpcaller}.fp.vcf"
    params:
        data = "custom",
        outdir = snpcall_dir
    conda:
        "config/conda_env.yaml"
    threads: threads
    shell:
        """
        python program/extract_TP_FP_SNPs.py {input.vcf} {input.genome_diff} {params.data} {params.outdir}
        """

rule snp_benchmark:
    input:
        vcfs = expand(snpcall_dir + "/{snpcaller}.filtered.vcf",
                      snpcaller=snpcallers),
        genome_diff = rules.gdiff.output.snps
    output:
        snp_benchmark_table = results_dir + "/final_tables/snpcall_benchmark.txt",
        snp_benchmark_figure = results_dir + "/final_figures/snpcall_benchmark.pdf",
        snp_venn_figure = results_dir + "/final_figures/snpcall_venn.pdf"
    params: snpcallers = snpcallers
    conda:
        "config/conda_env.yaml"
    script:
        "scripts/custom_snp_benchmark.R"
