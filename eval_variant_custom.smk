include: "rules/load_config_custom.smk"

try:
    vcfs = list(map(str.strip, config["vcfs"].split(",")))
except AttributeError as e:
    raise PathNotGiven(
        "The VCF files from SNP calling are not specified.")

snp_dir = results_dir + "/snp"
snpcall_dir = snp_dir + "/callers"
g1_name, g2_name = (os.path.splitext(os.path.basename(ref))[0] for ref in refs)

gdiff_name = g1_name + "_" + g2_name

labels = config['labels']
novenn = config['novenn']
callers = [os.path.splitext(os.path.basename(vcf))[0] for vcf in vcfs] if labels is None else [label for label in labels.split(',')]

caller_vcf_dict = dict(zip(callers, vcfs))

# def get_vcf(wc):
#     return caller_vcf_dict[wc.snpcaller]
    # for vcf_file in vcfs:
    #     if vcf_file.endswith(wc.snpcaller + ".vcf"):
    #         return vcf_file


rule all:
    input:
        snp_benchmark_figure = results_dir + "/final_figures/snpcall_benchmark.pdf",
        # snp_venn_figure = results_dir + "/final_figures/snpcall_venn.pdf",
        snp_benchmark_table = results_dir + "/final_tables/snpcall_benchmark.txt"

# The first given ref should be the ref used to generate the VCFs
rule gdiff:
    input:
        refs
    output:
        delta = snp_dir + "/nucmer/" + gdiff_name + ".delta",
        snps = snp_dir + "/nucmer/" + gdiff_name + ".maskrepeat.snps"
        # mask_repeat_variants = snp_dir + "/nucmer/" + gdiff_name + ".maskrepeat.variants",
        # variants_vcf = snp_dir + "/nucmer/" + gdiff_name + ".maskrepeat.variants.vcf",
        # variants_vcf_bgz = snp_dir + "/nucmer/" + gdiff_name + ".maskrepeat.variants.vcf.gz"
    conda:
        "config/conda_env.yaml"
    params:
        genome_diff_prefix = snp_dir + "/nucmer/" + gdiff_name
    shell:
        """
        nucmer --prefix={params.genome_diff_prefix} {input}
        show-snps -CTHIlr <(delta-filter -r -q {output.delta}) > {output.snps}
        """
        # python3 program/mummer2vcf.py -s {output.mask_repeat_variants} --output-header -n -g {input[0]} > \
        #     {output.variants_vcf}
        # bgzip -c {output.variants_vcf} > {output.variants_vcf_bgz}
        # tabix -p vcf {output.variants_vcf_bgz}

rule extract_TP:
    input:
        vcf = lambda wc: caller_vcf_dict[wc.snpcaller],
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
        python program/extract_TP_FP_SNPs.py {input.vcf} {input.genome_diff} {params.data} {params.outdir} {wildcards.snpcaller}
        """

rule snp_benchmark:
    input:
        vcfs = expand(snpcall_dir + "/{snpcaller}.filtered.vcf",
                      snpcaller=callers),
        genome_diff = rules.gdiff.output.snps
    output:
        snp_benchmark_table = results_dir + "/final_tables/snpcall_benchmark.txt",
        snp_benchmark_figure = results_dir + "/final_figures/snpcall_benchmark.pdf"
        
    params: 
        callers = callers,
        novenn = novenn,
        snp_venn_figure = results_dir + "/final_figures/snpcall_venn.pdf"
    conda:
        "config/conda_env.yaml"
    script:
        "scripts/custom_snp_benchmark.R"
