# The selected callers for comparison in Venn digram
venn_snpcallers = ["lofreq", "varscan", "clc"]

rule snp_evaluate:
    input:
        mix = expand(snpcall_dir + "/{snpcallers}/{sample_ref}.{snpcallers}.filtered.vcf",
                     sample_ref=['{}.{}'.format(sample, sample_refname_dict[sample])
                                 for sample in sample_list if not sample.endswith(("-1-0", "-0-1"))],
                     snpcallers=snpcallers),
        diff = expand(
            snp_dir + "/nucmer/{genome_diff}.maskrepeat.snps", genome_diff=genome_diff_list)
    output:
        snpcaller_performance_summary = results_dir + \
            "/final_tables/snpcaller_performance_summary.txt",
        snpcaller_performance_figure_pdf = results_dir + \
            "/final_figures/snpcaller_performance.pdf",
    conda:
        "../config/conda_env.yaml"
    params:
        mix_sample = list(
            filter(lambda x: not x.endswith(("-1-0", "-0-1")), sample_list)),
        venn_snpcallers = venn_snpcallers
    script:
        "../scripts/snpcaller_performance_compare.R"
