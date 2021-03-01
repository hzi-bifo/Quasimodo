fp_compared_snpcallers = ["lofreq", "clc", "varscan", "freebayes"]

rule compareFP:
    input:
        fp = expand(snpcall_dir + "/{snpcaller}/fp/{sample_ref}.{snpcaller}.fp.vcf",
                    snpcaller=snpcallers,
                    sample_ref=['{}.{}'.format(sample, sample_refname_dict[sample])
                                for sample in sample_list if not sample.endswith(("-1-0", "-0-1"))])
    output:
        fp_compare_figure = results_dir + "/final_figures/snpcaller_fp_snp_compare.pdf",
#        fp_compare_table = results_dir + "/final_tables/snpcaller_fp_snp_compare.txt"
    conda:
        "../config/conda_env.yaml"
    params:
        mix_sample = list(
            filter(lambda x: not x.endswith(("-1-0", "-0-1")), sample_list)),
        fp_compared_snpcallers = fp_compared_snpcallers
    script:
        "../scripts/snpcaller_fp_compare.R"
