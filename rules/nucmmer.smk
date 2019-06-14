rule nucmmer:
    input:
        ref = lambda x: merlin_ref if x.genome_diff_name.startswith("TM") else ad_ref,
        qry = tb_ref 
    output:
        delta = snp_dir + "/nucmer/{genome_diff_name}.delta",
        no_repeat_snp = snp_dir + "/nucmer/{genome_diff_name}.norepeat.snps",
        mask_repeat_snp = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.snps"
    conda:
        "../config/conda_env.yaml"
    params:
        genome_diff_prefix = snp_dir + "/nucmer/{genome_diff_name}"
    shell:
        """
        nucmer --prefix={params.genome_diff_prefix} {input.ref} {input.qry}
        show-snps -CTIHr {output.delta} > {output.no_repeat_snp}
        show-snps -CTIHlr <(delta-filter -r -q {output.delta}) > {output.mask_repeat_snp}
        """
