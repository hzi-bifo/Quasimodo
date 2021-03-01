rule nucmer:
    input:
        ref = lambda x: merlin_ref if x.genome_diff_name.startswith("TM") else ad_ref,
        qry = tb_ref 
    output:
        delta = temp(snp_dir + "/nucmer/{genome_diff_name}.delta"),
        #no_repeat_snp = snp_dir + "/nucmer/{genome_diff_name}.norepeat.snps",
        # mask_repeat_snp = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.snps",
        mask_repeat_variants = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.variants",
        variants_vcf = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.variants.vcf",
        variants_vcf_bgz = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.variants.vcf.gz"
    conda:
        "../config/conda_env.yaml"
    params:
        genome_diff_prefix = snp_dir + "/nucmer/{genome_diff_name}"
    shell:
        ## No repeat: show-snps -CTIHr {output.delta} > {output.no_repeat_snp}
        # show-snps -CTIHlr <(delta-filter -r -q {output.delta}) > {output.mask_repeat_snp}
        """
        nucmer --prefix={params.genome_diff_prefix} {input.ref} {input.qry}
        
        
        show-snps -CTHlr <(delta-filter -r -q {output.delta}) > {output.mask_repeat_variants}
        python3 program/mummer2vcf.py -s {output.mask_repeat_variants} --output-header -n -g {input.ref} > \
            {output.variants_vcf}
        bgzip -c {output.variants_vcf} > {output.variants_vcf_bgz}
        tabix -p vcf {output.variants_vcf_bgz}
        """



rule nucmer4:
    input:
        ref = lambda x: merlin_ref if x.genome_diff_name.startswith("TM") else ad_ref,
        qry = tb_ref 
    output:
        delta = temp(snp_dir + "/nucmer4/{genome_diff_name}.delta"),
        #no_repeat_snp = snp_dir + "/nucmer4/{genome_diff_name}.norepeat.snps",
        # mask_repeat_snp = snp_dir + "/nucmer4/{genome_diff_name}.maskrepeat.snps",
        mask_repeat_variants = snp_dir + "/nucmer4/{genome_diff_name}.maskrepeat.variants",
        variants_vcf = snp_dir + "/nucmer4/{genome_diff_name}.maskrepeat.variants.vcf",
        variants_vcf_bgz = snp_dir + "/nucmer4/{genome_diff_name}.maskrepeat.variants.vcf.gz"
    # conda:
    #     "../config/conda_env.yaml"
    params:
        genome_diff_prefix = snp_dir + "/nucmer4/{genome_diff_name}"
    shell:
        ## No repeat: show-snps -CTIHr {output.delta} > {output.no_repeat_snp}
        # show-snps -CTIHlr <(delta-filter -r -q {output.delta}) > {output.mask_repeat_snp}
        """
        nucmer --prefix={params.genome_diff_prefix} {input.ref} {input.qry}
        
        
        show-snps -CTHlr <(delta-filter -r -q {output.delta}) > {output.mask_repeat_variants}
        python program/mummer2vcf.py -s {output.mask_repeat_variants} --output-header -n -g {input.ref} > \
            {output.variants_vcf}
        bgzip -c {output.variants_vcf} > {output.variants_vcf_bgz}
        tabix -p vcf {output.variants_vcf_bgz}
        """


rule minimap2:
    input:
        ref = lambda x: merlin_ref if x.genome_diff_name.startswith("TM") else ad_ref,
        qry = tb_ref 
    output:
        paf = temp(snp_dir + "/minimap2_asm20/{genome_diff_name}.paf"),        
        vcf = snp_dir + "/minimap2_asm20/{genome_diff_name}.variants.vcf",
        vcf_bgz = snp_dir + "/minimap2_asm20/{genome_diff_name}.variants.vcf.gz"
    # conda:
    #     "../config/conda_env.yaml"
    shell:
        """
        minimap2 -cx asm20 --cs {input.ref} {input.qry} > {output.paf}
        sort -k6,6 -k8,8n {output.paf} | paftools.js call -f {input.ref} - > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}
        """



