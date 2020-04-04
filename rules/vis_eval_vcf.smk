# The selected callers for comparison in Venn digram
venn_snpcallers = ["lofreq", "varscan", "clc"]


rule rtg:
    input:
        ref = ref_seq,
        gs_vcf = lambda wc: snp_dir + "/nucmer/{mix}.maskrepeat.variants.vcf.gz".format(mix=wc.sample[:2]),
        caller_vcf = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.vcf.gz"
    output:
        snp_dir + "/rtg/{snpcaller, [^\.\/]+}/{sample}.{ref, [^\.\/]+}/summary.txt",
        snp_dir + "/rtg/{snpcaller}/{sample}.{ref, [^\.\/]+}/weighted_roc.tsv.gz"
    conda:
        "../config/conda_env.yaml"
    params:
        rtg_outdir = snp_dir + "/rtg/{snpcaller}/{sample}.{ref}"
    shell:
        """
        rm -rf {params.rtg_outdir}
        rtg vcfeval -b {input.gs_vcf} -c {input.caller_vcf}  -o {params.rtg_outdir} \
            -t {input.ref[2]} --squash-ploidy -f QUAL --sample ALT
        """


rule extract_snp:
    input:
        variants = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.vcf",
    output:
        vcf = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.xsnp.vcf",
        vcf_bgz = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.xsnp.vcf.gz",
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        awk -F"\t" '/^#.*/{{print}}$4~/^[actgACTG]$/&&$5~/^[actgACTG]$/' {input.variants} > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}
        """

rule extract_indel:
    input:
        variants = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.vcf",
    output:
        vcf = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.xindel.vcf",
        vcf_bgz = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.xindel.vcf.gz",
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        awk -F"\t" '/^#.*/{{print}}$4~/^[actgACTG]{{2,}}/||$5~/^[actgACTG]{{2,}}/' {input.variants} > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}
        """

rule extract_nucmer_snp:
    input:
        variants = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.variants.vcf"

    output:
        vcf = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.xsnp.vcf",
        vcf_bgz = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.xsnp.vcf.gz"
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        awk -F"\t" '/^#.*/{{print}}$4~/^[actgACTG]$/&&$5~/^[actgACTG]$/' {input.variants} > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}
        """

rule extract_nucmer_indel:
    input:
        variants = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.variants.vcf"
    output:
        vcf = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.xindel.vcf",
        vcf_bgz = snp_dir + "/nucmer/{genome_diff_name}.maskrepeat.xindel.vcf.gz"
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        awk -F"\t" '/^#.*/{{print}}$4~/^[actgACTG]{{2,}}/||$5~/^[actgACTG]{{2,}}/' {input.variants} > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}
        """

rule rtg_snp:
    input:
        ref = ref_seq,
        gs_vcf = lambda wc: snp_dir + "/nucmer/{mix}.maskrepeat.xsnp.vcf.gz".format(mix=wc.sample[:2]),
        caller_vcf = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.xsnp.vcf.gz"
    output:
        snp_dir + "/rtg/{snpcaller}/{sample}.{ref}.xsnp/summary.txt",
        snp_dir + "/rtg/{snpcaller}/{sample}.{ref}.xsnp/weighted_roc.tsv.gz"
    conda:
        "../config/conda_env.yaml"
    params:
        rtg_outdir = snp_dir + "/rtg/{snpcaller}/{sample}.{ref}.xsnp"
    shell:
        """
        rm -rf {params.rtg_outdir}
        rtg vcfeval -b {input.gs_vcf} -c {input.caller_vcf}  -o {params.rtg_outdir} \
            -t {input.ref[2]} --squash-ploidy -f QUAL --sample ALT
        """

rule rtg_indel:
    input:
        ref = ref_seq,
        gs_vcf = lambda wc: snp_dir + "/nucmer/{mix}.maskrepeat.xindel.vcf.gz".format(mix=wc.sample[:2]),
        caller_vcf = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.xindel.vcf.gz"
    output:
        snp_dir + "/rtg/{snpcaller}/{sample}.{ref}.xindel/summary.txt"
    conda:
        "../config/conda_env.yaml"
    params:
        rtg_outdir = snp_dir + "/rtg/{snpcaller}/{sample}.{ref}.xindel"
    shell:
        """
        rm -rf {params.rtg_outdir}
        rtg vcfeval -b {input.gs_vcf} -c {input.caller_vcf}  -o {params.rtg_outdir} \
            -t {input.ref[2]} --squash-ploidy -f QUAL --sample ALT
        """



rule snp_evaluate:
    input:
        snp = expand(snpcall_dir + "/{snpcaller}/{sample_ref}.{snpcaller}.filtered.vcf",
                    sample_ref=sample_ref,
                    snpcaller=snpcallers),
        diff = expand(
            snp_dir + "/nucmer/{genome_diff}.maskrepeat.variants.vcf", genome_diff=genome_diff_list),
        rtg_summary = expand(snp_dir + "/rtg/{snpcaller}/{sample_ref}/summary.txt",
                     sample_ref=mixed_strain_sample_ref,
                     snpcaller=snpcallers),
        rtg_roc = expand(snp_dir + "/rtg/{snpcaller}/{sample_ref}/weighted_roc.tsv.gz",
                     sample_ref=mixed_strain_sample_ref,
                     snpcaller=snpcallers),
        rtg_snp_roc = expand(snp_dir + "/rtg/{snpcaller}/{sample_ref}.xsnp/weighted_roc.tsv.gz",
                     sample_ref=mixed_strain_sample_ref,
                     snpcaller=snpcallers),
        rtg_snp_summary = expand(snp_dir + "/rtg/{snpcaller}/{sample_ref}.xsnp/summary.txt",
                     sample_ref=mixed_strain_sample_ref,
                     snpcaller=snpcallers)
    output:
        performance_table = results_dir + \
            "/final_tables/caller_performance.tsv",
        rtg_performance_table = results_dir + \
            "/final_tables/caller_rtg_performance.tsv",
        rtg_snp_performance_table = results_dir + \
            "/final_tables/caller_rtg_snp_performance.tsv",

        rtg_score20_performance_table = results_dir + \
            "/final_tables/caller_rtg_score20_performance.tsv",
        rtg_snp_score20_performance_table = results_dir + \
            "/final_tables/caller_rtg_snp_score20_performance.tsv",


        rtg_roc_table = results_dir + "/final_tables/caller_rtg_roc_performance.tsv",
        rtg_snp_roc_table = results_dir + "/final_tables/caller_rtg_snp_roc_performance.tsv",

        
        rtg_roc_figure = results_dir + "/final_figures/caller_rtg_recall_precision_curve.pdf",
        rtg_snp_roc_figure = results_dir + "/final_figures/caller_rtg_snp_recall_precision_curve.pdf",
        # fp_tp_pdf = results_dir + "/final_figures/rtg_fp_tp_curve.pdf",
        
        
        performance_figure = results_dir + \
            "/final_figures/caller_performance.pdf",
        rtg_performance_figure = results_dir + \
            "/final_figures/caller_rtg_bestf1snp_score20all_score20snp_performance.pdf",
        snp_venndiagram_figure = results_dir + \
            "/final_figures/caller_snp_venndiagram.pdf"
    conda:
        "../config/conda_env.yaml"
    params:
        mix_sample = list(
            filter(lambda x: not x.endswith(("-1-0", "-0-1")), sample_list)),
        venn_snpcallers = venn_snpcallers
    script:
        "../scripts/caller_performance_compare.R"




