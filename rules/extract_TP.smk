rule extractTP:
    input:
        vcf = snpcall_dir + "/{snpcallers}/{sample}.{ref}.{snpcallers}.vcf",
        genome_diff = get_genome_diff
    output:
        filtered = snpcall_dir + \
            "/{snpcallers}/{sample}.{ref}.{snpcallers}.filtered.vcf",
        # tp = snpcall_dir + \
        #     "/{snpcallers}/tp/{sample}.{ref}.{snpcallers}.tp.vcf",
        fp = snpcall_dir + \
            "/{snpcallers}/fp/{sample}.{ref}.{snpcallers}.fp.vcf"
    params:
        data = "hcmv",
        outdir = snpcall_dir + "/{snpcallers}"
    conda:
        "../config/conda_env.yaml"
    threads: threads
    shell:
        """
        python program/extract_TP_FP_SNPs.py {input.vcf} {input.genome_diff} {params.data} {params.outdir}
        """
