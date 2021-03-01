rule extractTP:
    input:
        vcf = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.vcf",
        genome_diff = lambda wc: snp_dir + "/nucmer/{mix}.maskrepeat.variants.vcf".format(mix=wc.sample[:2])
    output:
        filtered = snpcall_dir + \
            "/{snpcaller}/{sample}.{ref}.{snpcaller}.filtered.vcf",
        # tp = snpcall_dir + \
        #     "/{snpcaller}/tp/{sample}.{ref}.{snpcaller}.tp.vcf",
        fp = snpcall_dir + \
            "/{snpcaller}/fp/{sample}.{ref}.{snpcaller}.fp.vcf"
    params:
        data = "hcmv",
        outdir = snpcall_dir + "/{snpcaller}"
    conda:
        "../config/conda_env.yaml"
    threads: threads
    shell:
        """
        python program/extract_TP_FP_SNPs.py {input.vcf} {input.genome_diff} {params.data} {params.outdir}
        """
