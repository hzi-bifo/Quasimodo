rule extractTP:
    input:
        vcf = snpcall_dir + "/{snpcallers}/{sample}.{ref}.{snpcallers}.vcf",
        genome_diff = get_genome_diff
    output:
        tp = snpcall_dir + "/{snpcallers}/tp/{sample}.{ref}.{snpcallers}.tp.vcf",
        fp = snpcall_dir + "/{snpcallers}/fp/{sample}.{ref}.{snpcallers}.fp.vcf"
    conda:
        "config/conda_env.yaml"
    shell:
        """
        python program/extract_TP_FP_SNPs.py {input.vcf} {input.genome_diff} {output.tp} {output.fp}
        """

