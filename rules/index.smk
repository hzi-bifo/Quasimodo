rule build_idx:
    input:
        ref = "{ref_seq}"
    output:
        fa_idx = "{ref_seq}.fai",
        bwa_idx = "{ref_seq}.bwt"
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        samtools faidx {input.ref}
        bwa index {input.ref}
        """
