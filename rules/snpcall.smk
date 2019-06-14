rule lofreq:
    input:
        rmdupbam = rules.rmdup.output.rmdupbam, 
        ref = ref_seq
    output:
        vcf = snpcall_dir + "/lofreq/{sample}.{ref_name}.lofreq.vcf"
    conda:
        "config/conda_env.yaml"
    threads: 16
    shell:
        """
        lofreq call-parallel --pp-threads {threads} -q 20 -Q 20 -m 20 \
            -f {input.ref[0]} {input.rmdupbam} -o {output.vcf}
        """

rule varscan:
    input:
        rmdupbam = rules.rmdup.output.rmdupbam, 
        ref = ref_seq
    output:
        pileup = seq_dir + "/pileup/{sample}.onvirus.{ref_name}.pileup",
        var = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.out",
        vcf = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.vcf"
    conda:
        "config/conda_env.yaml"
    threads: 16
    shell:
        """
        samtools mpileup -f {input.ref[0]} {input.rmdupbam} > {output.pileup} 
        varscan pileup2snp {output.pileup} --min-avg-qual 20 --p-value 0.01 > {output.var}
        python program/varscan2vcf.py {output.var} > {output.vcf}
        """

rule freebayes:
    input:
        rmdupbam = rules.rmdup.output.rmdupbam, 
        ref = ref_seq
    output:
        vcf = snpcall_dir + "/freebayes/{sample}.{ref_name}.freebayes.vcf"
    threads: 16
    conda:
        "config/conda_env.yaml"
    shell:
        """
        freebayes -p 1 -m 20 -q 20 -F 0.01 --min-coverage 10 -f {input.ref[0]} \
            {input.rmdupbam} | awk -F"\t" '/^#/{{print}}$6>=20{{print}}' > {output.vcf}
        """

rule bcftools:
    input:
        rmdupbam = rules.rmdup.output.rmdupbam, 
        ref = ref_seq
    output:
        vcf = snpcall_dir + "/bcftools/{sample}.{ref_name}.bcftools.vcf"
    threads: 16
    conda:
        "config/conda_env.yaml"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref[0]} {input.rmdupbam} |\
            bcftools call -p 0.01 --ploidy 1 -mv -Ob |\
            bcftools view -i '%QUAL>=20' - > {output.vcf}
        """
