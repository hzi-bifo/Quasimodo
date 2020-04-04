# Call INDEls only work with version 2.1.3.1
rule lofreq:
    input:
        rmdupbam = rules.rmdup.output.rmdupbam, 
        ref = ref_seq
    output:
        indel_quality_bam = temp(seq_dir + "/bam/{sample}.{ref_name}.rmdup.indelq.bam"),
        indel_quality_bai = temp(seq_dir + "/bam/{sample}.{ref_name}.rmdup.indelq.bam.bai"),
        vcf = snpcall_dir + "/lofreq/{sample}.{ref_name}.lofreq.vcf",
        vcf_bgz = snpcall_dir + "/lofreq/{sample}.{ref_name}.lofreq.vcf.gz"
    benchmark:
        report_dir + "/benchmarks/{sample}.{ref_name}.lofreq.benchmark.txt"
    conda:
        "../config/conda_env.yaml"
    threads: threads
    shell:
        """
        lofreq indelqual --dindel -f {input.ref[0]} {input.rmdupbam} -o {output.indel_quality_bam}
        samtools index {output.indel_quality_bam}
        lofreq call-parallel --pp-threads {threads} --call-indels -q 20 -Q 20 -m 20 \
            -f {input.ref[0]} {output.indel_quality_bam} -o {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}
        """

rule mpileup:
    input:
        rmdupbam = rules.rmdup.output.rmdupbam, 
        ref = ref_seq
    output:
        seq_dir + "/pileup/{sample}.{ref_name}.mpileup"
    benchmark:
        report_dir + "/benchmarks/{sample}.{ref_name}.mpileup.benchmark.txt"
        
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        samtools mpileup -f {input.ref[0]} {input.rmdupbam} > {output} 
        """

rule varscan:
    input:
        # rmdupbam = rules.rmdup.output.rmdupbam, 
        mpileup = rules.mpileup.output,
        ref = ref_seq
    output:
        # snp_vcf = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.snp.vcf",
        # indel_vcf = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.indel.vcf",
        pileup2snp_var = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.pileup2snp.var",
        pileup2indel_var = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.pileup2indel.var",
        # vcf = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.vcf",
        # vcf_bgz = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.vcf.gz"
    benchmark:
        report_dir + "/benchmarks/{sample}.{ref_name}.varscan.benchmark.txt"
    conda:
        "../config/conda_env.yaml"
    shell:
        """        
        varscan pileup2snp {input.mpileup} --min-avg-qual 20 --p-value 0.01 > {output.pileup2snp_var}
        varscan pileup2indel {input.mpileup} --min-avg-qual 20 --p-value 0.01 > {output.pileup2indel_var}
        """

rule varscan2vcf:
    input:
        snp = rules.varscan.output.pileup2snp_var,
        indel = rules.varscan.output.pileup2indel_var
    output:
        vcf = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.vcf",
        vcf_bgz = snpcall_dir + "/varscan/{sample}.{ref_name}.varscan.vcf.gz",
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        python3 program/varscan2vcf.py -s {input.snp} -i {input.indel} -o {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}
        """

rule freebayes:
    input:
        rmdupbam = rules.rmdup.output.rmdupbam, 
        ref = ref_seq
    output:
        raw_vcf = snpcall_dir + "/freebayes/{sample}.{ref_name}.freebayes.raw.vcf",
        vcf = snpcall_dir + "/freebayes/{sample}.{ref_name}.freebayes.vcf",
        vcf_bgz = snpcall_dir + "/freebayes/{sample}.{ref_name}.freebayes.vcf.gz"
    benchmark:
        report_dir + "/benchmarks/{sample}.{ref_name}.freebayes.benchmark.txt"
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        freebayes -p 1 -m 20 -q 20 -F 0.01 --min-coverage 10 -f {input.ref[0]} \
            {input.rmdupbam} > {output.raw_vcf}
        bcftools norm -m -both {output.raw_vcf} | vt  decompose_blocksub - -o {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}  
        """

rule bcftools:
    input:
        rmdupbam = rules.rmdup.output.rmdupbam, 
        ref = ref_seq
    output:
        vcf = snpcall_dir + "/bcftools/{sample}.{ref_name}.bcftools.vcf",
        vcf_bgz = snpcall_dir + "/bcftools/{sample}.{ref_name}.bcftools.vcf.gz"
    benchmark:
        report_dir + "/benchmarks/{sample}.{ref_name}.bcftools.benchmark.txt"
    threads: threads
    conda:
        "../config/conda_env.yaml"
    shell:
        """
        bcftools mpileup --threads {threads} -Ou -f {input.ref[0]} {input.rmdupbam} |\
            bcftools call --threads {threads} -p 0.01 --ploidy 1 -mv -Ob |\
            bcftools view -i 'INFO/DP>=10' - > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}  
        """


# add read group information for GATK Mutect2 and Haplotypecaller
rule add_read_group:
    input:
        rules.rmdup.output.rmdupbam,
    output:
        seq_dir + "/bam/{sample}.{ref_name}.rmdup.RG.bam"

    shell:
        """
        gatk AddOrReplaceReadGroups -I {input} \
            -O {output} \
            -RGLB lib1 \
            -RGPL ILLUMINA \
            -RGPU unit1 \
            -RGSM {wildcards.sample}
        samtools index {output}
        """

# run the variants calling with Mutect2
rule mutect2:
    input:
        rmdupbam = rules.add_read_group.output,
        ref = ref_seq
    output:
        vcf = snpcall_dir + "/mutect2/{sample}.{ref_name}.mutect2.vcf",
        vcf_bgz = snpcall_dir + "/mutect2/{sample}.{ref_name}.mutect2.vcf.gz"
    benchmark:
        report_dir + "/benchmarks/{sample}.{ref_name}.mutect2.benchmark.txt"
    shell:
        """
        gatk Mutect2 --input {input.rmdupbam} \
            --reference {input.ref[0]} --output {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz} 
        """


# GATK haplotypecaller
rule gatk:
    input:
        rmdupbam = rules.add_read_group.output, 
        ref = ref_seq
    output:
        vcf = snpcall_dir + "/gatk/{sample}.{ref_name}.gatk.vcf",
        vcf_bgz = snpcall_dir + "/gatk/{sample}.{ref_name}.gatk.vcf.gz"
    benchmark:
        report_dir + "/benchmarks/{sample}.{ref_name}.gatk.benchmark.txt"
    # threads: threads
    # conda:
    #     "../config/conda_env.yaml"
    shell:
        """
        gatk HaplotypeCaller --input {input.rmdupbam} --reference {input.ref[0]} \
            --min-base-quality-score 20 -ploidy 1 --output {output.vcf} 
        bgzip -c {output.vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}  
        """