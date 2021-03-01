include: "rules/load_config.smk"

snp_dir = "/".join([results_dir, "snp"])
snpcall_dir = "/".join([snp_dir, "callers"])

# Two mixtures
genome_diff_list = ["TM", "TA"]

# SNP callers to evaluate
snpcallers = ["lofreq", "varscan", "clc", "bcftools", "freebayes", "gatk"]#, "mutect2"]


mixed_strain_sample_ref = [sample for sample in sample_ref if not sample.split('.')[0].endswith(('-1-0', '-0-1'))]


sample_ref_tp = ['{}.{}'.format(sample, sample_refname_dict[sample]) for sample in sample_list if
                 sample not in ["TM-1-0", "TM-0-1", "TA-1-0", "TA-0-1"]]



# Selected caller for mutation context analysis
selected_snpcaller = "lofreq"

ruleorder: snp_evaluate > mutationcontext


# The final output of SNP calling
rule all:
    input:
        # snpcaller_performance_summary = results_dir + \
        #     "/final_tables/snpcaller_performance_summary.txt",
        expand(snp_dir + "/rtg/{snpcaller}/{sample_ref}.xindel/summary.txt", 
            sample_ref=mixed_strain_sample_ref, 
            snpcaller=snpcallers),
        expand(snp_dir + "/rtg/{snpcaller}/{sample_ref}.xsnp/summary.txt",
            sample_ref=mixed_strain_sample_ref, 
            snpcaller=snpcallers),

        performance_table = results_dir + \
            "/final_tables/caller_performance.tsv",
        rtg_roc_table = results_dir + "/final_tables/caller_rtg_roc_performance.tsv",
        rtg_snp_roc_table = results_dir + "/final_tables/caller_rtg_snp_roc_performance.tsv",

        rtg_roc_figure = results_dir + "/final_figures/caller_rtg_recall_precision_curve.pdf",
        rtg_snp_roc_figure = results_dir + "/final_figures/caller_rtg_snp_recall_precision_curve.pdf",
        # fp_tp_pdf = results_dir + "/final_figures/rtg_fp_tp_curve.pdf",
        
        
        performance_figure = results_dir + \
            "/final_figures/caller_performance.pdf",
        rtg_performance_figure = results_dir + \
            "/final_figures/caller_rtg_bestf1snp_score20all_score20snp_performance.pdf",

        fp_compare_figure = results_dir + "/final_figures/snpcaller_fp_snp_compare.pdf",


# Build the BWA and fai index for reference
include: "rules/index.smk"

# If not run on reads, copy the resulting VCF provided within the software for benchmarking
if not run_on_reads:
    rule cp_vcf:
        input: 
            expand(cd + "/data/snp/vcf/{{snpcaller}}/{{sample}}.{{ref}}.{{snpcaller}}.{ext}", 
                    ext=["vcf", "vcf.gz", "vcf.gz.tbi"])

        output: 
            expand(snpcall_dir + "/{{snpcaller}}/{{sample}}.{{ref}}.{{snpcaller}}.{ext}", 
                    ext=["vcf", "vcf.gz", "vcf.gz.tbi"])
            
            # vcf_bgz = snpcall_dir + "/{snpcaller}/{sample}.{ref}.{snpcaller}.vcf.gz",
        params:
            target_dir = snpcall_dir + "/{snpcaller}/"
        shell:
            """
            cp {input} {params.target_dir}
            """

    rule cp_genome_diff:
        input: 
            cd + "/data/snp/nucmer/{mix}.maskrepeat.variants.vcf",
            cd + "/data/snp/nucmer/{mix}.maskrepeat.variants.vcf.gz",
            cd + "/data/snp/nucmer/{mix}.maskrepeat.variants.vcf.gz.tbi"
        output: 
            snp_dir + "/nucmer/{mix}.maskrepeat.variants.vcf",
            snp_dir + "/nucmer/{mix}.maskrepeat.variants.vcf.gz",
            snp_dir + "/nucmer/{mix}.maskrepeat.variants.vcf.gz.tbi"
        params:
            target_dir = snp_dir + "/nucmer/"
        shell:
            """
            cp {input} {params.target_dir}
            """

else:
    rule cp_clc:
        input: 
            cd + "/data/snp/vcf/clc/{sample}.{ref}.clc.vcf"
            # cd + "/data/snp/clc/{sample}.{ref}.clc.vcf.gz",
            # cd + "/data/snp/clc/{sample}.{ref}.clc.vcf.gz.tbi"
        output: 
            snpcall_dir + "/clc/{sample}.{ref}.clc.vcf",
            # vcf_bgz = snpcall_dir + "/clc/{sample}.{ref}.clc.vcf.gz",
        params:
            target_dir = snpcall_dir + "/clc/"
        shell:
            """
            cp {input}* {params.target_dir}
            """

    # Remove remaining Phix, host read
    include: "rules/decontamination.smk"

    # BWA alignment
    include: "rules/bwa.smk"

    # Remove duplicate from BAM file using picard
    include: "rules/rmdup.smk"

    # Perform SNP calling
    include: "rules/vcfcall.smk"

    # Compare genome differences using NUCmer, minimap2
    include: "rules/genome_diff.smk"

# Extract the TP, FP SNPs
include: "rules/extract_TP.smk"

# Evaluate the SNP calling results
include: "rules/vis_eval_vcf.smk"

# Mutation context analysis for called SNPs
include: "rules/mutationcontext.smk"

# Compare the FP SNps
include: "rules/compare_FP.smk"

onsuccess:
    print("The SNPs calling evaluation is done!")
#    shell("mail -s 'The SNPs calling evaluation is done' youremail@provider.com")
