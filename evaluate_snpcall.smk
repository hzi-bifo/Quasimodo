include: "rules/load_config.smk"

snp_dir = "/".join([results_dir, "snp"])
snpcall_dir = "/".join([snp_dir, "callers"])

# Two mixtures
genome_diff_list = ["TM", "TA"]

# SNP callers to evaluate
snpcallers = ["lofreq", "varscan", "clc", "bcftools", "freebayes"]


# Reference for SNP calling of each sample
sample_refname_dict = {"TM-0-1": "Merlin", "TM-1-1": "Merlin",
                       "TM-1-10": "Merlin", "TM-1-50": "Merlin", "TM-1-0": "TB40E",
                       "TA-1-0": "TB40E", "TA-1-1": "AD169", "TA-1-10": "AD169",
                       "TA-1-50": "AD169", "TA-0-1": "AD169"}

sample_ref = ['{}.{}'.format(sample, sample_refname_dict[sample])
              for sample in sample_list]
sample_ref_tp = ['{}.{}'.format(sample, sample_refname_dict[sample]) for sample in sample_list if
                 sample not in ["TM-1-0", "TM-0-1", "TA-1-0", "TA-0-1"]]


# Get the correct reference for each mixture sample
def ref_seq(w):
    if sample_refname_dict[w.sample] == "Merlin":
        ref = merlin_ref
    elif sample_refname_dict[w.sample] == "TB40E":
        ref = tb_ref
    else:
        ref = ad_ref
    return [ref, ref + ".bwt"]

# Get the dir where the reference file is


def get_ref_dir(w):
    if sample_refname_dict[w.sample] == "Merlin":
        ref = merlin_ref
    elif sample_refname_dict[w.sample] == "TB40E":
        ref = tb_ref
    else:
        ref = ad_ref
    return os.path.dirname(ref) + "/"


# Selected caller for mutation context analysis
selected_snpcaller = "lofreq"

ruleorder: snp_evaluate > mutationcontext


def get_genome_diff(wc):
    if wc.sample.startswith("TM"):
        genome_diff = snp_dir + "/nucmer/TM.maskrepeat.snps"
    else:
        genome_diff = snp_dir + "/nucmer/TA.maskrepeat.snps"
    return genome_diff


# The final output of SNP calling
rule all:
    input:
        snpcaller_performance_summary = results_dir + \
            "/final_tables/snpcaller_performance_summary.txt",
        snpcaller_performance_figure = results_dir + \
            "/final_figures/snpcaller_performance.pdf",
        vcf = expand(snpcall_dir + "/{snpcallers}/{sample_ref}.{snpcallers}.vcf",
                     sample_ref=sample_ref, snpcallers=snpcallers),
        snps = expand(
            snp_dir + "/nucmer/{genome_diff}.maskrepeat.snps", genome_diff=genome_diff_list),
        fp = expand(snpcall_dir + "/{snpcallers}/fp/{sample_ref}.{snpcallers}.fp.vcf",
                    snpcall_dir=snpcall_dir, snpcallers=snpcallers, sample_ref=sample_ref),
        snp_profile_figure = expand(results_dir + "/final_figures/{mix}.{selected_snpcaller}.snp.profile.pdf",
                                    mix=["TM", "TA"], selected_snpcaller=selected_snpcaller),
        mutationcontext_figure = expand(results_dir + "/final_figures/{mix}.{selected_snpcaller}.mutationcontext.total.tp.fp.pdf",
                                        mix=["TM", "TA"], selected_snpcaller=selected_snpcaller),
        averaged_mc_figure = expand(results_dir + "/final_figures/{mix}.{selected_snpcaller}.mutationcontext.total.tp.fp.averaged.pdf",
                                    mix=["TM", "TA"], selected_snpcaller=selected_snpcaller),
        fp_compare_figure = results_dir + "/final_figures/snpcaller_fp_snp_compare.pdf"

rule cp_clc:
    input: "data/clc/{sample}.{ref}.clc.vcf"
    output: snpcall_dir + "/clc/{sample}.{ref}.clc.vcf"
    shell:
        """
        cp {input} {output}
        """

# Build the BWA and fai index for reference
include: "rules/index.smk"

# Remove remaining Phix reads
include: "rules/rm_phix.smk"

# BWA alignment
include: "rules/bwa.smk"

# Remove duplicate from BAM file using picard
include: "rules/rmdup.smk"

# Perform SNP calling
include: "rules/snpcall.smk"

# Compare genome differences using NUCmer
include: "rules/nucmmer.smk"

# Extract the TP, FP SNPs
include: "rules/extract_TP.smk"

# Evaluate the SNP calling results
include: "rules/snpEvaluate.smk"

# Mutation context analysis for called SNPs
include: "rules/mutationcontext.smk"

# Compare the FP SNps
include: "rules/compare_FP.smk"

onsuccess:
    print("The SNPs calling evaluation is done!")
#    shell("mail -s 'The SNPs calling evaluation is done' youremail@provider.com")
