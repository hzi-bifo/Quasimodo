import pandas as pd
import os

# Load the configuration
configfile: "config/config.yaml"

# References for Merlin, TB40E, AD169
merlin_ref = os.path.abspath(config["MerlinRef"])
tb_ref = os.path.abspath(config["TB40ERef"])
ad_ref = os.path.abspath(config["AD169Ref"])
phix_ref = os.path.abspath(config["PhixRef"])

run_on_reads = config["runOnReads"]

cd = os.path.dirname(os.path.abspath(__name__))
# Define the folder to put the outputs
project_dir = config["outpath"].rstrip("/")


sample_refname_dict = {"TM-0-1": "Merlin", "TM-1-1": "Merlin",
                       "TM-1-10": "Merlin", "TM-1-50": "Merlin", "TM-1-0": "TB40E",
                       "TA-1-0": "TB40E", "TA-1-1": "AD169", "TA-1-10": "AD169",
                       "TA-1-50": "AD169", "TA-0-1": "AD169"}



# Sample list
if not run_on_reads:
    # sample_list, ref_names = glob_wildcards(os.path.join(cd,
    #                                                      "data/snp/vcf/clc") + "/{sample}.{ref_name}.clc.vcf")

    # Extract variants
    if not os.path.exists(cd + "/data/snp"):
        shell("tar -xzvf {} -C {}".format(cd + "/data/snp.tar.gz",
                                          cd + "/data/"))
    
    sample_list, ref_names = glob_wildcards(os.path.join(cd,
                                                         "data/snp/vcf/clc") + "/{sample}.{ref_name}.clc.vcf")

else:
    samples = pd.read_csv(config["samplesDesc"],
                          dtype=str, sep="\t").set_index(["sample"], drop=False)
    sample_list = sorted(list(samples["sample"]))

sample_ref = ['{}.{}'.format(sample, sample_refname_dict[sample])
              for sample in sample_list]


# Reports, QC sequences, results, assemblies, metaquast dir
report_dir = "/".join([project_dir, "reports"])
seq_dir = "/".join([project_dir, "data/seqs"])
results_dir = "/".join([project_dir, "results"])

# Number of threads to use
threads = config["threads"]

# Wildcards constrains
wildcard_constraints:
    sample = "[^\.\/]+",
    ref_name = "[^\.\/]+"

# Load read pairs

rm_human_ecoli = config["rmHumanEcoli"]
if rm_human_ecoli:
    human_ref_idx = config["HumanRefBWAIdx"]

    ecoli_ref_idx = config["EcoliRefBWAIdx"]

def get_fastq(w):
    return samples.loc[(w.sample), ["r1", "r2"]].dropna()


# Get the correct reference for each mixture sample
def ref_seq(w):
    if sample_refname_dict[w.sample] == "Merlin":
        ref = merlin_ref
    elif sample_refname_dict[w.sample] == "TB40E":
        ref = tb_ref
    else:
        ref = ad_ref
    return (ref, ref + ".bwt", ref + ".sdf")