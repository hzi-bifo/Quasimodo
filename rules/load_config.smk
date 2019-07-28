import pandas as pd
import os

# Load the configuration
configfile: "config/config.yaml"

# This rule requires the R 3.5.1 conda environment

# Sample list
samples = pd.read_csv(config["samplesDesc"],
                      dtype=str, sep="\t").set_index(["sample"], drop=False)
#sample_list = list(samples.index)
sample_list = sorted(list(samples["sample"]))

# References for Merlin, TB40E, AD169
merlin_ref = os.path.abspath(config["MerlinRef"])
tb_ref = os.path.abspath(config["TB40ERef"])
ad_ref = os.path.abspath(config["AD169Ref"])
phix_ref = os.path.abspath(config["PhixRef"])

run_on_reads = config["runOnReads"]

cd = os.path.dirname(os.path.abspath(__name__))
# Define the folder to put the outputs
project_dir = config["outpath"].rstrip("/")

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


def get_fastq(w):
    return samples.loc[(w.sample), ["r1", "r2"]].dropna()
