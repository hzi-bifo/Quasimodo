#!/usr/bin/env python
# -- coding: utf-8 --
# @Author: ZL Deng
# @Date:   2018-08-7 11:54:40

import argparse
import subprocess
import os
from argparse import RawTextHelpFormatter


def extract_tp_fp_snp(vcf_file, snp_file):
    """ Extract the TP SNPs from the VCF file of caller
    @param vcf_file: The input VCF file for extracting the true positive (TP) and false positive (FP) SNPS.
    @param snp_file: The SNPs between two reference genomes due to sequence differences.
    @param tp_out: The output of TP SNPs identified by caller which are originated from genome differences. 
    @param fp_out: The output of FP SNPs identified by caller by false calling. 
    """
    dirname = os.path.dirname(vcf_file)
    fname_wo_ext = os.path.basename(vcf_file)[:-4]
    filtered_out = vcf_file[:-4] + ".filtered.vcf"
    fp_out = os.path.join(dirname, "fp", fname_wo_ext + ".fp.vcf")

    snp_by_caller = r'''awk -F"\t" '$4~/^[ACGT]$/&&$5~/^[ACGT]$/&&$6>=20' {}'''.format(
        vcf_file)

    filter_qual_cmd = r'''(grep -E "^#" {};{}) > {}'''.format(vcf_file,
                                                              snp_by_caller, filtered_out)
    # print(filter_qual_cmd)
    filtered = subprocess.Popen(
        filter_qual_cmd, shell=True, executable="/bin/bash")
    filtered.communicate()
    if os.path.basename(vcf_file).split(".")[0].endswith(("-1-0", "-0-1")):
        fp_cmd = r'''cp {} {}'''.format(filtered_out, fp_out)
        fp = subprocess.Popen(fp_cmd, shell=True, executable="/bin/bash")
        fp.communicate()

    else:
        if not os.path.exists(os.path.join(dirname, "tp")):
            os.makedirs(os.path.join(dirname, "tp"))
        tp_out = os.path.join(dirname, "tp", fname_wo_ext + ".tp.vcf")
        genome_snp = r'''awk -F"\t" '$2!="."&&$3!="."{{print $1, ".", $2, $3}}' OFS="\t" {}'''.format(
            snp_file)

        tp_cmd = r'''(grep -E "^#" {};fgrep -wf <({}) <({})) > {}'''.format(vcf_file,
                                                                            genome_snp, snp_by_caller, tp_out)
        fp_cmd = r'''(grep -E "^#" {};fgrep -wvf <({}) <({})) > {}'''.format(vcf_file,
                                                                             genome_snp, snp_by_caller, fp_out)

        tp = subprocess.Popen(tp_cmd, shell=True, executable="/bin/bash")
        fp = subprocess.Popen(fp_cmd, shell=True, executable="/bin/bash")
        fp.communicate()


def extract_tp_fp_custom_snp(vcf_file, snp_file, outdir):
    """ Extract the TP SNPs from the VCF file of caller
    @param vcf_file: The input VCF file for extracting the true positive (TP) and false positive (FP) SNPS.
    @param snp_file: The SNPs between two reference genomes due to sequence differences.
    @param tp_out: The output of TP SNPs identified by caller which are originated from genome differences. 
    @param fp_out: The output of FP SNPs identified by caller by false calling. 
    @param outdir: the output directory for filtered and TP, FP SNP VCF files
    """
    #dirname = os.path.dirname(vcf_file)
    fname_wo_ext = os.path.basename(vcf_file)[:-4]
    filtered_out = os.path.join(outdir, fname_wo_ext + ".filtered.vcf")
    fp_out = os.path.join(outdir, "fp", fname_wo_ext + ".fp.vcf")

    snp_by_caller = r'''awk -F"\t" '$4~/^[ACGT]$/&&$5~/^[ACGT]$/&&$6>=20' {}'''.format(
        vcf_file)

    filter_qual_cmd = r'''(grep -E "^#" {};{}) > {}'''.format(vcf_file,
                                                              snp_by_caller, filtered_out)
    # print(filter_qual_cmd)
    filtered = subprocess.Popen(
        filter_qual_cmd, shell=True, executable="/bin/bash")
    filtered.communicate()
    if os.path.basename(vcf_file).split(".")[0].endswith(("-1-0", "-0-1")):
        fp_cmd = r'''cp {} {}'''.format(filtered_out, fp_out)
        fp = subprocess.Popen(fp_cmd, shell=True, executable="/bin/bash")
        fp.communicate()

    else:
        if not os.path.exists(os.path.join(outdir, "tp")):
            os.makedirs(os.path.join(outdir, "tp"))
        tp_out = os.path.join(outdir, "tp", fname_wo_ext + ".tp.vcf")
        genome_snp = r'''awk -F"\t" '$2!="."&&$3!="."{{print $1, ".", $2, $3}}' OFS="\t" {}'''.format(
            snp_file)

        tp_cmd = r'''(grep -E "^#" {};fgrep -wf <({}) <({})) > {}'''.format(vcf_file,
                                                                            genome_snp, snp_by_caller, tp_out)
        fp_cmd = r'''(grep -E "^#" {};fgrep -wvf <({}) <({})) > {}'''.format(vcf_file,
                                                                             genome_snp, snp_by_caller, fp_out)

        tp = subprocess.Popen(tp_cmd, shell=True, executable="/bin/bash")
        fp = subprocess.Popen(fp_cmd, shell=True, executable="/bin/bash")
        fp.communicate()


if __name__ == "__main__":
    usage = r'''
    extract_TP_FP_SNPs.py --- extract the TP and FP SNPs from VCF input and output them as VCF files. 
        This program requires Python modules argparse.

    Usage:
    python extract_TP_FP_SNPs.py <VCF input file> <SNPs file (genome differences)> 
            <VCF output file for TP> <VCF output file for FP>
    
    This program output the converted format content to the stdin. You can redirect it into a file
    '''

    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=RawTextHelpFormatter,
                                     epilog="   by ZL Deng")

    parser.add_argument("vcffile", type=str,
                        help="the input VCF file")
    parser.add_argument("snpfile", type=str,
                        help="the input SNP file of genome differences generated by MUMmer")
    parser.add_argument("data", type=str,
                        choices=["hcmv", "custom"],
                        help="the source of input data")
    parser.add_argument("outdir", type=str,
                        help="the output dir")
    args = parser.parse_args()
    if args.data == "hcmv":
        extract_tp_fp_snp(args.vcffile, args.snpfile)
    else:
        extract_tp_fp_custom_snp(args.vcffile, args.snpfile, args.outdir)
