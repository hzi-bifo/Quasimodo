# -*- coding: utf-8 -*-
# @Author: ZL Deng
# @Date:   2018-08-7 11:54:40

import argparse
from argparse import RawTextHelpFormatter

# Convert the normal VarScan output to VCF


def print_header():
    print(r'''##fileformat=VCFv4.1
##source=NUCmer
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth of quality bases">
##INFO=<ID=REF1,Number=1,Type=String,Description="The name of the 1st reference sequence">
##INFO=<ID=REF2,Number=1,Type=String,Description="The name of the 2nd reference sequence">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO''')


def snp2vcf(nucmer_snp_file):
    print_header()
    # header in SNP file: [P1]        [SUB]   [SUB]   [P2]    [BUFF]  [DIST]
    # [LEN R] [LEN Q] [FRM]   [TAGS]
    with open(nucmer_snp_file, "r") as snp_fh:
        # snp_fh.readline() # no header line
        for line in snp_fh:

            line = line.strip()
            line_entries = line.split("\t")
            p1, ref, alt, ref1_name, ref2_name = [
                line_entries[i] for i in [0, 1, 2, 10, 11]]
            if ref in list("ATCG") and alt in list("ATCG"):
                sid = "."
                filtered = "PASS"
                dp = "30"
                qual = "30"
                info = "DP={};REF1={};REF2={}".format(dp, ref1_name, ref2_name)
                out_line = "\t".join(
                    [ref1_name, p1, sid, ref, alt, qual, filtered, info])
                print(out_line)


if __name__ == "__main__":
    usage = r'''
    snp2vcf.py --- convert the normal SNP output from NUCmer to standard VCF format. 


    Usage:
    python snp2vcf.py <NUCmer SNP output file>
    
    This program output the converted format content to the stdout. You can redirect it into a file
    '''

    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=RawTextHelpFormatter,
                                     epilog="   by ZL Deng")

    parser.add_argument("snpfile", type=str,
                        help="the output file of NUCmer")
    args = parser.parse_args()
    snp2vcf(args.snpfile)
