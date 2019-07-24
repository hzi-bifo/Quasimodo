#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: ZL Deng
# @Date:   2018-08-7 11:54:40

import argparse
from argparse import RawTextHelpFormatter

## Convert the normal VarScan output to VCF

def print_header():
    print("##fileformat=VCFv4.1\n"
          "##source=VarScan2\n"
          "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth of quality bases\">\n"
          "##INFO=<ID=PV,Number=1,Type=Float,Description=\"Significance of variant read count vs. expected baseline error\">\n"
          "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n"
          "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n"
          "##FILTER=<ID=str10,Description=\"Depth over 8, PV below 0.001, QUAL over 20, MQUAL over 20\">\n"
          "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO")

def varscan2vcf(varscan_out_file):
    print_header()
    with open(varscan_out_file, "r") as varscan_fh:
        varscan_fh.readline()
        for line in varscan_fh:

            line = line.strip()
            line_entries = line.split("\t")
            chrom, pos, ref, alt, qual1, qual2, pvalue, af, reff, refr, altf, altr = [line_entries[i] for i in [0,1,2,18,9,10,11,6,14,15,16,17]]
            if len(ref) == 1 and len(alt) == 1:
                sid = "."
                filtered = "PASS"

                dp = str(sum(map(int,[reff, refr, altf, altr])))
                qual = str((int(qual1)+int(qual2))/2)
                af = str(float(af.strip("%"))/100)
                info = "DP={};PV={};AF={};DP4={},{},{},{}".format(dp, pvalue, af, reff, refr, altf, altr) 
                out_line = "\t".join([chrom, pos, sid, ref, alt, qual, filtered, info])
                print(out_line)


if __name__ == "__main__":
    usage = r'''
    varscan2vcf.py --- convert the normal VarScan output to standard VCF format. 
        This program requires Python modules argparse.

    Usage:
    python varscan2vcf.py <VarScan output file> >
    
    This program output the converted format content to the stdin. You can redirect it into a file
    '''

    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=RawTextHelpFormatter,
                                     epilog="   by ZL Deng")

    parser.add_argument("varscanfile", type=str,
                        help="the output file of VarScan")
    args = parser.parse_args()
    varscan2vcf(args.varscanfile)
