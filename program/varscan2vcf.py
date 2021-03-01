#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: ZL Deng
# @Date:   2018-08-7 11:54:40

import sys
import click
# from pathlib import Path
from os.path import abspath

# Convert the normal VarScan output to VCF


def header():
    return("##fileformat=VCFv4.1\n"
           "##source=VarScan2\n"
           "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth of quality bases\">\n"
           "##INFO=<ID=PV,Number=1,Type=Float,Description=\"Significance of variant read count vs. expected baseline error\">\n"
           "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n"
           "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n"
           "##FILTER=<ID=str10,Description=\"Depth over 8, PV below 0.001, QUAL over 20, MQUAL over 20\">\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")


@click.command()
@click.option('-s', '--snp', type=click.Path(exists=True), help='The VarScan output for SNPs')
@click.option('-i', '--indel', type=click.Path(exists=True), help='The VarScan output for INDELs')
@click.option('-o', '--output', type=str, help='The output VCF file of SNPs or/and INDELs')
def to_vcf_file(snp, indel, output):
    if not snp and not indel:
        raise Exception(
            'You need to specify at least one of SNP or INDEL file')
    else:
        out_fh = open(abspath(output), 'w') if output else sys.stdout
        out_fh.write(header() + '\n')
        if snp:
            snp_list = varscan2vcf_list(snp)
            if indel:
                indel_list = varscan2vcf_list(indel)
                concat_vcf_list = sorted(
                    snp_list + indel_list, key=lambda x: int(x[1]))
                for cols in concat_vcf_list:
                    out_fh.write('\t'.join(cols) + '\n')
            else:
                for cols in snp_list:
                    out_fh.write('\t'.join(cols) + '\n')

        else:
            indel_list = varscan2vcf_list(indel)
            for cols in indel_list:
                out_fh.write('\t'.join(cols) + '\n')

        out_fh.close()


def varscan2vcf_list(varscan_file):
    vcf_list = []
    with open(abspath(varscan_file), 'r') as fh:
        fh.readline()
        for line in fh:
            line = line.strip()
            line_entries = line.split("\t")
            chrom, pos, ref, alt, qual1, qual2, pvalue, af, reff, refr, altf, altr = [
                line_entries[i] for i in [0, 1, 2, 18, 9, 10, 11, 6, 14, 15, 16, 17]]
            # if len(ref) == 1 and len(alt) == 1:
            sid = "."
            filtered = "PASS"
            if alt.startswith('-'):
                ref, alt = ref + alt[1:], ref

            elif alt.startswith('+'):
                alt = ref + alt[1:]
            else:
                pass

            dp = str(sum(map(int, [reff, refr, altf, altr])))
            qual = str((int(qual1) + int(qual2)) / 2)
            af = str(round(float(af.strip("%")) / 100, 5))
            info = "DP={};PV={};AF={};DP4={},{},{},{}".format(
                dp, pvalue, af, reff, refr, altf, altr)
            out_cols = (chrom, pos, sid, ref, alt, qual, filtered, info)
            vcf_list.append(out_cols)
    return vcf_list


if __name__ == "__main__":
    to_vcf_file()
