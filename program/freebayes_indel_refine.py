#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
File: freebayes_refine.py
Created Date: March 1st 2020
Author: ZL Deng <dawnmsg(at)gmail.com>
---------------------------------------
Last Modified: 1st March 2020 5:46:08 pm
'''

import sys
import click
from os.path import abspath


@click.command()
@click.argument('vcf', type=click.Path(exists=True))
@click.option('-o', '--output', type=str, help='The output refined VCF file for FreeBayes')
def refine_vcf(vcf, output):
    out_fh = open(abspath(output), 'w') if output else sys.stdout
    with open(abspath(vcf), 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                out_fh.write(line)
            else:
                cols = line.split('\t')
                ref = cols[3].strip()
                len_ref = len(ref)
                alt = cols[4].strip()
                len_alt = len(alt)
                if len_ref == len_alt == 1:
                    out_fh.write(line)
                elif len(ref) > len(alt):
                    for i in range(1, len_alt + 1):
                        if ref[-i] == alt[-i]:
                            continue
                        else:
                            ref_idx = len_ref - i + 1
                            alt_idx = len_alt - i + 1
                            if ref[:alt_idx] == alt[:alt_idx]:
                                out_line = '\t'.join([cols[0], cols[1], str(int(
                                    cols[1]) + alt_idx - 1), cols[2], ref[alt_idx - 1:-i + 1], alt[alt_idx - 1]]) + '\n'  # + cols[5:])
                                # cols[1] = str(int(cols[1]) + alt_idx - 1)
                                # cols[3] = ref[alt_idx - 1:-i + 1]
                                # cols[4] = alt[alt_idx - 1]
                                # out_fh.write('\t'.join(cols))
                                out_fh.write(out_line)
                            break
                else:
                    for i in range(1, len_ref + 1):
                        if ref[-i] == alt[-i]:
                            continue
                        else:
                            ref_idx = len_ref - i + 1
                            alt_idx = len_alt - i + 1
                            if ref[:ref_idx] == alt[:ref_idx]:
                                out_line = '\t'.join([cols[0], cols[1], str(int(
                                    cols[1]) + ref_idx - 1), cols[2], ref[ref_idx - 1], alt[ref_idx - 1:-i + 1]]) + '\n'  # + cols[5:])
                                # cols[1] = str(int(cols[1]) + ref_idx - 1)
                                # cols[3] = ref[ref_idx - 1]
                                # cols[4] = alt[ref_idx - 1:-i + 1]
                                # out_fh.write('\t'.join(cols))
                                out_fh.write(out_line)
                            break
    out_fh.close()


if __name__ == "__main__":
    refine_vcf()
