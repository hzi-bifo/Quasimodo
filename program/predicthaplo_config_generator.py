#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
File: predicthaplo_config_generator.py
Created Date: February 28th 2020
Author: ZL Deng <dawnmsg(at)gmail.com>
---------------------------------------
Last Modified: 28th February 2020 10:39:28 pm
'''

import click
from os import path


@click.command()
@click.argument("prefix", type=str)
@click.argument("ref", type=click.Path(exists=True))
@click.argument("sam", type=click.Path(exists=True))
@click.option("--haplotypes", type=click.Choice(("0", "1")), default="0",
              help='True haplotype available? (0: no, 1: yes)')
@click.option('-m', '--msa', type=str, default='haplotypes.msa.fasta',
              help='The MSA of true haplotypes in fasta format. (fill in any dummy filename)')
@click.option('--start', type=int, default=1, help='reconstruction start')
@click.option('--stop', type=int, default=None, help='reconstruction stop')
@click.option('--readlen', type=int, default=100, help='minimum read length')
def config_generator(prefix, ref, sam, haplotypes, msa, start, stop, readlen):
    if stop is None:
        from Bio import SeqIO
        total_len = 0
        for record in SeqIO.parse(ref, "fasta"):
            record_len = len(record.seq)
            if record_len >= 1000:
                total_len += record_len
        stop = total_len

    (ref, sam, msa) = map(path.abspath, [ref, sam, msa])

    out_config = '''% configuration file for the HIVhaplotyper
% prefix
{prefix}_
% filename of reference sequence (FASTA)
{ref}
% do_visualize (1 = true, 0 = false)
1
% filname of the aligned reads (sam format)
{sam}
% have_true_haplotypes  (1 = true, 0 = false)
{haplotypes}
% filname of the true haplotypes (MSA in FASTA format) (fill in any dummy filename if there is no "true" haplotypes)
{msa}
% do_local_analysis  (1 = true, 0 = false) (must be 1 in the first run)
1
% max_reads_in_window;
10000
% entropy_threshold
4e-2
%reconstruction_start
{start}
%reconstruction_stop
{stop}
%min_mapping_qual
30
%min_readlength
{readlen}
%max_gap_fraction (relative to alignment length)
0.05
%min_align_score_fraction (relative to read length)
0.35
%alpha_MN_local (prior parameter for multinomial tables over the nucleotides)
25
%min_overlap_factor (reads must have an overlap with the local reconstruction window of at least this factor times the window size)
0.85
%local_window_size_factor (size of  local reconstruction window relative to the median of the read lengths)
0.7
% max number of clusters (in the truncated Dirichlet process)
25
% MCMC iterations
501
% include deletions (0 = no, 1 = yes)
1'''.format(prefix=prefix, ref=ref, sam=sam, haplotypes=haplotypes, msa=msa, start=start, stop=stop, readlen=readlen)

    print(out_config)


if __name__ == '__main__':
    config_generator()
