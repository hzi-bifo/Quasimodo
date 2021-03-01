#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# This script was modified based on the script in this repo:
# https://github.com/MatteoSchiavinato/Utilities
# ----------------------------------------------------------------------------

import argparse
import sys
from time import strftime
from Bio import SeqIO
from operator import itemgetter

if len(sys.argv) <= 1:
    sys.argv.append("-h")

if sys.argv[1] in ["-h", "--help", "getopt", "usage", "-help", "help"]:
    sys.exit('''

USAGE:  mummer2vcf.py [ options ]

Convert Mummer SNP/indel output as produced by the 'show-snps -T' command to
a VCF format.

OPTIONS:

    -s|--snps		    	Mummer "snps" file (output of 'show-snps -T')
				(mandatory)

       --input-header		Use this flag if the file has a header

    -n|--no-Ns		    	Exclude variants featuring Ns
				(default: off)

    -t|--type   SNP|INDEL|ALL   Restrict the output to just SNPs or indels
				(default: ALL)

       --output-header		Add a newly-generated VCF header to the output 
				(Requires also --reference)

    -g|--reference		Reference genome FASTA file for header generation

''')


# parser
p = argparse.ArgumentParser()
p.add_argument("-s", "--snps", required=True)
p.add_argument("--input-header", action="store_true")
p.add_argument("-n", "--no-Ns", action="store_true")
p.add_argument("-t", "--type", choices=["SNP", "INDEL", "ALL"], default="ALL")
p.add_argument("--output-header", action="store_true")
p.add_argument("-g", "--reference")
args = p.parse_args()


### conditions ###

if (args.output_header) and not args.reference:
    sys.exit("ERROR: --add-vcf-header requires --reference as well\n\n")

# functions


def remove_header(Lines):
    Lines_no_header = Lines[4:]
    return Lines_no_header


def convert_snps_to_vcf(line):
    #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    lst = line.rstrip("\n\r\b").split("\t")
    scaffold = lst[10]
    pos = lst[0]
    ID = "."
    ref = lst[1]
    alt = lst[2]
    qual = "30"
    filt = "PASS"
    # DP=30;REF1=Merlin_1555;REF2=TB40E-GFP_MM_ann
    info = "DP=30;REF1={};REF2={}".format(scaffold, lst[11])
    var_type = "."
    orig_pos = str(lst[11]) + ":" + str(lst[3])
    vcf_line = "\t".join([scaffold, pos, ID, ref, alt,
                          qual, filt, info, var_type, orig_pos])
    return vcf_line


def exclude_Ns(Vcf_lines):
    Vcf_lines_no_Ns = []
    for line in Vcf_lines:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")

        if (("N" not in ref) and ("n" not in ref) and
                ("N" not in alt) and ("n" not in alt)):

            Vcf_lines_no_Ns.append(
                "\t".join((scaffold, pos, ID, ref, alt, qual, filt, info, var_type, orig_pos)))

    return Vcf_lines_no_Ns


def infer_var_type(Vcf_lines):
    Vcf_lines_w_types = []
    for line in Vcf_lines:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")

        if ((len(ref) == 1) and ("." not in ref) and
                (len(alt) == 1) and ("." not in alt)):
            var_type = "."

        else:
            var_type = "INDEL"

        Vcf_lines_w_types.append(
            "\t".join((scaffold, pos, ID, ref, alt, qual, filt, info, var_type, orig_pos)))

    return Vcf_lines_w_types


def collapse_snps(Sorted_snps):
    Collapsed_snps = []
    for line in Sorted_snps:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")

        if len(Collapsed_snps) == 0:
            Collapsed_snps.append(
                [scaffold, pos, ID, ref, alt, qual, filt, info, var_type, orig_pos])

        elif int(pos) == pos_old:
            alt_lst = str(Collapsed_snps[-1][4]).split(",")
            if str(alt) not in alt_lst:
                alt_lst.append(str(alt))
                Collapsed_snps[-1][4] = ",".join(alt_lst)

        else:
            Collapsed_snps.append(
                [scaffold, pos, ID, ref, alt, qual, filt, info, var_type, orig_pos])

        pos_old = int(pos)

    return Collapsed_snps


def collapse_indels(Indels):
    Collapsed_indels = []
    for line in Indels:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")

        if len(Collapsed_indels) == 0:
            Collapsed_indels.append(
                [scaffold, pos, ID, ref, alt, qual, filt, info, var_type, orig_pos])

        elif ((int(pos) == pos_old) and (ref == ".")) or \
             ((int(pos) == pos_old + 1) and (alt == ".")):

            # insertions
            if ref == ".":
                Tmp_insertions = Collapsed_indels[-1][4].split(",")

                if orig_pos != orig_pos_old:
                    x = ",".join([str(indel) + str(alt)
                                  for indel in Tmp_insertions])
                    Collapsed_indels[-1][4] = x

                elif orig_pos == orig_pos_old:
                    ### NEW INDELS ###
                    indel = Collapsed_indels[-1][4]
                    Collapsed_indels[-1][4] = ",".join([indel,
                                                        str(indel[0:len(indel) - 1]) + str(alt)])
            # deletions
            elif alt == ".":
                Collapsed_indels[-1][3] = str(Collapsed_indels[-1]
                                              [3]) + str(ref)

        else:
            Collapsed_indels.append(
                [scaffold, pos, ID, ref, alt, qual, filt, info, var_type, orig_pos])

        pos_old = int(pos)
        orig_pos_old = orig_pos

    return Collapsed_indels


def add_indel_first_base(Collapsed_indels, Sequences):
    Indels_w_first_base = []
    Joined_indels = ["\t".join(x) + "\n" for x in Collapsed_indels]
    for line in Joined_indels:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")
        # double -1, one is to make it 0-based, one is to take the position before the listed one
        ref_position = str(Sequences[scaffold][int(pos) - 1 - 1])
        if ref == ".":
            ref = ref_position
            alt = ",".join([str(ref_position) + str(indel)
                            for indel in alt.split(",")])
        elif alt == ".":
            ref = ref_position + ref
            alt = ref_position

        # scaling pos within the file to the position before according to the guidelines
        pos = str(int(pos) - 1)
        Indels_w_first_base.append(
            [scaffold, pos, ID, ref, alt, qual, filt, info, var_type, orig_pos])

    return Indels_w_first_base


def parse_sequences(Reference):
    Sequences = {}
    for record in SeqIO.parse(Reference, "fasta"):
        name = str(record.id)
        seq = str(record.seq)
        Sequences[name] = seq

    return Sequences


def separate_vars_by_type(Vcf_lines):
    Indels = []
    Snps = []
    for line in Vcf_lines:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")
        if var_type == ".":
            Snps.append(line)
        elif var_type == "INDEL":
            Indels.append(line)

    return (Indels, Snps)


def collapse_variants(Reference, Vcf_lines):
    Sequences = parse_sequences(Reference)
    (Indels, Snps) = separate_vars_by_type(Vcf_lines)

    ### snps ###
    Sorted_snps = []
    for line in Snps:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")
        pos = int(pos)
        Sorted_snps.append((scaffold, pos, ID, ref, alt,
                            qual, filt, info, var_type, orig_pos))

    Sorted_snps = ["\t".join([str(y) for y in x])
                   for x in sorted(Sorted_snps, key=itemgetter(1))]
    Collapsed_snps = collapse_snps(Sorted_snps)

    ### indels ###
    Indel_tuples = []
    for line in Indels:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")
        pos = int(pos)
        Indel_tuples.append((scaffold, pos, ID, ref, alt,
                             qual, filt, info, var_type, orig_pos))

    Indels = ["\t".join([str(y) for y in x]) for x in Indel_tuples]
    Collapsed_indels = collapse_indels(Indels)
    ### add prior base to indels ###
    Collapsed_indels = add_indel_first_base(Collapsed_indels, Sequences)

    Vcf_lines = []
    for x in Collapsed_snps:
        Vcf_lines.append(x)
    for y in Collapsed_indels:
        Vcf_lines.append(y)

    ### sort ###
    Vcf_lines = sorted(Vcf_lines, key=itemgetter(1))

    ### return ###
    Vcf_lines = ["\t".join([str(y) for y in x]) for x in Vcf_lines]
    return Vcf_lines


def select_by_type(Vcf_lines, selected_type):
    Selected_vcf_lines = []
    for line in Vcf_lines:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = line.rstrip("\n\r\b").split("\t")

        if ((selected_type == "ALL") or
            (selected_type == "SNP" and var_type == ".") or
                (selected_type == "INDEL" and var_type == "INDEL")):

            Selected_vcf_lines.append(line)

    return Selected_vcf_lines


def sort_vcf_lines(Vcf_lines):
    Lines = [x.rstrip("\n\r\b").split("\t") for x in Vcf_lines]
    Lines_w_int = []
    for lst in Lines:
        lst[1] = int(lst[1])
        Lines_w_int.append(lst)

    Sorted_lines = sorted(Lines_w_int, key=itemgetter(0, 1))
    Vcf_lines = []
    for lst in Sorted_lines:
        (scaffold, pos, ID, ref, alt, qual, filt, info, var_type,
         orig_pos) = lst
        var_type = 'INDEL' if var_type == 'INDEL' else 'SNV'
        info = '{};ORIG={};TYPE={}'.format(info, orig_pos, var_type)

        pos = str(pos)
        out_lst = [scaffold, pos, ID, ref, alt, qual, filt, info]
        line = "\t".join(out_lst)
        Vcf_lines.append(line)

    return Vcf_lines


def add_header(Reference, Vcf_lines):
    file_format = "##fileformat=VCFv4.2"
    file_date = "##fileDate={0}".format(strftime("%Y%m%d"))
    source = "##source=mummer2vcf.py"
    reference_file = "##reference={0}".format(str(Reference))
    Header = [file_format, file_date, source, reference_file]

    Contigs_in_vars = list(
        set([line.rstrip("\n\r\b").split("\t")[0] for line in Vcf_lines]))
    Contigs = []
    for record in SeqIO.parse(Reference, "fasta"):
        ID = str(record.id)
        if ID in Contigs_in_vars:
            Length = str(len(str(record.seq)))
            Contigs.append("##contig=<ID={0},length={1}>".format(ID, Length))

    if len(Contigs) > 0:
        contig_lines = "\n".join(Contigs)
        Header.append(contig_lines)

    dp_line = '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth of quality bases">'
    ref1_line = '##INFO=<ID=REF1,Number=1,Type=String,Description="The name of the 1st reference sequence">'
    ref2_line = '##INFO=<ID=REF2,Number=1,Type=String,Description="The name of the 2nd reference sequence">'
    orig_line = '##INFO=<ID=ORIG,Number=1,Type=String,Description="The original position of variant at 2nd reference sequence">'

    # indel_line = '##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">'
    var_type_line = '##INFO=<ID=TYPE,Number=1,Type=String,Description="Indicates that the variant is an INDEL or SNV.">'
    Header.extend([dp_line, ref1_line, ref2_line, orig_line, var_type_line])

    final_line = "\t".join(
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    Header.append(final_line)
    Header = "\n".join(Header)
    Final_lines = [Header]
    for line in Vcf_lines:
        Final_lines.append(line)

    return Final_lines


# main script
if __name__ == "__main__":
    infile = args.snps
    INPUT = open(infile, "r")
    Lines = [line for line in INPUT]
    INPUT.close()

    if args.input_header:
        Lines = remove_header(Lines)

    Vcf_lines = [convert_snps_to_vcf(line) for line in Lines]

    if args.no_Ns:
        Vcf_lines = exclude_Ns(Vcf_lines)

    Vcf_lines = infer_var_type(Vcf_lines)
    Vcf_lines = collapse_variants(args.reference, Vcf_lines)
    Vcf_lines = select_by_type(Vcf_lines, args.type)
    Vcf_lines = sort_vcf_lines(Vcf_lines)

    if args.output_header:
        Vcf_lines = add_header(args.reference, Vcf_lines)

    for line in Vcf_lines:
        sys.stdout.write(line + "\n")
