#!/usr/bin/env python

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Fasta file")
parser.add_argument("-o", "--out", type=str, default = "", help="Corrected output file")
args = parser.parse_args()


with open(args.input, "r") as ifile:
    line = ifile.readline()
    line = ifile.readline()
    buffer = ""
    feat = ""
    while line != "":
        spl = line.strip().split("\t")
        lt, contig, start, stop, type, prod = spl[1:]
        if contig != feat:
            buffer += F">Feature {contig}\n"
            feat = contig
        buffer += F"{start}\t{stop}\tgene\n"
        buffer += F"\t\t\tlocus_tag\t{lt}\n"
        buffer += F"{start}\t{stop}\t{type}\n"
        if "internal transcribed spacer" in prod:
            buffer += F"\t\t\tncRNA_class\tother\n"
        buffer += F"\t\t\tproduct\t{prod}\n"
        line = ifile.readline()
        
with open(args.out, "w") as ofile:
    ofile.write(buffer)