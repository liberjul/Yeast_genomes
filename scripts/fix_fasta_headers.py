#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input FASTA path")
parser.add_argument("-o", "--output", type=str, help="Output FASTA path")
args = parser.parse_args()

with open(args.input, "r") as ifile:
    with open(args.output, "w") as ofile:
        line = ifile.readline()
        while line != "":
            if line[0] == ">":
                ofile.write(line.split(" ")[0] + "\n")
            else:
                ofile.write(line)
            line = ifile.readline()