#!/usr/bin/env python

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Fasta file")
parser.add_argument("-o", "--out", type=str, default = "", help="Corrected output file")
parser.add_argument("-m", "--contam", default= "", type=str, help="Contamination to remove, in the form of <contig1>:start..stop;<contig2>:start..stop (1-indexed).")
parser.add_argument("-p", "--print", action="store_true", help="Only print contig lengths.")
args = parser.parse_args()

record_dict = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

if args.contam != "":
    spl = args.contam.split(";")
    # Format: {contig: [start, stop]}
    contam_dict = {x.split(":")[0] : x.split(":")[1].split("..") for x in spl}

if args.print:
    for i in record_dict:
        print(F"{i}:{len(record_dict[i])}")
else:
    for i in record_dict:
        if i in contam_dict:
            start, stop = int(contam_dict[i][0]), int(contam_dict[i][1]) # these are 1-indexed
            print(i)
            print(record_dict[i])
            print("Length before trim: ", len(record_dict[i]))
            if start == 1:
                record_dict[i] = record_dict[i][stop:]
            elif stop == len(record_dict[i]):
                record_dict[i] = record_dict[i][:start-1]
            else:
                print("Removing middle region, beware of index shifts")
                record_dict[i] = record_dict[i][:start-1] + record_dict[i][stop:]
            print(record_dict[i])
            print("Length after trim: ", len(record_dict[i]))

    with open(args.out, "w") as ofile:
        for i in record_dict:
            ofile.write(F">{i}\n{str(record_dict[i].seq)}\n")