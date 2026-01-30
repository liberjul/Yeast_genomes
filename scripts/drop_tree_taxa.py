#!/usr/bin/env python
import glob, argparse, os
from Bio.Nexus import Nexus

'''
python drop_tree_taxa.py \
    -t "Atractiella_rhizophila" \
    -i ../data/phylo/aligned_seqs/20250904_class.yeast_seqs. \
    -s _plus_new.excl.ready.nex \
    -n taxa_dropout_003 \
    -p ../data/phylo/aligned_seqs/20250904_partition_no_introns.nex
'''

def remove_taxa(file, taxa):
    nex_file = Nexus.Nexus(file)
    print("Before trim: ", nex_file.ntax, nex_file.nchar)
    nex_file.matrix = nex_file.crop_matrix(delete = taxa)
    nex_file.taxlabels = [x for x in nex_file.taxlabels if x not in taxa]
    nex_file.ntax = len(nex_file.taxlabels)
    gaps = nex_file.gaponly()
    nex_file.matrix = nex_file.crop_matrix(exclude = gaps)
    gaps = set(gaps)
    nex_file.nchar = len(nex_file.matrix[nex_file.taxlabels[0]])
    for i in nex_file.charpartitions["UNTITLED"]:
        nex_file.charpartitions["UNTITLED"][i] = [x for x in nex_file.charpartitions["UNTITLED"][i] if x not in gaps]
        
    print("After trim: ", nex_file.ntax, nex_file.nchar)
    return nex_file

#    with open(file, "r") as ifile:
#        lines = 

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="prefix for input nexus")
parser.add_argument("-s", "--suff", type=str, help="suffix for input nexus")
parser.add_argument("-t", "--taxa", type=str, help="taxa to drop, separated by semi-colons")
parser.add_argument("-n", "--name", type=str, help="prefix for output alignments and partition")
parser.add_argument("-p", "--partition", type=str, help="partition file in nex format")
args = parser.parse_args()

taxa = args.taxa.split(";")

files = glob.glob(args.input + "*" + args.suff)

spl_in = args.input.split("/")
outdir = "/".join(spl_in[:-1])

print(files)
partition_count = 1
partition_buffer = "#nexus\nbegin sets;\n"
for i in files:
    fn = F"{outdir}/{args.name}.{os.path.basename(i)}" # add prefix before filename
    print(fn)
    nex_file = Nexus.Nexus(i)
    taxa_to_drop = [x for x in taxa if x in nex_file.taxlabels]
    nex_file.write_nexus_data(filename = open(fn, "w"), delete = taxa_to_drop) # remove taxa
    nex_file = Nexus.Nexus(fn) # reimport for correct formatting
    nex_file.write_nexus_data(filename = open(fn, "w"), exclude = nex_file.gaponly()) # remove gap-only characters
    with open(fn, "r") as ifile: # read in for partition information
        lines = ifile.readlines()
        # charpartition UNTITLED = 1: 1-202\3 207-432\3, 2: 2-203\3 205-433\3, 3: 3-204\3 206-434\3;
        part_lines = [x for x in lines if "charpartition" in x]
        if len(part_lines) == 0:
            partition_buffer += F"\tcharset part{partition_count} = {fn}: 1-{nex_file.nchar - len(nex_file.gaponly())};\n"
            partition_count += 1
        else:
            codonpart = part_lines[0].split(";")[0].split(" = ")[1]
            for i in codonpart.split(", "):
                chars = i.split(": ")[1]
                partition_buffer += F"\tcharset part{partition_count} = {fn}: {chars};\n"
                partition_count += 1
partition_buffer += "end;\n"

print(partition_buffer)

with open(F"{outdir}/{args.name}.{os.path.basename(args.partition)}", "w") as ofile:
    ofile.write(partition_buffer)