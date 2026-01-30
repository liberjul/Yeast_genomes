#!/usr/bin/env python

import primer3, argparse, json, sys, subprocess, os

'''
cd /mnt/c/Users/jal138/Duke\ Bio_Ea\ Dropbox/He\ Lab/Julian_Liber/Yeast_genomes
python ./scripts/correct_gene_names_tbl.py \
  -t ./data/annotation/JL201_funannotate_annotate_18march25/annotate_results/Aimania_JL201.tbl \
  -l OIV83 \
  -o ./data/annotation/JL201_funannotate_annotate_18march25/annotate_corrected/Aimania_JL201.tbl

python ./scripts/correct_gene_names_tbl.py \
  -t ./data/annotation/JL221_funannotate_annotate_18march25/annotate_results/Aimania_JL221.tbl \
  -l OIO90 \
  -o ./data/annotation/JL221_funannotate_annotate_18march25/annotate_corrected/Aimania_JL221.tbl

python ./scripts/correct_gene_names_tbl.py \
  -t ./data/annotation/NB124_funannotate_annotate_18march25/annotate_results/Aimania_NB124.tbl \
  -l ACM66B \
  -o ./data/annotation/NB124_funannotate_annotate_18march25/annotate_corrected/Aimania_NB124.tbl

'''

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tbl", type=str, help="Feature table file to correct")
parser.add_argument("-l", "--locus_tag", type=str, help="Locus tag prefixes to use")
parser.add_argument("-o", "--out", type=str, help="Corrected output file")
parser.add_argument("-c", "--contigs", default= "", type=str, help="Tab-delimited file to change contig names from first col to second.")
# parser.add_argument("-m", "--contam", default= "", type=str, help="Contamination to remove, in the form of <contig1>:start..stop;<contig2>:start..stop (1-indexed).")
args = parser.parse_args()

feat_hier_dict = {}
cds_dict = {}
mrna_dict = {}
text_blocks = {}
loc_dict = {}
with open(args.tbl, "r") as ifile:
    line = ifile.readline()
    feature_name = ""
    while line != "":
        if ">Feature" in line:
            feature_name = line.strip().replace(">Feature ", "")
            line = ifile.readline()
        else:
            spl = line.strip().split("\t")
            if len(spl) > 2:
                # print(spl)
                if spl[2] == "gene":
                    block_text = line
                    start, stop = spl[:2]
                    # print(spl)
                    while spl[0] != "locus_tag":
                        line = ifile.readline()
                        block_text += line
                        # print(line)
                        spl = line.strip().split("\t")
                    # print(spl)
                    id = spl[1]
                    feat_hier_dict[id] = {"feature" : feature_name,
                                          "loc" : F"{start}-{stop}",
                                          "transcript_id" : set(),
                                          "protein_id" : set()} # sets because we only need unique values
                    text_blocks[id] = block_text
                    line = ifile.readline()
                    spl = line.strip().split("\t")
                elif spl[2] == "mRNA":
                    block_text, feat_id, loc = "", "", ""
                    # print(line[-4:])
                    while line != "" and line[0] != ">" and line[-4:] != "CDS\n":
                        block_text += line
                        # print(line)
                        if line[0] != "\t":
                            loc += F"{spl[0]}-{spl[1]},"
                        if spl[0] == "transcript_id":
                            feat_id = spl[1]
                        line = ifile.readline()
                        spl = line.strip().split("\t")
                    text_blocks[feat_id] = block_text
                    loc_dict[feat_id] = loc
                    feat_hier_dict[id]["transcript_id"].add(feat_id)
                elif spl[2] == "CDS":
                    block_text, feat_id = "", ""
                    while line != "" and line[0] != ">" and (line[-5:] != "mRNA\n" and line[-5:] != "gene\n"):
                        block_text += line
                        # print(line)
                        if line[0] != "\t":
                            loc += F"{spl[0]}-{spl[1]},"
                        if spl[0] == "protein_id":
                            feat_id = spl[1]
                        line = ifile.readline()
                        spl = line.strip().split("\t")
                    text_blocks[feat_id] = block_text
                    loc_dict[feat_id] = loc
                    feat_hier_dict[id]["protein_id"].add(feat_id)
                else:
                    line = ifile.readline()
            else:
                # print(spl)
                line = ifile.readline()
            

gene_ids = feat_hier_dict.keys()
old_new_genes_dict = {}
old_new_genes_dict = {x: F"{args.locus_tag}_{i+1:06d}" for i, x in enumerate(gene_ids)}


old_new_trans_dict = {}
old_new_prot_dict = {}
feat_id_list = []
feature = ""
for x in feat_hier_dict:
    # print(x, feat_hier_dict[x])
    if feat_hier_dict[x]["feature"] != feature:
        feature = feat_hier_dict[x]["feature"]
        feat_id_list.append(F">Feature {feature}\n")
    feat_id_list.append(x)
    locus_tag_new = old_new_genes_dict[x]
    proteins = list(feat_hier_dict[x]["protein_id"])
    proteins.sort()
    # print(proteins)
    excluded_prots = []
    for i in proteins:
        for j in proteins:
            if i != j and (loc_dict[i] == loc_dict[j]):
                pair = [i, j]
                pair.sort()
                excluded_prots.append(pair[1])
    included_prots = [p for p in proteins if p not in excluded_prots]
    # print(excluded_prots)
    for i, y in enumerate(included_prots):
        old_new_trans_dict[F"{y}_mrna"] = F"gnl|ncbi|{locus_tag_new}-T{i+1}_mrna"        
        old_new_prot_dict[y] = F"gnl|ncbi|{locus_tag_new}-T{i+1}"
        feat_id_list.extend([F"{y}_mrna", y])
        
if args.contigs != "":
    with open(args.contigs, "r") as ifile:
        lines = ifile.readlines()
        contig_dict = {x.split()[0] : x.strip().split()[1] for x in lines}

# print(old_new_trans_dict)

text_blocks_new = {}
for i in text_blocks:
    lines = text_blocks[i].split("\n")
    lines = [F"\t\t\tlocus_tag\t{old_new_genes_dict[x.strip().split()[1]]}" if "locus_tag" in x else x for x in lines]
    lines = [F"\t\t\tprotein_id\t{old_new_prot_dict[x.strip().split()[1]]}" if "protein_id" in x else x for x in lines]
    lines = [F"\t\t\ttranscript_id\t{old_new_trans_dict[x.strip().split()[1]]}" if "transcript_id" in x else x for x in lines]
    text_blocks_new[i] = "\n".join(lines)

# print(text_blocks_new)
# print(old_new_genes_dict["gene-1882"])
with open(args.out, "w") as ofile:
    buffer = ""
    for i in feat_id_list:
        if i[0] == ">":
            if buffer != "":
                ofile.write(buffer)
                buffer = ""
            if args.contigs != "":
                new_contig = contig_dict[i.strip().split(" ")[1]]
                buffer += F">Feature {new_contig}\n"
            else:
                buffer += line
        else:
            buffer += text_blocks_new[i]
    ofile.write(buffer)




# buffer = ""
# with open(args.tbl, "r") as ifile, open(args.out, "w") as ofile:
    # line = ifile.readline()
    # buffer = ""
    # while line != "":
        # spl = line.strip().split("\t")
        # if spl[0] == "locus_tag":
            # buffer += F"\t\t\tlocus_tag\t{old_new_genes_dict[spl[1]]}\n"
        # elif spl[0] ==  "transcript_id":
            # buffer += F"\t\t\ttranscript_id\t{old_new_trans_dict[spl[1]]}\n"
        # elif spl[0] ==  "protein_id":
            # buffer += F"\t\t\tprotein_id\t{old_new_prot_dict[spl[1]]}\n"
        # elif line[:8] == ">Feature":
            # if args.contigs != "":
                # new_contig = contig_dict[line.strip().split(" ")[1]]
                # buffer += F">Feature {new_contig}\n"
            # else:
                # buffer += line
            # ofile.write(buffer)
            # buffer = ""
        # else:
            # buffer += line
        # line = ifile.readline()
    # ofile.write(buffer)
        
        