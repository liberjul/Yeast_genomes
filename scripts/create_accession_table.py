#!/usr/bin/env python
import os, argparse, json, subprocess
import pandas as pd

'''
python create_accession_table.py \
    -i ../data/phylo/unaligned_seqs/class.yeast_seqs \
    -m ../data/phylo/info/accession_metadata.json \
    -n ../data/phylo/info/manual.accession_metadata.json \
    -r ../data/phylo/info/refseq_nucleotide_pairs.json \
    -s ../data/phylo/info/cleaned_syns_dict.json \
    -o ../data/phylo/info/
'''

def clean_strain(strain):
    strain = strain.replace(" ", "").rstrip("T")
    for hit in ["strain", "culture", "isolate", "voucher"]:
        strain_list = strain.split(hit)
        if hit in strain and strain_list[0] != "":
            strain = strain_list[0]
            break
        else:
            strain = strain.replace(hit, "")
    return strain
def strain_set_sort(strain_list):
    strain_list = list(strain_list)
    out_list = []
    for i in ["AS", "ATCC", "CBS", "CGMCC", "JCM", "NRRL", "PYCC", "DSM", "AFTOL"]:
        out_list.extend([x for x in strain_list if i in x])
        strain_list = [x for x in strain_list if i not in x]
    return("/".join(out_list+strain_list))



parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="prefix for input FASTAs")
parser.add_argument("-m", "--metadata", type=str, help="input json with accession metadata")
parser.add_argument("-n", "--manual_metadata", type=str, help="input json with accession metadata, manually inputted")
parser.add_argument("-r", "--refseq", type=str, help="input json with refseq-nucleotide accession pairs", default="")
parser.add_argument("-s", "--synonym", type=str, help="input json with type strain synonyms of a binomial")
parser.add_argument("-o", "--output", type=str, help="output file path")
args = parser.parse_args()


if os.path.exists(args.metadata):
    with open(args.metadata, "r") as ifile:
        metadata_dict = json.load(ifile)
else:
    raise IOError(F"Metadata file {args.metadata} does not exist.")

if os.path.exists(args.manual_metadata):
    with open(args.manual_metadata, "r") as ifile:
        metadata_dict.update(json.load(ifile))
else:
    raise IOError(F"Metadata file {args.manual_metadata} does not exist.")

if args.refseq != "" and os.path.exists(args.refseq):
    use_refseq = True
    with open(args.refseq, "r") as ifile:
        refseq_dict = json.load(ifile)
else:
    use_refseq = False

with open(args.synonym, "r") as ifile:
    cleaned_syns_dict = json.load(ifile)
    for i in cleaned_syns_dict:
        cleaned_syns_dict[i] = set(cleaned_syns_dict[i])


for i in ["SSU", "ITS1", "5_8S", "ITS2", "LSU"]:
    with open(F"{args.input}.{i}.fasta", "r") as ifile:
        lines = ifile.readlines()
        accessions = [x[1:].split("|")[0] for x in lines if x[0] == ">"]
    for acc in accessions:
        metadata_dict[acc].append(i)
for i in ["RPB1", "RPB2", "TEF1", "CYTB"]:
        with open(F"{args.input}.{i}.fasta", "r") as ifile:
            lines = ifile.readlines()
            accessions = [x[1:].split(" ")[0] for x in lines if x[0] == ">"]
        for acc in accessions:
            metadata_dict[acc].append(i)
table_dict = {}
unique_refseq_accs = set()
no_strain_accs = ""
for i in metadata_dict:
    gen, sp, strain = metadata_dict[i][:3]
    binom = F"{gen} {sp}"
    strain = clean_strain(strain)
    print(i, gen, sp, strain)
    if strain == "":
        no_strain_accs += F"{i},"
    else:
        metadata_dict[i][2] = strain

subprocess.run(F'efetch -db nuccore -id "{no_strain_accs.strip(",")}" -api_key=95c4f8e7f8558cac9bbb733cb87ea362ff09 -format gb > temp.gb', shell=True)
with open("temp.gb", "r") as ifile:
    lines = ifile.readlines()
    lines = [x for x in lines if "strain=" in x or "VERSION" in x]
previous_is_acc = False
for line in lines:
    if "VERSION" in line and not previous_is_acc:
        acc = line.strip().split()[-1]
        previous_is_acc = True
    elif "strain=" in line:
        previous_is_acc = False
        metadata_dict[acc][2] = line.strip().split("=")[-1].strip('"').split(";")[0]
    else:
        strain = ""


for i in metadata_dict:
    if "NR_" in i:
        unique_refseq_accs.add(i)
    else:
        gen, sp, strain = metadata_dict[i][:3]
        if sp != "sp.":
            binom = F"{gen} {sp}"
            strain = clean_strain(strain)
            if strain in cleaned_syns_dict[binom]:
                strain = strain_set_sort(cleaned_syns_dict[binom])
            name = " ".join([gen, sp, strain])
            if name not in table_dict:
                table_dict[name] = {}
            loci = metadata_dict[i][3:]
            for locus in loci:
                if ":" not in i:
                    table_dict[name][locus] = i.split(".")[0]
                else:
                    table_dict[name][locus] = i

df = pd.DataFrame.from_dict(table_dict, orient="index")
df.to_csv(args.output + "loci_accessions.txt", sep = "\t")
with open(args.output + "refseq_accessions.txt", "w") as ofile:
    ofile.write("\n".join(list(unique_refseq_accs)))

print(no_strain_accs.strip(","))
# missing_ssu_or_lsu = ""
# for i in table_dict:
#     if ("ITS1" in table_dict[i] and "5_8S" in table_dict[i] and "ITS2" in table_dict[i]) and ("SSU" not in table_dict[i] or "LSU" not in table_dict[i]):
#         if len(set([table_dict[i]["ITS1"], table_dict[i]["5_8S"], table_dict[i]["ITS2"]])) == 1:
#             # print(i, table_dict[i]["ITS1"])
#             missing_ssu_or_lsu += table_dict[i]["ITS1"] + ","
# print(missing_ssu_or_lsu)
