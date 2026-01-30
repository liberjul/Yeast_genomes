import sys
import pandas as pd
from Bio.Nexus import Nexus
from Bio import SeqIO
from find_types_class import exclude_chars_and_update_partitions

non_missing_char_dict = {}
regions = ["5_8S", "SSU", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
for reg in regions:
    print(reg)
    a = Nexus.Nexus(F"../data/phylo/aligned_seqs/class.yeast_seqs.{reg}.dedup_types_plus_new.nex")
nexi = [[reg, Nexus.Nexus(F"../data/phylo/aligned_seqs/class.yeast_seqs.{reg}.dedup_types_plus_new.nex")] for reg in regions]
# taxa_set = set()
# seqs_to_add = {k : {} for k in regions}
# for i in nexi: # get all taxon names
    # taxa_set.update(set(i[1].taxlabels))
# for reg in regions: # get sequences of new records
    # records_dict = SeqIO.to_dict(SeqIO.parse(F"../data/phylo/aligned_seqs/class.yeast_seqs.{reg}.dedup_types_plus_new.fasta", "fasta"))
    # new_recs = [x.strip("'") for x in records_dict if x.strip("'") not in taxa_set]
    # for rec in new_recs:
        # seqs_to_add[reg][rec] = str(records_dict[rec].seq).upper()
        # print(rec)
# for i in nexi: #add sequences to the alignments
    # for rec in seqs_to_add[i[0]]:
        # i[1].add_sequence(name=rec, sequence=seqs_to_add[i[0]][rec])
partition_buffer = "#nexus\nbegin sets;\n"
partition_count = 1
for i in nexi:
    print(F"Total characters before exclusion, {i[0]}: ", i[1].nchar)
    trimmed = exclude_chars_and_update_partitions(i[1]) # removed excluded chars
    print(trimmed)
    print(i[0])
    file = F"../data/phylo/aligned_seqs/class.yeast_seqs.{i[0]}_plus_new.excl.nex"
    i[1].write_nexus_data(filename=open(file, "w"), exclude = i[1].gaponly())
    matrix = i[1].crop_matrix(exclude = i[1].gaponly())
    for taxon in matrix:
        if taxon not in non_missing_char_dict:
            non_missing_char_dict[taxon] = {i[0] : len(matrix[taxon]) - matrix[taxon].count("?")}
        else:
            non_missing_char_dict[taxon][i[0]] = len(matrix[taxon]) - matrix[taxon].count("?")
    print(i[1].nchar)
    print(F"Total characters after exclusion, {i[0]}: ", i[1].nchar)
    if trimmed != {}:
        with open(file, "r") as ifile:
            lines = ifile.readlines()
            # charpartition UNTITLED = 1: 1-202\3 207-432\3, 2: 2-203\3 205-433\3, 3: 3-204\3 206-434\3;
            codonpart = [x for x in lines if "charpartition" in x][0].split(";")[0].split(" = ")[1]
            for i in codonpart.split(", "):
                chars = i.split(": ")[1]
                partition_buffer += F"\tcharset part{partition_count} = {file}: {chars};\n"
                partition_count += 1
    else:
        partition_buffer += F"\tcharset part{partition_count} = {file}: 1-{i[1].nchar - len(i[1].gaponly())};\n"
        partition_count += 1

partition_buffer += "end;\n"

with open(F"../data/phylo/aligned_seqs/partition.nex", "w") as ofile:
    ofile.write(partition_buffer)
print(partition_buffer)
df = pd.DataFrame.from_dict(non_missing_char_dict, orient="index")
df.to_csv("../data/phylo/info/non-missing_characters_per_taxon.txt", sep = "\t")
