from Bio.Nexus import Nexus
import re

for i in ["RPB1", "RPB2", "TEF1", "CYTB", "5_8S", "SSU", "LSU"]:
    pref = F"../data/phylo/aligned_seqs/class.yeast_seqs.{i}"
    a = Nexus.Nexus()
    a.read(input=F"{pref}.dedup_types.nex")
    a.terminal_gap_to_missing()
    a.export_fasta(F"{pref}.dedup_types.fasta")
    with open(F"{pref}.dedup_types.fasta", "r") as ifile, open(F"{pref}.dedup_types_no_miss.fasta", "w") as ofile:
        text = ifile.read()
        ofile.write(text.replace("?", "-"))
        # ofile.write(re.sub("Y|R|S|K|M|W|B|D|H|V", "N", text.upper().replace("?", "-")))
