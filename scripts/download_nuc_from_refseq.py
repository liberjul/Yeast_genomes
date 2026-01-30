import json, subprocess

with open("../data/phylo/info/refseq_nucleotide_pairs.json") as ifile:
    refseq_dict = json.load(ifile)

nucleotide_set = set()
for i in refseq_dict:
    nucleotide_set.add(refseq_dict[i])

nucl_str = ",".join(list(nucleotide_set))
subprocess.run(F'efetch -db nuccore -id "{nucl_str}" -api_key=95c4f8e7f8558cac9bbb733cb87ea362ff09 -format fasta > ../data/phylo/unaligned_seqs/all_refseq_nuccore_accs.fasta', shell=True)
