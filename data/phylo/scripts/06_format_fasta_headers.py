import json, os

with open("../info/accession_metadata.json", "r") as ifile:
    metadata_dict = json.load(ifile)

if os.path.exists("../info/new.accession_metadata.json"):
    with open("../info/new.accession_metadata.json", "r") as ifile:
        metadata_dict.update(json.load(ifile))

strain_set = set()
for acc in metadata_dict:
    metadata_dict[acc][2] = metadata_dict[acc][2].strip(" ")
    strain_set.add(" ".join(metadata_dict[acc]))
for i in range(10):
    metadata_dict[F"NS00000{i}.1"] = ["Novel", "sp.", "strain JL201"]
strain_set.add("Novel sp. strain JL201")

strain_list = list(strain_set)
print(strain_list)

regions = ["SSU", "ITS1", "5_8S", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
for reg in regions:
    with open(F"../aligned_seqs/yeast_seqs.{reg}.mafft.clipkit.fasta", "r") as ifile:
        seq_dict = {}
        min_length = 9999
        max_length = 0
        seq = ""
        lines = ifile.readlines()
        for line in lines:
            if line != "" and line[0] == ">":
                if seq != "":
                    seq_dict[acc] = seq
                    if len(seq) > max_length:
                        max_length = len(seq)
                    if len(seq) < min_length:
                        min_length = len(seq)
                header_line = line
                if "|F|" in header_line:
                    acc = header_line[1:].split("|")[0]
                else:
                    acc = header_line[1:].split(" ")[0]
                seq = ""
            elif line != "" and line != "\n":
                seq += line.strip()
        if seq != "":
            seq_dict[acc] = seq
            if len(seq) > max_length:
                max_length = len(seq)
            if len(seq) < min_length:
                min_length = len(seq)
    if max_length != min_length:
        raise ValueError(F"Aligned sequences for {reg} are not all equal length!")
    seq_dict_strain = {}
    for acc in seq_dict:
        seq_dict_strain[" ".join(metadata_dict[acc])] = [acc, seq_dict[acc]]
    buffer = ""
    for i in strain_list:
        if i in seq_dict_strain:
            buffer = F"{buffer}>{i}|{reg}\n{seq_dict_strain[i][1]}\n"
        else:
            buffer = F"{buffer}>{i}|{reg}\n{'-'*max_length}\n"
    with open(F"../aligned_seqs/yeast_seqs.{reg}.mafft.ordered.fasta", "w") as ofile:
        ofile.write(buffer)
