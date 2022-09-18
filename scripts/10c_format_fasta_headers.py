import json, os

data_dir = "../data/phylo"
with open(F"{data_dir}/info/accession_metadata.json", "r") as ifile:
    metadata_dict = json.load(ifile)

if os.path.exists(F"{data_dir}/info/new.accession_metadata.json"):
    with open(F"{data_dir}/info/new.accession_metadata.json", "r") as ifile:
        metadata_dict.update(json.load(ifile))

strain_set = set()
for acc in metadata_dict:
    metadata_dict[acc][2] = metadata_dict[acc][2].strip(" ")
    strain_set.add(" ".join(metadata_dict[acc]))
for i in [0, 2, 4, 6, 8, 10, 11, 12, 13]:
    metadata_dict[F"NS0000{str(i).zfill(2)}.1"] = ["Novel", "sp.", "strain JL201"]
strain_set.add("Novel sp. strain JL201")
for i in [1, 3, 5, 7, 9, 14, 15, 16, 17]:
    metadata_dict[F"NS0000{str(i).zfill(2)}.1"] = ["Novel", "sp.", "strain JL221"]
strain_set.add("Novel sp. strain JL221")

with open(F"{data_dir}/info/updated.accession_metadata.json", "w") as ofile:
    json.dump(metadata_dict, ofile)

strain_list = list(strain_set)
print(strain_list)

regions = ["SSU", "ITS1", "5_8S", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
for reg in regions:
    with open(F"{data_dir}/aligned_seqs/yeast_seqs.{reg}.mafft.clipkit.fasta", "r") as ifile:
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
    with open(F"{data_dir}/aligned_seqs/yeast_seqs.{reg}.mafft.ordered.fasta", "w") as ofile:
        ofile.write(buffer)
