import glob, json

with open("../info/accession_metadata.json", "r") as ifile:
    metadata_dict = json.dump(ifile)

strain_set = set()
for acc in metadata_dict:
    strain_set.add(" ".join(metadata_dict[acc]))

strain_list = list(strain_set)

regions = ["SSU", "ITS1", "5_8S", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
for reg in regions:
    with open(F"../aligned_seqs/yeast_seqs.{reg}.fasta", "r") as ifile:
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
        seq_dict_strain[metadata_dict[acc]] = [acc, seq_dict[acc]]
    buffer = ""
    for i in strain_list:
        if i in seq_dict_strain:
            buffer = F"{buffer}>{seq_dict_strain[acc][0]}|{i}|{reg}\n{seq_dict_strain[acc][1]}\n"
        else:
            buffer = F"{buffer}>{seq_dict_strain[acc][0]}|{i}|{reg}\n{'-'*max_length}\n"
    with open(F"../aligned_seqs/yeast_seqs.{reg}.ordered.fasta", "w") as ofile:
        ofile.write(buffer)
