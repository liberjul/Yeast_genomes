import glob, json

fastas = glob.glob("../unaligned_seqs/li_etal_accs_*.fasta")

seq_dict = {}

for i in fastas:
    with open(i, "r") as ifile:
        seq = ""
        lines = ifile.readlines()
        for line in lines:
            if line != "" and line[0] == ">":
                if seq != "":
                    seq_dict[header_line] = seq
                header_line = line
                seq = ""
            elif line != "" and line != "\n":
                seq += line
        if seq != "":
            seq_dict[header_line] = seq

buf_ssu = ""
buf_its = ""
buf_lsu = ""
buf_rpb1 = ""
buf_rpb2 = ""
buf_tef1 = ""
buf_cytb = ""
metadata_dict = {}
for header in seq_dict:
    if "UNVERIFIED" in header:
        acc, foo, gen, sp, first = header[1:].split(" ")[0:5]
    else:
        acc, gen, sp, first = header[1:].split(" ")[0:4]
    if "small subunit ribosomal" in header:
        buf_ssu = F"{buf_ssu}{header}{seq_dict[header]}"
        if "18S" in header:
            strain = header.split("18S")[0].split(sp)[1]
        elif "small subunit" in header:
            strain = header.split("small subunit")[0].split(sp)[1]
    if "internal transcribed" in header:
        buf_its = F"{buf_its}{header}{seq_dict[header]}"
        if "18S" in header:
            if first == "18S":
                strain = ""
            else:
                strain = header.split("18S")[0].split(sp)[1]
        elif "small subunit" in header:
            if first == "small":
                strain = ""
            else:
                strain = header.split("small")[0].split(sp)[1]
        else:
            if first == "internal":
                strain = ""
            else:
                strain = header.split("internal")[0].split(sp)[1]
    if "large subunit ribosomal" in header:
        buf_lsu = F"{buf_lsu}{header}{seq_dict[header]}"
        if "18S" in header:
            if first == "18S":
                strain = ""
            else:
                strain = header.split("18S")[0].split(sp)[1]
        elif "small subunit" in header:
            if first == "small":
                strain = ""
            else:
                strain = header.split("small")[0].split(sp)[1]
        elif "internal" in header:
            if first == "internal":
                strain = ""
            else:
                strain = header.split("internal")[0].split(sp)[1]
        elif "26S" in header:
            if first == "26S":
                strain = ""
            else:
                strain = header.split("26S")[0].split(sp)[1]
        elif "large subunit" in header:
            if first == "large":
                strain = ""
            else:
                strain = header.split("large subunit")[0].split(sp)[1]
    if "RNA polymerase II largest" in header:
        buf_rpb1 = F"{buf_rpb1}{header}{seq_dict[header]}"
        if first == "RNA":
            strain = ""
        else:
            strain = header.split("RNA")[0].split(sp)[1]
    if "RNA polymerase II second" in header:
        buf_rpb2 = F"{buf_rpb2}{header}{seq_dict[header]}"
        if first == "RNA":
            strain = ""
        else:
            strain = header.split("RNA")[0].split(sp)[1]
    if "translation elongation factor" in header:
        buf_tef1 = F"{buf_tef1}{header}{seq_dict[header]}"
        if first == "translation":
            strain = ""
        else:
            strain = header.split("translation elongation factor")[0].split(sp)[1]
    if "cytochrome b" in header:
        buf_cytb = F"{buf_cytb}{header}{seq_dict[header]}"
        if first == "cytochrome":
            strain = ""
        else:
            strain = header.split("cytochrome")[0].split(sp)[1]
    if strain != " ":
        metadata_dict[acc] = [gen, sp, strain.replace(":", "")]
    else:
        metadata_dict[acc] = [gen, sp, ""]


with open("../unaligned_seqs/ssu.fasta", "w") as ofile:
    ofile.write(buf_ssu)
with open("../unaligned_seqs/its1_58s_its2.fasta", "w") as ofile:
    ofile.write(buf_its)
with open("../unaligned_seqs/lsu.fasta", "w") as ofile:
    ofile.write(buf_lsu)
with open("../unaligned_seqs/rpb1.fasta", "w") as ofile:
    ofile.write(buf_rpb1)
with open("../unaligned_seqs/rpb2.fasta", "w") as ofile:
    ofile.write(buf_rpb2)
with open("../unaligned_seqs/tef1.fasta", "w") as ofile:
    ofile.write(buf_tef1)
with open("../unaligned_seqs/cytb.fasta", "w") as ofile:
    ofile.write(buf_cytb)

with open("../info/accession_metadata.json", "w") as ofile:
    json.dump(metadata_dict, ofile)
