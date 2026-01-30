import glob, json, argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--new", action="store_true", help="process and print new sequences")
args = parser.parse_args()
data_dir = "../data/phylo/unaligned_seqs"
mbm_genera = (["Rosettozyma", "Trigonosporomyces", "Heterogastridium", "Spencerozyma", "Vonarxula", "Oberwinklerozyma", "Yamadamyces", "Kriegeria", "Phenoliferia",
              "Fellozyma", "Bannozyma", "Chrysozyma", "Hamamotoa", "Yurkovia", "Udeniozyma", "Colacogloea", "Glaciozyma", "Pseudohyphozyma", "Slooffia", "Leucosporidiella"
              "Reniforma", "Yunzhangia", "Ustilentyloma", "Microbotryozyma", "Microbotryum", "Sphacelotheca", "Sampaiozyma", "Leucosporidium", "Pseudoleucosporidium",
              "Curvibasidium", "Heitmania", "Rhodosporidiobolus", "Rhodotorula", "Sporobolomyces", "Camptobasidium", "Psychromyces", "Cryolevonia", "Meredithblackwellia",
              "Sporidiobolus", "Libkindia", "Crucella", "[Rhodotorula]", "Pycnopulvinus"])

fastas = glob.glob(F"{data_dir}/mbm.*.fasta")
if args.new and len(fastas) == 0:
    raise ValueError("Cannot run in --new mode because there are not already existing mbm.*.fasta files.")
old_seqs = {}
for i in fastas:
    with open(i, "r") as ifile:
        seq = ""
        lines = ifile.readlines()
        for line in lines:
            if line != "" and line[0] == ">":
                if seq != "":
                    old_seqs[header] = seq
                header = line[1:].split(" ")[0].split("|")[0]
                seq = ""
            elif line != "" and line != "\n":
                seq += line
        if seq != "":
            old_seqs[header] = seq


fastas = glob.glob(F"{data_dir}/*MBM_types*fasta") #Changed to find Microbotryomycete seqs
fastas += glob.glob(F"{data_dir}/yeast_seqs.*.fasta")
fastas += glob.glob(F"{data_dir}/missing_ssu_lsu_accs*.fasta")
fastas += glob.glob(F"{data_dir}/all_refseq_nuccore_accs.fasta")
fastas += glob.glob(F"{data_dir}/2022_11_13_accs.fasta")

seq_dict = {}

for i in fastas:
    with open(i, "r") as ifile:
        seq = ""
        lines = ifile.readlines()
        for line in lines:
            if line != "" and line[0] == ">":
                if seq != "":
                    if args.new:
                        header = header_line[1:].split(" ")[0].split("|")[0]
                        if header not in old_seqs:
                            seq_dict[header_line] = seq
                    else:
                        seq_dict[header_line] = seq
                header_line = line
                seq = ""
            elif line != "" and line != "\n":
                seq += line
        if seq != "":
            if args.new:
                header = header_line[1:].split(" ")[0].split("|")[0]
                if header not in old_seqs:
                    seq_dict[header_line] = seq
            else:
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
    strain = ""
    if "UNVERIFIED" in header:
        acc, foo, gen, sp, first = header[1:].split(" ")[0:5]
    else:
        acc, gen, sp, first = header[1:].split(" ")[0:4]
    if gen in mbm_genera:
        if "small subunit ribosomal" in header or "18S" in header:
            buf_ssu = F"{buf_ssu}{header}{seq_dict[header]}"
            header_new = header.replace("genes for ", "").replace("genomic DNA containing ", "")
            if "18S" in header:
                strain = header_new.split("18S")[0].split(sp)[1]
            elif "small subunit" in header:
                strain = header_new.split("small subunit")[0].split(sp)[1]
            elif "genes for" in header:
                strain = header.split("genes")[0].split(sp)[1]
        if "internal transcribed" in header or "ITS" in header:
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
            elif "genes for" in header:
                strain = header.split("genes")[0].split(sp)[1]
            else:
                if first == "internal" or first == "ITS":
                    strain = ""
                else:
                    if "ITS" in header:
                        strain = header.split("ITS")[0].split(sp)[1]
                    else:
                        strain = header.split("internal")[0].split(sp)[1]
        if "large subunit ribosomal" in header or "26S" in header or "28S" in header:
            buf_lsu = F"{buf_lsu}{header}{seq_dict[header]}"
            header_new = header.replace("genes for ", "").replace("genomic DNA containing ", "")
            if "18S" in header:
                if first == "18S":
                    strain = ""
                else:
                    strain = header_new.split("18S")[0].split(sp)[1]
            elif "small subunit" in header:
                if first == "small":
                    strain = ""
                else:
                    strain = header_new.split("small")[0].split(sp)[1]
            elif "genes for" in header:
                strain = header.split("genes")[0].split(sp)[1]
            elif "internal" in header:
                if first == "internal":
                    strain = ""
                else:
                    strain = header_new.split("internal")[0].split(sp)[1]
            elif "26S" in header:
                if first == "26S":
                    strain = ""
                else:
                    strain = header_new.split("26S")[0].split(sp)[1]
            elif "28S" in header:
                if first == "28S":
                    strain = ""
                else:
                    strain = header_new.split("28S")[0].split(sp)[1]
            elif "large subunit" in header:
                if first == "large":
                    strain = ""
                else:
                    strain = header_new.split("large subunit")[0].split(sp)[1]
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
        if strain != " " and strain != "":
            metadata_dict[acc] = [gen, sp, strain.replace(":", "")]
        else:
            metadata_dict[acc] = [gen, sp, ""]

if args.new:
    pref = "mbm.new." #Changed for MBM
else:
    pref = "mbm." #Changed for MBM
with open(F"{data_dir}/{pref}ssu.fasta", "w") as ofile:
    ofile.write(buf_ssu)
with open(F"{data_dir}/{pref}its1_58s_its2.fasta", "w") as ofile:
    ofile.write(buf_its)
with open(F"{data_dir}/{pref}lsu.fasta", "w") as ofile:
    ofile.write(buf_lsu)
with open(F"{data_dir}/{pref}rpb1.fasta", "w") as ofile:
    ofile.write(buf_rpb1)
with open(F"{data_dir}/{pref}rpb2.fasta", "w") as ofile:
    ofile.write(buf_rpb2)
with open(F"{data_dir}/{pref}tef1.fasta", "w") as ofile:
    ofile.write(buf_tef1)
with open(F"{data_dir}/{pref}cytb.fasta", "w") as ofile:
    ofile.write(buf_cytb)

with open(F"{data_dir}/../info/{pref}accession_metadata.json", "w") as ofile:
    json.dump(metadata_dict, ofile)
