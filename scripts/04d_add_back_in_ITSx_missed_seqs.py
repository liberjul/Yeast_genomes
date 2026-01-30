import glob, json, os, sys
import pandas as pd

data_dir = "../data/phylo/unaligned_seqs"

if os.path.exists("../data/phylo/info/new_names.json"):
    with open("../data/phylo/info/new_names.json", "r") as ifile:
        new_name_dict = json.load(ifile)

mbm_genera = (["Rosettozyma", "Trigonosporomyces", "Heterogastridium", "Spencerozyma", "Vonarxula", "Oberwinklerozyma", "Yamadamyces", "Kriegeria", "Phenoliferia",
              "Fellozyma", "Bannozyma", "Chrysozyma", "Hamamotoa", "Yurkovia", "Udeniozyma", "Colacogloea", "Glaciozyma", "Pseudohyphozyma", "Slooffia",
              "Reniforma", "Yunzhangia", "Ustilentyloma", "Microbotryozyma", "Microbotryum", "Sphacelotheca", "Sampaiozyma", "Leucosporidium", "Pseudoleucosporidium",
              "Curvibasidium", "Heitmania", "Rhodosporidiobolus", "Rhodotorula", "Sporobolomyces", "Camptobasidium", "Psychromyces", "Cryolevonia", "Meredithblackwellia",
              "Sporidiobolus", "Libkindia", "Crucella", "[Rhodotorula]", "Pycnopulvinus"])

fastas = glob.glob(F"{data_dir}/mbm.*su.fasta")
fastas += glob.glob(F"{data_dir}/outgroup.*su.fasta")
all_seq_dict = {}
for i in fastas:
    if "ssu" in i:
        region = "SSU"
    else:
        region = "LSU"
    with open(i, "r") as ifile:
        seq = ""
        lines = ifile.readlines()
        for line in lines:
            if line != "" and line[0] == ">":
                if seq != "":
                    all_seq_dict[header_line] = [region, seq]
                header_line = line
                seq = ""
            elif line != "" and line != "\n":
                seq += line
        if seq != "":
            all_seq_dict[header_line] = [region, seq]

fastas = glob.glob(F"{data_dir}/mbm.yeast_seqs.*SU.fasta")
fastas += glob.glob(F"{data_dir}/outgroup.yeast_seqs.*SU.fasta")
included_seq_dict = {}

for i in fastas:
    if "SSU" in i:
        region = "SSU"
    else:
        region = "LSU"
    with open(i, "r") as ifile:
        seq = ""
        lines = ifile.readlines()
        for line in lines:
            if line != "" and line[0] == ">":
                if seq != "":
                    if header_line not in included_seq_dict:
                        included_seq_dict[header_line] = {region : seq}
                    else:
                        included_seq_dict[header_line][region] = seq
                header_line = line
                seq = ""
            elif line != "" and line != "\n":
                seq += line
        if seq != "":
            if header_line not in included_seq_dict:
                included_seq_dict[header_line] = {region : seq}
            else:
                included_seq_dict[header_line][region] = seq

ssu_headers = {}
lsu_headers = {}
for header in all_seq_dict:
    strain = ""
    acc, gen, sp, first = header[1:].split(" ")[0:4]
    if "small subunit ribosomal" in header or "18S" in header:
        header_new = header.replace("genes for ", "").replace("genomic DNA containing ", "")
        if "18S" in header:
            strain = header_new.split("18S")[0].split(sp)[1]
        elif "small subunit" in header:
            strain = header_new.split("small subunit")[0].split(sp)[1]
        if strain != " " and strain != "":
            if F"{gen} {sp}" in new_name_dict:
                header_str = " ".join([new_name_dict[F"{gen} {sp}"], strain.replace(":", "").strip(" ")])
            else:
                header_str = " ".join([gen, sp, strain.replace(":", "").strip(" ")])
            if header_str in ssu_headers:
                ssu_headers[header_str].append(acc)
            else:
                ssu_headers[header_str] = [acc]
        else:
            if F"{gen} {sp}" in new_name_dict:
                header_str = " ".join([new_name_dict[F"{gen} {sp}"], ""])
            else:
                header_str = " ".join([gen, sp, ""])
            if header_str in ssu_headers:
                ssu_headers[header_str].append(acc)
            else:
                ssu_headers[header_str] = [acc]
    if "internal transcribed" in header or "ITS region" in header:
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
            if first == "internal" or first == "ITS":
                strain = ""
            else:
                if "ITS" in header:
                    strain = header.split("ITS")[0].split(sp)[1]
                else:
                    strain = header.split("internal")[0].split(sp)[1]
    if "large subunit ribosomal" in header or "26S" in header or "28S" in header:
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
        if strain != " " and strain != "":
            if F"{gen} {sp}" in new_name_dict:
                header_str = " ".join([new_name_dict[F"{gen} {sp}"], strain.replace(":", "").strip(" ")])
            else:
                header_str = " ".join([gen, sp, strain.replace(":", "").strip(" ")])
            if header_str in lsu_headers:
                lsu_headers[header_str].append(acc)
            else:
                lsu_headers[header_str] = [acc]
        else:
            if F"{gen} {sp}" in new_name_dict:
                header_str = " ".join([new_name_dict[F"{gen} {sp}"], ""])
            else:
                header_str = " ".join([gen, sp, ""])
            if header_str in lsu_headers:
                lsu_headers[header_str].append(acc)
            else:
                lsu_headers[header_str] = [acc]

data_dir = "../data/phylo"
with open(F"{data_dir}/info/mbm.accession_metadata.json", "r") as ifile:
    metadata_dict = json.load(ifile)
with open(F"{data_dir}/info/outgroup.accession_metadata.json", "r") as ifile:
    metadata_dict.update(json.load(ifile))

with open(F"{data_dir}/info/class.updated.accession_metadata.json", "r") as ifile:
    metadata_dict.update(json.load(ifile))

if os.path.exists(F"{data_dir}/info/manual.accession_metadata.json"):
    with open(F"{data_dir}/info/manual.accession_metadata.json", "r") as ifile:
        metadata_dict.update(json.load(ifile))

with open(F"{data_dir}/info/accession_metadata.json", "w") as ofile:
    json.dump(metadata_dict, ofile)

if os.path.exists(F"{data_dir}/info/acc_strain_pairs.json"):
    with open(F"{data_dir}/info/acc_strain_pairs.json", "r") as ifile:
        acc_strain_pairs_dict = json.load(ifile)

synonym_dict = {}

with open(F"{data_dir}/info/type_strains.txt", "r") as ifile:
    for i in ifile.readlines():
        spl = i.strip().split(' = ')
        spl = [x.replace(" ", "") for x in spl]
        spl[0] = spl[0].strip("T")
        synonym_dict[spl[0]] = set(spl[1:])
        for j in range(1, len(spl)):
            # print(spl[0:j] + spl[j+1:])
            synonym_dict[spl[j]] = set(spl[0:j] + spl[j+1:])

gen_sp_dict = {}
for i in metadata_dict:
    gen_sp = " ".join(metadata_dict[i][0:2])
    if gen_sp in gen_sp_dict:
        gen_sp_dict[gen_sp].add(metadata_dict[i][2])
    else:
        gen_sp_dict[gen_sp] = set([metadata_dict[i][2]])

# print(gen_sp_dict)
species_w_mult_strains = []
cleaned_syns_dict = {}
for i in gen_sp_dict:
    strain_list = list(gen_sp_dict[i])
    strain_set = set()
    for j in strain_list:
        strain_str = j.split(" ")
        strain_str = [x for x in strain_str if x not in ["culture", "strain", "isolate", "gene", "for"]]
        strain_set.add("".join(strain_str))
    strain_set = [x for x in strain_set if x != ""]
    cleaned_syns_dict[i] = set(strain_set)
    if len(strain_set) > 1:
        same = False
        for j in strain_set[1:]:
            if strain_set[0] in synonym_dict and j in synonym_dict[strain_set[0]]:
                same = True
        if not same:
            species_w_mult_strains.append(i)

if os.path.exists(F"{data_dir}/info/cleaned_syns_dict.json"):
    with open(F"{data_dir}/info/cleaned_syns_dict.json", "r") as ifile:
        new_dict = json.load(ifile)
        for i in new_dict:
            new_dict[i] = set(new_dict[i])
        cleaned_syns_dict.update(new_dict)
with open(F"{data_dir}/info/cleaned_syns_dict.json", "w") as ofile:
    out_dict = {}
    for i in cleaned_syns_dict:
        out_dict[i] = list(cleaned_syns_dict[i])
    json.dump(out_dict, ofile)

buffer = ""
count_in = 0
ssu_name_to_accs = {}
for i in ssu_headers:
    name_str = " ".join(i.split(" ")[:2])
    strain_list = i.replace("specimen-voucher", "").replace("voucher ", "").split(" ")[2:]
    strain_list = [x.replace(":", "").replace(" ", "") for x in strain_list if "&" not in x and "," not in x]
    strain_list = [x for x in strain_list if x not in ["culture", "strain", "isolate", "gene", "for"]] # remove descriptors
    strain_str = "".join(strain_list)
    if name_str not in cleaned_syns_dict:
        print(name_str, "Missing")
    elif (strain_str == "" or strain_str not in cleaned_syns_dict[name_str]) and os.path.exists(F"{data_dir}/info/acc_strain_pairs.json") and name_str in cleaned_syns_dict and ssu_headers[i][0] in acc_strain_pairs_dict:
        strain_list = "/".join(acc_strain_pairs_dict[ssu_headers[i][0]])
        strain_list = strain_list.replace("specimen-voucher", "").replace("voucher ", "").split("/")
        strain_list = [x.replace(":", "").replace(" ", "").strip("T") for x in strain_list if "&" not in x and "," not in x]
        in_syn_dict = [x for x in strain_list if x in cleaned_syns_dict[name_str]]
        strain_str = "/".join(strain_list)
        if len(in_syn_dict) > 0:
            cleaned_syns_dict[name_str].update(set(strain_list))
            print(ssu_headers[i], name_str, strain_str, cleaned_syns_dict[name_str], True)
            count_in += 1
            ssu_name_to_accs[name_str] = ssu_name_to_accs[name_str] + ssu_headers[i] if name_str in ssu_name_to_accs else ssu_headers[i]
        else:
            print(ssu_headers[i], name_str, strain_str, cleaned_syns_dict[name_str], False)
    elif name_str in cleaned_syns_dict:
        print(ssu_headers[i], name_str, strain_str, cleaned_syns_dict[name_str], strain_str in cleaned_syns_dict[name_str])
        if strain_str in cleaned_syns_dict[name_str]:
            count_in += 1
            ssu_name_to_accs[name_str] = ssu_name_to_accs[name_str] + ssu_headers[i] if name_str in ssu_name_to_accs else ssu_headers[i]
        else:
            buffer += F"{ssu_headers[i][0]}\n"
    else:
        print(ssu_headers[i], name_str, strain_str, False)
print(F"Total headers: {len(ssu_headers)}\nIn synonym dict: {count_in}")

count_in = 0
lsu_name_to_accs = {}
for i in lsu_headers: # each accession/record for which we have an LSU seq
    name_str = " ".join(i.split(" ")[:2]) # binomial name
    strain_list = i.replace("specimen-voucher", "").replace("voucher ", "").split(" ")[2:] # Remove extra wording from strain names
    strain_list = [x.replace(":", "").replace(" ", "") for x in strain_list if "&" not in x and "," not in x] # Remove soaces, colons, and non-strain elements
    strain_list = [x for x in strain_list if x not in ["culture", "strain", "isolate", "gene", "for"]] # remove descriptors
    strain_str = "".join(strain_list) # a single strain name, without spaces
    # First check if 1) no strain name or strain not in synonym dictionary for binomial AND 2) if accession-strain data exists AND 3) binomial in synonym dictionary keys
    if name_str not in cleaned_syns_dict:
        print(name_str, "Missing")
    elif (strain_str == "" or strain_str not in cleaned_syns_dict[name_str]) and os.path.exists(F"{data_dir}/info/acc_strain_pairs.json") and name_str in cleaned_syns_dict and lsu_headers[i][0] in acc_strain_pairs_dict:
        if lsu_headers[i][0] not in acc_strain_pairs_dict: # No data for accession-strain matching
            buffer += F"{lsu_headers[i][0]}\n"
        else:
            strain_list = "/".join(acc_strain_pairs_dict[lsu_headers[i][0]]) # Separate strains with slashes for later split
            strain_list = strain_list.replace("specimen-voucher", "").replace("voucher ", "").split("/")
            strain_list = [x.replace(":", "").replace(" ", "").strip("T") for x in strain_list if "&" not in x and "," not in x]
            in_syn_dict = [x for x in strain_list if x in cleaned_syns_dict[name_str]] # Strain names which are known synonyms
            strain_str = "/".join(strain_list) # rejoin as a single strain word with multiple synonyms
            if len(in_syn_dict) > 0: # If at least one strain name is a known synonym
                cleaned_syns_dict[name_str].update(set(strain_list)) # Add other strain names
                print(lsu_headers[i], name_str, strain_str, cleaned_syns_dict[name_str], True)
                count_in += 1
                lsu_name_to_accs[name_str] = lsu_name_to_accs[name_str] + lsu_headers[i] if name_str in lsu_name_to_accs else lsu_headers[i]
            else: # None of the strains are known synonyms
                print(lsu_headers[i], name_str, strain_str, cleaned_syns_dict[name_str], False)
    elif name_str in cleaned_syns_dict: # If the binomial is still in the strain dictionary
        print(lsu_headers[i], name_str, strain_str, cleaned_syns_dict[name_str], strain_str in cleaned_syns_dict[name_str])
        if strain_str in cleaned_syns_dict[name_str]: # If the strain name is a known synonym
            count_in += 1
            lsu_name_to_accs[name_str] = lsu_name_to_accs[name_str] + lsu_headers[i] if name_str in lsu_name_to_accs else lsu_headers[i]
        else:
            buffer += F"{lsu_headers[i][0]}\n"
    else:
        print(lsu_headers[i], name_str, strain_str, False)
        buffer += F"{lsu_headers[i][0]}\n"

print(F"Total headers: {len(lsu_headers)}\nIn synonym dict: {count_in}")
print(F"SSU unique names: {len(ssu_name_to_accs.keys())}")
print(F"LSU unique names: {len(lsu_name_to_accs.keys())}")

with open(F"{data_dir}/info/accs_wo_matching_strain.txt", "w") as ofile:
    ofile.write(buffer)
with open(F"{data_dir}/info/cleaned_syns_dict.json", "w") as ofile:
    out_dict = {}
    for i in cleaned_syns_dict:
        out_dict[i] = list(cleaned_syns_dict[i])
    json.dump(out_dict, ofile)
# print(lsu_headers)
max_accs = 0
for i in ssu_name_to_accs:
    if max_accs < len(ssu_name_to_accs[i]):
        max_accs = len(ssu_name_to_accs[i])
print(F"Maximum number of accessions for SSU: {max_accs}")
max_accs = 0
for i in lsu_name_to_accs:
    if max_accs < len(lsu_name_to_accs[i]):
        max_accs = len(lsu_name_to_accs[i])
print(F"Maximum number of accessions for LSU: {max_accs}")

acc_to_header = {}
for i in all_seq_dict:
    acc = i[1:].split(" ")[0]
    acc_to_header[acc] = i

long_acc_table = "Name\tAccession\tRegion\tHeader\tLength_Uncut\tSeq_Uncut\tLength_ITSx_Cut\tSeq_ITSx_Cut\n"
all_acc_to_name = {}

included_seqs_acc_to_header = {}
for header in included_seq_dict:
    acc = header[1:].split("|")[0]
    included_seqs_acc_to_header[acc] = header

for i in ssu_name_to_accs:
    for acc in ssu_name_to_accs[i]:
        header, reg = acc_to_header[acc].replace('\n', ''), "SSU"
        seq = all_seq_dict[acc_to_header[acc]][1].replace('\n', '')
        cut_seq = included_seq_dict[included_seqs_acc_to_header[acc]]["SSU"].replace('\n', '') if (acc in included_seqs_acc_to_header and "SSU" in included_seq_dict[included_seqs_acc_to_header[acc]]) else ""
        long_acc_table += F"{i}\t{acc}\t{reg}\t{header}\t{len(seq)}\t{seq}\t{len(cut_seq)}\t{cut_seq}\n"
        all_acc_to_name[acc] = i
for i in lsu_name_to_accs:
    for acc in lsu_name_to_accs[i]:
        header, reg = acc_to_header[acc].replace('\n', ''), "LSU"
        seq = all_seq_dict[acc_to_header[acc]][1].replace('\n', '')
        cut_seq = included_seq_dict[included_seqs_acc_to_header[acc]]["LSU"].replace('\n', '') if (acc in included_seqs_acc_to_header and "LSU" in included_seq_dict[included_seqs_acc_to_header[acc]]) else ""
        long_acc_table += F"{i}\t{acc}\t{reg}\t{header}\t{len(seq)}\t{seq}\t{len(cut_seq)}\t{cut_seq}\n"
        all_acc_to_name[acc] = i

with open("../data/phylo/info/SSU_LSU_accessions.txt", "w") as ofile:
    ofile.write(long_acc_table)

# Columns: "Name\tAccession\tRegion\tHeader\tLength_Uncut\tSeq_Uncut\tLength_ITSx_Cut\tSeq_ITSx_Cut\n"
acc_table_df = pd.read_csv("../data/phylo/info/SSU_LSU_accessions.txt", sep = "\t")
use_uncut_lsu = acc_table_df.query('Region == "LSU" and Length_ITSx_Cut < 100')
use_uncut_lsu = use_uncut_lsu[use_uncut_lsu["Header"].str.contains("internal|ITS") == False]
use_uncut_lsu = use_uncut_lsu[use_uncut_lsu["Accession"].str.contains("_") == False]
# print(use_uncut_lsu)
uncut_lsu_names = {}
for i in use_uncut_lsu['Accession']:
    row = use_uncut_lsu.loc[use_uncut_lsu['Accession'] == i,:]
    if row.Name.iloc[0] in uncut_lsu_names:
        if row.Length_Uncut.iloc[0] > uncut_lsu_names[row.Name.iloc[0]][0]:
            uncut_lsu_names[row.Name.iloc[0]] = [row.Length_Uncut.iloc[0], i, row.Seq_Uncut.iloc[0]]
    else:
        uncut_lsu_names[row.Name.iloc[0]] = [row.Length_Uncut.iloc[0], i, row.Seq_Uncut.iloc[0]]

use_cut_lsu = acc_table_df.query('Region == "LSU" and Length_ITSx_Cut >= 100')
use_cut_lsu = use_cut_lsu[use_cut_lsu["Accession"].str.contains("_") == False]
# print(use_cut_lsu)
cut_lsu_names = {}
for i in use_cut_lsu['Accession']:
    row = use_cut_lsu.loc[use_cut_lsu['Accession'] == i,:]
    if row.Name.iloc[0] in cut_lsu_names:
        if row.Length_ITSx_Cut.iloc[0] > cut_lsu_names[row.Name.iloc[0]][0]:
            cut_lsu_names[row.Name.iloc[0]] = [row.Length_ITSx_Cut.iloc[0], i, row.Seq_ITSx_Cut.iloc[0]]
    else:
        cut_lsu_names[row.Name.iloc[0]] = [row.Length_ITSx_Cut.iloc[0], i, row.Seq_ITSx_Cut.iloc[0]]

use_uncut_ssu = acc_table_df.query('Region == "SSU" and Length_ITSx_Cut < 100')
use_uncut_ssu = use_uncut_ssu[use_uncut_ssu["Header"].str.contains("internal|ITS") == False]
use_uncut_ssu = use_uncut_ssu[use_uncut_ssu["Accession"].str.contains("_") == False]
# print(use_uncut_ssu)
uncut_ssu_names = {}
for i in use_uncut_ssu['Accession']:
    row = use_uncut_ssu.loc[use_uncut_ssu['Accession'] == i,:]
    if row.Name.iloc[0] in uncut_ssu_names:
        if row.Length_Uncut.iloc[0] > uncut_ssu_names[row.Name.iloc[0]][0]:
            uncut_ssu_names[row.Name.iloc[0]] = [row.Length_Uncut.iloc[0], i, row.Seq_Uncut.iloc[0]]
    else:
        uncut_ssu_names[row.Name.iloc[0]] = [row.Length_Uncut.iloc[0], i, row.Seq_Uncut.iloc[0]]

use_cut_ssu = acc_table_df.query('Region == "SSU" and Length_ITSx_Cut >= 100')
use_cut_ssu = use_cut_ssu[use_cut_ssu["Accession"].str.contains("_") == False]
# print(use_cut_ssu)
cut_ssu_names = {}
for i in use_cut_ssu['Accession']:
    row = use_cut_ssu.loc[use_cut_ssu['Accession'] == i,:]
    if row.Name.iloc[0] in cut_ssu_names:
        if row.Length_ITSx_Cut.iloc[0] > cut_ssu_names[row.Name.iloc[0]][0]:
            cut_ssu_names[row.Name.iloc[0]] = [row.Length_ITSx_Cut.iloc[0], i, row.Seq_ITSx_Cut.iloc[0]]
    else:
        cut_ssu_names[row.Name.iloc[0]] = [row.Length_ITSx_Cut.iloc[0], i, row.Seq_ITSx_Cut.iloc[0]]

lsu_seqs_to_use = {}
lsu_uncut_only = uncut_lsu_names.keys() - cut_lsu_names.keys()
lsu_cut_only = cut_lsu_names.keys() - uncut_lsu_names.keys()
lsu_inter = cut_lsu_names.keys() & uncut_lsu_names.keys()

lsu_seqs_to_use.update({k : uncut_lsu_names[k] for k in lsu_uncut_only})
lsu_seqs_to_use.update({k : cut_lsu_names[k] for k in lsu_cut_only})
for i in lsu_inter:
    if uncut_lsu_names[i][0] >= cut_lsu_names[i][0]:
        lsu_seqs_to_use[i] = uncut_lsu_names[i]
    else:
        lsu_seqs_to_use[i] = cut_lsu_names[i]

ssu_seqs_to_use = {}
ssu_uncut_only = uncut_ssu_names.keys() - cut_ssu_names.keys()
ssu_cut_only = cut_ssu_names.keys() - uncut_ssu_names.keys()
ssu_inter = cut_ssu_names.keys() & uncut_ssu_names.keys()

ssu_seqs_to_use.update({k : uncut_ssu_names[k] for k in ssu_uncut_only})
ssu_seqs_to_use.update({k : cut_ssu_names[k] for k in ssu_cut_only})
for i in ssu_inter:
    if uncut_ssu_names[i][0] >= cut_ssu_names[i][0]:
        ssu_seqs_to_use[i] = uncut_ssu_names[i]
    else:
        ssu_seqs_to_use[i] = cut_ssu_names[i]

buffer = ""
for name in lsu_seqs_to_use:
    seq_len, acc, seq = lsu_seqs_to_use[name]
    buffer += F">{acc}|F|LSU Extracted LSU sequence ({seq_len} bp)\n{seq}\n"

with open("../data/phylo/unaligned_seqs/class.yeast_seqs.LSU.fasta", "w") as ofile:
    ofile.write(buffer)

buffer = ""
for name in ssu_seqs_to_use:
    seq_len, acc, seq = ssu_seqs_to_use[name]
    buffer += F">{acc}|F|SSU Extracted SSU sequence ({seq_len} bp)\n{seq}\n"

with open("../data/phylo/unaligned_seqs/class.yeast_seqs.SSU.fasta", "w") as ofile:
    ofile.write(buffer)
