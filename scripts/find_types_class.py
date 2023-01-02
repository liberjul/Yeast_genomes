import json, os
from datetime import date
from Bio.Nexus import Nexus

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
def remove_duplicate(seq_dict, dup_accs):
    dup_groups = {}
    excluded_headers = []
    for i in dup_accs:
        binom = "_".join(i[1:3])
        if binom in dup_groups:
            dup_groups[binom].append(i)
        else:
            dup_groups[binom] = [i]
    for i in dup_groups:
        accs = [x[0] for x in dup_groups[i]]
        max_info_chars = 0
        sub_seq_dict = {}
        for header in seq_dict:
            for acc in accs:
                if acc in header:
                    sub_seq_dict[header] = seq_dict[header]
                    seq = seq_dict[header][1]
                    seq.replace("-", "").replace("?", "")
                    if max_info_chars < len(seq):
                        max_info_chars = len(seq)
        included_dict = {}
        for header in sub_seq_dict:
            seq = sub_seq_dict[header][1]
            seq.replace("-", "").replace("?", "")
            too_short, refseq, wrong_name, no_sp = 0, 0, 0, 0
            if len(seq) < max_info_chars:
                too_short += 1
                excluded_headers.append(header)
            elif "NR_" in header:
                refseq += 1
                excluded_headers.append(header)
            elif reg in ["RPB1", "RPB2", "CYTB", "TEF1"] and i.replace("|", "_") not in header:
                wrong_name += 1
                excluded_headers.append(header)
            elif i != "Novel_sp." and "_sp." in i:
                no_sp += 1
                excluded_headers.append(header)
            else:
                included_dict[header] = sub_seq_dict[header]
        if len(included_dict) > 1:
            # print("More than one header for species ", i, included_dict.keys())
            chosen_key = list(included_dict.keys())[0]
            excluded_headers.extend([x for x in included_dict if chosen_key != x])
        elif len(included_dict) == 0:
            print("No included record for species ", i, F"Too short: {too_short}, Refseq acc: {refseq}, Wrong name: {wrong_name}, No species: {no_sp}")
    return excluded_headers
def remove_others(seq_dict):
    excluded_headers = []
    for header in seq_dict:
        if "NR_" in header:
            excluded_headers.append(header)
    return excluded_headers

def remove_excluded_chars(seq_dict, exclusion_str):
    exclusion_str = exclusion_str.strip().strip(";").split(" = ")[1]
    split_list = exclusion_str.replace(" -  ", "-").split()
    split_list = [[x,x] if "-" not in x else x.split("-") for x in split_list]
    # print(split_list)
    exclusion_list = []
    for i in split_list:
        exclusion_list.append([int(x) for x in i if x != ""])
    out_dict = {}
    for i in seq_dict:
        seq = list(seq_dict[i][1])
        for x in exclusion_list:
            min, max = x
            seq[min-1: max] =  ["-"]*len(seq[min-1: max])
        out_dict[i] = seq_dict[i][0], "".join(seq)
    return out_dict

def exclude_chars_and_update_partitions(nex_file):
    gaponly = nex_file.gaponly()
    # nex_file.matrix = nex_file.crop_matrix(exclude = gaponly)
    if "UNTITLED" in nex_file.charpartitions:
        partitions = nex_file.charpartitions["UNTITLED"]
        gaps_sparse = [1 if x in gaponly else 0 for x in range(nex_file.nchar)]
        cumul_gaps = [0]*nex_file.nchar
        gaps_count = 0
        for i in range(nex_file.nchar):
            gaps_count += gaps_sparse[i]
            cumul_gaps[i] = gaps_count
        new_partitions = {"NEW_POS" : {x : [] for x in partitions}, "OLD_POS" : {x : [] for x in partitions}}
        for i in partitions:
            for j in partitions[i]:
                if j not in gaponly:
                    new_partitions["NEW_POS"][i].append(j - cumul_gaps[j])
                    new_partitions["OLD_POS"][i].append(j)
        # nex_file.charpartitions["UNTITLED"] = new_partitions["NEW_POS"]
        # nex_file.charpartitions["OLD_POS"] = new_partitions["OLD_POS"]
    # nex_file.nchar = nex_file.nchar - len(gaponly)
        return new_partitions
    else:
        return {}


data_dir = "../data/phylo"

# types = set()
# with open(F"{data_dir}/info/type_strains.txt", "r") as ifile:
#     for i in ifile.readlines():
#         spl = i.strip().split()
#         spl[0] = spl[0].strip("T")
#         for j in spl:
#             types.add(j)
#             types.add(j.replace(" ", ""))


with open(F"{data_dir}/info/accession_metadata.json", "r") as ifile:
    metadata_dict = json.load(ifile)

cleaned_syns_dict = {}
if os.path.exists(F"{data_dir}/info/cleaned_syns_dict.json"):
    with open(F"{data_dir}/info/cleaned_syns_dict.json", "r") as ifile:
        new_dict = json.load(ifile)
        for i in new_dict:
            new_dict[i] = set(new_dict[i])
        cleaned_syns_dict.update(new_dict)
if os.path.exists(F"{data_dir}/info/acc_strain_pairs.json"):
    with open(F"{data_dir}/info/acc_strain_pairs.json", "r") as ifile:
        acc_strain_dict = json.load(ifile)


# type_dict = {}
# species_set = set()
# type_metadata_dict = {}
# for i in metadata_dict:
#     if metadata_dict[i][1] != "sp.":
#         species = " ".join(metadata_dict[i][0:2])
#         species_set.add(species)
#         if metadata_dict[i][2] != "":
#             strain_str = metadata_dict[i][2].split()
#             if len(strain_str) > 2:
#                 strain_str = [x for x in strain_str if x not in ["culture", "strain", "isolate"]]
#                 strain = strain_str[0]
#             elif strain_str[0] in ["culture", "strain", "isolate"]:
#                 strain = strain_str[1]
#             else:
#                 strain = strain_str[0]
#             # test if the strain is a type
#             if strain in types:
#                 type_dict[species] = strain
#                 type_metadata_dict[i] = metadata_dict[i]
#     # elif metadata_dict[i][0] == "Novel":
#     #     print(metadata_dict[i])
# type_keys_set = set(type_dict.keys())
# print("No type strain found:")
# print("\n".join(list(species_set - type_keys_set)))
# print("Number of type strains: ", len(type_keys_set))
#
# type_headers = set(["Novel sp. strain JL201", "Novel sp. strain JL221"])
# for i in type_metadata_dict:
#     type_headers.add(" ".join(type_metadata_dict[i]))
# print(type_headers)

### Filter alignments to include only types and novel strains

all_non_type_accs = set()
# regions = ["SSU"]#, "ITS1", "5_8S", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
# regions = ["5_8S"]#, "SSU", "ITS1", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
# regions = ["5_8S", "SSU", "ITS1", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
regions = ["5_8S", "SSU", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
for reg in regions:
    exclusion_str = ""
    print("REGION = ", reg)
    buf1, buf2, buf3 = "", "", ""
    with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.mafft.manual.fas.nex", "r") as ifile:
        line = ifile.readline()
        while "TAXLABELS" not in line:
            buf1 += line
            line = ifile.readline()
        buf1 += line
        line = ifile.readline()
        taxa = line.strip().strip("\t").split("' '")
        if "|" in taxa[0]:
            accs = [x.strip("'").split("|")[0] for x in taxa]
        else:
            accs = [x.strip("'").split("_")[0] for x in taxa]
        accs_not_in_dict = [x for x in accs if x not in metadata_dict]
        print(accs_not_in_dict)
        metadata_strs = [([x] + metadata_dict[x][:2]+[clean_strain(metadata_dict[x][2])]) if metadata_dict[x][2] != "" else ([x] + metadata_dict[x][:2]+[clean_strain(acc_strain_dict[x][0])]) for x in accs]
        binoms_set = set()
        dup_set = set()
        header_binom = {}
        for x, header in zip(metadata_strs, taxa):
            if "NR_" not in x[0]:
                binom = "_".join(x[1:3])
                if binom == "Novel_sp." or binom == "Microbotryomycetes_sp.":
                    header_binom[header.strip("'")] = binom.strip(".") + "_" + x[3]
                else:
                    header_binom[header.strip("'")] = binom
                if binom not in binoms_set or binom == "Novel_sp.":
                    binoms_set.add(binom)
                else:
                    dup_set.add(binom)
        if reg == "5_8S":
            all_binoms = list(set(header_binom.values()))
        # print(all_binoms)
        dup_accs = [x for x in metadata_strs if "|".join(x[1:3]) in dup_set]
        type_accs = [x[0] for x in metadata_strs if x[3] in cleaned_syns_dict[" ".join(x[1:3])]]
        non_type_accs = ["|".join(x) for x in metadata_strs if x[3] not in cleaned_syns_dict[" ".join(x[1:3])]]
        # print("Non-type accessions:", non_type_accs)
        # print("Duplicate species accessions:", dup_accs)
        all_non_type_accs.update(set(non_type_accs))

        while "\tMATRIX" not in line:
            line = ifile.readline()
            buf2 += line
        line = ifile.readline()
        seq_dict = {}
        while line != "\n":
            # print([line])
            header, seq = line.strip().split("' ")[:2]
            seq = seq.strip(" ")
            header = header.strip("'")
            acc = header.split("|")[0] if "|" in header else header.split("_")[0]
            seq_dict[header] = [acc, seq]
            line = ifile.readline()
        while line != "" and "BEGIN MESQUITECHARMODELS;" not in line:
            buf3 += line
            if "EXSET" in line:
                exclusion_str = line
            line = ifile.readline()
    ## Resolve extra_seqs
    excluded_headers = remove_duplicate(seq_dict, dup_accs)
    excluded_headers.extend(remove_others(seq_dict))
    if exclusion_str != "":
        print(F"Changing excluded characters to gaps for {reg}")
        seq_dict = remove_excluded_chars(seq_dict, exclusion_str)
    remaining_headers = {header_binom[x] : seq_dict[x] for x in seq_dict if x not in excluded_headers}
    nchar = int(buf2.split("NCHAR=")[1].split(";")[0])
    exported_count = 0
    seq_buffer = ""
    for header in all_binoms:
        if header in remaining_headers:
            seq_buffer += F"\t'{header}'            {remaining_headers[header][1]}\n"
            exported_count += 1
        # else:
        #     buffer += F"\t'{header}'            {'?'*nchar}\n"
    buffer = buf1.replace(F"NTAX={len(seq_dict)}", F"NTAX={exported_count}")
    buffer += "\t'" + "' '".join(remaining_headers) + "'\n" + buf2.replace("TITLE  Character_Matrix", F"TITLE  {reg}")
    buffer += seq_buffer + buf3
    print(F"{reg}: count remaining headers: {len(remaining_headers)}, count exported: {exported_count}")
    with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.dedup_types.nex", "w") as ofile:
        ofile.write(buffer)
    with open(F"{data_dir}/aligned_seqs/{reg}.nex", "w") as ofile:
        ofile.write(buffer)

partition_buffer = "#nexus\nbegin sets;\n"
partition_count = 1

nexi = [[reg, Nexus.Nexus(F"../data/phylo/aligned_seqs/class.yeast_seqs.{reg}.dedup_types.nex")] for reg in regions]
for i in nexi:
    print(F"Total characters before exclusion, {i[0]}: ", i[1].nchar)
    # print(F"Number of gap only characters before exclusion, {i[0]}: ", len(i[1].gaponly()))
    trimmed = exclude_chars_and_update_partitions(i[1])
    # print(F"Total characters after exclusion, {i[0]}: ", trimmed.nchar)
    # print(F"Number of gap only characters after exclusion, {i[0]}: ", len(trimmed.gaponly()))
    print(trimmed)
    # print("Rhodotorula mucilaginosa in taxa: ")
    # if "Rhodotorula_mucilaginosa" in i[1].matrix:
    #     print(i[1].matrix["Rhodotorula_mucilaginosa"])
    print(i[0])
    file = F"../data/phylo/aligned_seqs/class.yeast_seqs.{i[0]}.excl.nex"
    i[1].write_nexus_data(filename=open(file, "w"), exclude = i[1].gaponly())
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
#charset part1 = dna.phy:



buffer = "\n".join([x.split("|")[0] for x in all_non_type_accs])
with open("../data/phylo/info/accs_wo_matching_strain.txt", "w") as ofile:
    ofile.write(buffer)

    #     seq = ""
    #     line = ifile.readline()
    #     while line != "":
    #         if line != "" and line[0] == ">":
    #             if seq != "" and header_short in type_headers:
    #                 buffer += F"{header}{seq}\n"
    #             header = line
    #             header_short = line[1:].split("|")[0]
    #             seq = ""
    #         elif line != "":
    #             seq += line.strip()
    #         line = ifile.readline()
    #     if seq != "" and header_short in type_headers:
    #         buffer += F"{header}{seq}\n"
    # with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.mafft.ordered.types.fasta", "w") as ofile:
    #     ofile.write(buffer)
