import json, os, sys, subprocess
import pandas as pd

data_dir = "../data/phylo"

with open(F"{data_dir}/info/class.updated.accession_metadata.json", "r") as ifile:
    metadata_dict = json.load(ifile)

if os.path.exists(F"{data_dir}/info/manual.accession_metadata.json"):
    with open(F"{data_dir}/info/manual.accession_metadata.json", "r") as ifile:
        metadata_dict.update(json.load(ifile))

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

species_w_mult_strains = []
cleaned_syns_dict = {}
for i in gen_sp_dict:
    strain_list = list(gen_sp_dict[i])
    strain_set = set()
    for j in strain_list:
        strain_str = j.split(" ")
        strain_str = [x for x in strain_str if x not in ["culture", "strain", "isolate"]]
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
            # print(i)
            # print(strain_set)
# print(cleaned_syns_dict)

ofile = "strain_name_synonyms_NCBI_taxon_query.txt"
## Takes a while to run, so uncomment if needed
# subprocess.run(F'> {data_dir}/info/{ofile}', shell=True)
# for i in species_w_mult_strains:
#     subprocess.run(F'esearch -db taxonomy -query "{i}[ORGN]" | efetch -db taxonomy -mode xml | xtract -pattern Taxon -element TaxId ScientificName DispName | tail -n1 >> {data_dir}/info/{ofile}', shell=True)
# print(gen_sp_dict)
taxon_query_dict = {}
if os.path.exists(F"{data_dir}/info/{ofile}"):
    with open(F"{data_dir}/info/{ofile}", "r") as ifile:
        line = ifile.readline()
        while line != "":
            spl = line.split("\t")
            if len(spl) > 2:
                strain_syns = set([x.strip("\n").replace(":", "").replace(" ", "").replace("specimen-voucher", "") for x in spl[3:] if "&" not in x and "," not in x])
                taxon_query_dict[spl[1]] = strain_syns
                cleaned_syns_dict[spl[1]] = strain_syns
                strain_syns_list = list(strain_syns)
                if len(strain_syns) > 1 and strain_syns_list[0] not in synonym_dict:
                    synonym_dict[strain_syns_list[0]] = set(strain_syns_list[1:])
                    for j in range(1, len(strain_syns_list)):
                        # print(strain_syns_list[0:j] + strain_syns_list[j+1:])
                        if strain_syns_list[j] not in synonym_dict:
                            synonym_dict[strain_syns_list[j]] = set(strain_syns_list[0:j] + strain_syns_list[j+1:])
            line = ifile.readline()
# print(taxon_query_dict)
excluded_strains = []
for i in species_w_mult_strains:
    if i not in taxon_query_dict:
        # print("No taxon hit with automated query: ", i)
        strains_from_seqs = cleaned_syns_dict[i]
        strains_from_seqs = [x.replace("CGMCCAS", "CGMCC").replace("putative","") for x in strains_from_seqs]
        # print("\tFrom sequence: ",  strains_from_seqs)
        if i != "Novel sp.":
            exclude = [x for x in strains_from_seqs if x not in synonym_dict]
            excluded_strains += exclude
    else:
        all_match = True
        strains_from_seqs = cleaned_syns_dict[i]
        strains_from_seqs = [x.replace("CGMCCAS", "CGMCC").replace("putative", "") for x in strains_from_seqs]
        strains_from_taxon = taxon_query_dict[i]
        for strain in strains_from_seqs:
            if strain not in strains_from_taxon:
                all_match = False
        if not all_match:
            exclude = [x for x in strains_from_seqs if x not in synonym_dict]
            excluded_strains += exclude
            # print("Not all strains are synonyms: ", i,)
            # print("\tFrom sequence: ",  strains_from_seqs)
            # print("\tFrom taxonomy: ",  strains_from_taxon)
excluded_strains = set(excluded_strains)
print("Excluded strains:", excluded_strains)



species_gap_dict = {}
regions = ["SSU", "ITS1", "5_8S", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
for reg in regions:
    # buffer = ""
    with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.mafft.ordered.fasta", "r") as ifile:
        seq = ""
        line = ifile.readline()
        while line != "":
            if line != "" and line[0] == ">":
                if seq != "" and strain not in excluded_strains:
                    # buffer += F"{header}{seq}\n"
                    if species in species_gap_dict:
                        species_gap_dict[species] += seq.count("-")
                    else:
                        species_gap_dict[species] = seq.count("-")
                header = line
                header_short = line[1:].split("|")[0]
                strain = "".join(header_short.split(" ")[2:])
                for i in ["culture", "isolate", "putative", "strain"]:
                    strain = strain.replace(i, "")
                # print(strain)
                species = " ".join(header_short.split(" ")[0:2])
                seq = ""
            elif line != "":
                seq += line.strip()
            line = ifile.readline()
        if seq != "" and strain not in excluded_strains:
            # buffer += F"{header}{seq}\n"
            if species in species_gap_dict:
                species_gap_dict[species] += seq.count("-")
            else:
                species_gap_dict[species] = seq.count("-")
    # with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.mafft.ordered.types.fasta", "w") as ofile:
    #     ofile.write(buffer)

buffer = "Genus,Species,Gap_count\n"
for i in species_gap_dict:
    gen, sp = i.split(" ")
    buffer += F"{gen},{sp},{species_gap_dict[i]}\n"
with open(F"{data_dir}/info/species_gap_counts.csv", "w") as ofile:
    ofile.write(buffer)

species_gap_df = pd.read_csv(F"{data_dir}/info/species_gap_counts.csv")
species_included = set(["Novel sp."])
for genus in species_gap_df.Genus.unique():
    # print(genus)
    sub_df = species_gap_df.copy()[species_gap_df.Genus == genus]
    sub_df = sub_df.copy()[sub_df.Species != "sp."]
    sub_df["Name"] = sub_df.Genus + " " + sub_df.Species
    sub_df.sort_values(by = "Gap_count", inplace = True)
    if len(sub_df) > 4:
        species_included.update(set(sub_df.Name[0:4]))
    else:
        species_included.update(set(sub_df.Name[0:]))
    # print(sub_df)
print(species_included)
print(len(species_included))

# print(synonym_dict)
strains_used = set()
regions = ["SSU", "ITS1", "5_8S", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
for reg in regions:
    reg_strain_used = set()
    buffer = ""
    with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.mafft.ordered.fasta", "r") as ifile:
        seq = ""
        line = ifile.readline()
        while line != "":
            if line != "" and line[0] == ">":
                if seq != "" and strain not in excluded_strains and strain not in reg_strain_used and species in species_included:
                    buffer += F">{species} {strain}|{reg}\n{seq}\n"
                    reg_strain_used.add(strain)
                header = line
                header_short = line[1:].split("|")[0]
                temp_strain = "".join(header_short.split(" ")[2:])
                for i in ["culture", "isolate", "putative", "strain"]:
                    temp_strain = temp_strain.replace(i, "")
                if temp_strain not in strains_used and temp_strain in synonym_dict:
                    strain_decided = False
                    for i in synonym_dict[temp_strain]:
                        if i in strains_used:
                            strain = i
                            strain_decided = True
                    if not strain_decided:
                        strain = temp_strain
                else:
                    strain = temp_strain
                strains_used.add(strain)
                # print(strain)
                species = " ".join(header_short.split(" ")[0:2])
                seq = ""
            elif line != "":
                seq += line.strip()
            line = ifile.readline()
        if seq != "" and strain not in excluded_strains and strain not in reg_strain_used and species in species_included:
            buffer += F">{species} {strain}|{reg}\n{seq}\n"
            reg_strain_used.add(strain)
    with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.mafft.ordered.types.fasta", "w") as ofile:
        ofile.write(buffer)
# for i in strains_used:
#     print(i)
#     print(i in synonym_dict)
#     if i in synonym_dict:
#         print(synonym_dict[i])
