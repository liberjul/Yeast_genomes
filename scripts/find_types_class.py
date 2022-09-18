import json, os

data_dir = "../data/phylo"

types = set()
with open(F"{data_dir}/info/type_strains.txt", "r") as ifile:
    for i in ifile.readlines():
        spl = i.strip().split()
        spl[0] = spl[0].strip("T")
        for j in spl:
            types.add(j)
            types.add(j.replace(" ", ""))


with open(F"{data_dir}/info/class.updated.accession_metadata.json", "r") as ifile:
    metadata_dict = json.load(ifile)

if os.path.exists(F"{data_dir}/info/manual.accession_metadata.json"):
    with open(F"{data_dir}/info/manual.accession_metadata.json", "r") as ifile:
        metadata_dict.update(json.load(ifile))

type_dict = {}
species_set = set()
type_metadata_dict = {}
for i in metadata_dict:
    if metadata_dict[i][1] != "sp.":
        species = " ".join(metadata_dict[i][0:2])
        species_set.add(species)
        if metadata_dict[i][2] != "":
            strain_str = metadata_dict[i][2].split()
            if len(strain_str) > 2:
                strain_str = [x for x in strain_str if x not in ["culture", "strain", "isolate"]]
                strain = strain_str[0]
            elif strain_str[0] in ["culture", "strain", "isolate"]:
                strain = strain_str[1]
            else:
                strain = strain_str[0]
            # test if the strain is a type
            if strain in types:
                type_dict[species] = strain
                type_metadata_dict[i] = metadata_dict[i]
    # elif metadata_dict[i][0] == "Novel":
    #     print(metadata_dict[i])
type_keys_set = set(type_dict.keys())
print("No type strain found:")
print("\n".join(list(species_set - type_keys_set)))
print("Number of type strains: ", len(type_keys_set))

type_headers = set(["Novel sp. strain JL201", "Novel sp. strain JL221"])
for i in type_metadata_dict:
    type_headers.add(" ".join(type_metadata_dict[i]))
print(type_headers)
### Filter alignments to include only types and novel strains

# regions = ["SSU", "ITS1", "5_8S", "ITS2", "LSU", "RPB1", "RPB2", "TEF1", "CYTB"]
# for reg in regions:
#     buffer = ""
#     with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.mafft.ordered.fasta", "r") as ifile:
#         seq = ""
#         line = ifile.readline()
#         while line != "":
#             if line != "" and line[0] == ">":
#                 if seq != "" and header_short in type_headers:
#                     buffer += F"{header}{seq}\n"
#                 header = line
#                 header_short = line[1:].split("|")[0]
#                 seq = ""
#             elif line != "":
#                 seq += line.strip()
#             line = ifile.readline()
#         if seq != "" and header_short in type_headers:
#             buffer += F"{header}{seq}\n"
#     with open(F"{data_dir}/aligned_seqs/class.yeast_seqs.{reg}.mafft.ordered.types.fasta", "w") as ofile:
#         ofile.write(buffer)
