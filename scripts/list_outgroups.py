outgroup_genera = set(["Pseudobensingtonia","Bensingtonia", "Pseudohyphozyma", "Ruinenia", "Crittendenia", "Pseudosterigmatospora", "Sterigmatospora", "Meniscomyces"])
buffer = ""
with open("../data/phylo/alignments/20220820_infile.phy", "r") as ifile:
    line = ifile.readline()
    while line != "":
        if line[0] == ">" and line[1:].split("_")[0] in outgroup_genera:
            buffer += line[1:].strip() + ","
        line = ifile.readline()
print(buffer)
