name_list = []
with open("../alignments/2022_07_19_JL201_221_257_ITS_refseq_blast_agbm.fas", "r") as ifile:
    line = ifile.readline()
    count = 1
    while line != "":
        if line[0] == ">":
            name = line[1:].strip()
            name_list.append([F"t{count}", name])
            count += 1
        line = ifile.readline()

with open("20220719_raxml_JL201-221-257_ITS_agbm.raxml.support", "r") as ifile:
    text = ifile.read()
    for i in name_list:
        text = text.replace(F"{i[0]}:", F"{i[1]}:")

with open("20220719_raxml_JL201-221-257_ITS_agbm_names.raxml.support", "w") as ofile:
    ofile.write(text)
