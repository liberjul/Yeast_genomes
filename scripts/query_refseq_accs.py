import sys, subprocess, json

ofile = "../data/phylo/info/refseq_nucleotide_pairs.json"
acc_dict = {}
with open(sys.argv[1], "r") as ifile:
    line = ifile.readline()
    while line != "":
        print(line)
        print(F'esearch -db nucleotide -query "{line.strip()}" -api_key=95c4f8e7f8558cac9bbb733cb87ea362ff09 | efetch -db nucleotide -mode xml | xtract -pattern Seqdesc -element User-field_data_str | grep -oP "[A-Z]{{2}}[0-9]{{6}}" > temp.out')
        subprocess.run(F'esearch -db nucleotide -query "{line.strip()}" -api_key=95c4f8e7f8558cac9bbb733cb87ea362ff09 | efetch -db nucleotide -mode xml | xtract -pattern Seqdesc -element User-field_data_str | grep -oP "[A-Z]{{2}}[0-9]{{6}}" > temp.out', shell=True)
        with open("temp.out", "r") as res:
            acc = res.readline().strip()
        acc_dict[line.strip()] = acc
        print(acc)
        line = ifile.readline()

with open(ofile, "w") as output:
    json.dump(acc_dict, output)
