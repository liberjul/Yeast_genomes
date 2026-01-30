import sys, subprocess, json
ofile = "../data/phylo/info/acc_strain_pairs.json"

with open(ofile, "r") as input:
    acc_dict = json.load(input)

with open(sys.argv[1], "r") as ifile:
    line = ifile.readline()
    while line != "":
        print(line)
        print(F'efetch -db nuccore -id "{line.strip()}" -api_key=95c4f8e7f8558cac9bbb733cb87ea362ff09 -format gb > temp.out')
        subprocess.run(F'efetch -db nuccore -id "{line.strip()}" -api_key=95c4f8e7f8558cac9bbb733cb87ea362ff09 -format gb > temp.out', shell=True)
        with open("temp.out", "r") as res:
            lines = res.readlines()
            lines = [x for x in lines if "strain=" in x or "isolate=" in x or "culture_collection=" in x or "specimen_voucher=" in x]
            strains = [x.split("=")[1].strip() for x in lines]
            strain_str = ";".join(strains)
            strains = [x.strip(" ").strip("'").strip('"') for x in strain_str.split(";")]
        acc_dict[line.strip()] = strains
        print(strains)
        line = ifile.readline()

with open(ofile, "w") as output:
    json.dump(acc_dict, output)
