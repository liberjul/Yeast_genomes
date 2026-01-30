import argparse
import numpy as np

def gap_analyzer(array):
    if array[0] == 0 or array[-1] == 0 or np.sum(array[1:-1] != 0) != 0:
        raise ValueError("Array must be zeroes in the middle bounded by -1 or 1 as the first and last elements.")
    else:
        if array[0] == -1 and array[-1] == 1:
            dir = "divergent"
        elif array[0] == 1 and array[-1] == -1:
            dir = "convergent"
        elif (array[0] == 1 and array[-1] == 1) or (array[0] == -1 and array[-1] == -1):
            dir = "tandem"
        else:
            raise ValueError("Array must be zeroes in the middle bounded by -1 or 1 as the first and last elements.")
        return(dir, len(array)-2)

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help="genome FASTA file")
parser.add_argument("-g", "--gff", type=str, help="genome GFF file, with gene annotations")
parser.add_argument("-o", "--output", type=str, help="output files path and prefix")
args = parser.parse_args()

scaf_len_dict = {}
scaf_arr_dict = {}
with open(args.fasta, "r") as ifile:
    seq = ""
    line = ifile.readline()
    while line != "":
        if line[0] == ">":
            if seq != "" and len(seq) > 1000:
                scaf_len_dict[header] = len(seq)
                scaf_arr_dict[header] = np.zeros(len(seq), dtype=np.int8)
            header = line[1:].strip()
            seq = ""
        else:
            seq += line.strip()
        line = ifile.readline()
    if seq != "" and len(seq) > 1000:
        scaf_len_dict[header] = len(seq)
        scaf_arr_dict[header] = np.zeros(len(seq), dtype=np.int8)

# for i in scaf_len_dict:
#     print(i, scaf_len_dict[i])

dir_d = {"+" : 1, "-" : -1}
with open(args.gff, "r") as ifile:
    line = ifile.readline()
    while line != "":
        if "RNA\t" in line:
            dir = dir_d[line.split("\t")[6]]
            scaf = line.split("\t")[0]
            start, stop = line.split("\t")[3:5]
            if scaf in scaf_len_dict:
                scaf_arr_dict[scaf][int(start)-1:int(stop)] = dir
        line = ifile.readline()

buffer = "scaffold,start,stop,orientation,gap_size\n"
for i in scaf_arr_dict:
    # print(i, np.sum(scaf_arr_dict[i] != 0))
    z_loc = np.where(scaf_arr_dict[i] != 0)[0]
    x = 0
    while x+1 < len(z_loc):
        # print(z_loc[x])
        if z_loc[x] + 1 != z_loc[x+1]:
            roi = scaf_arr_dict[i][z_loc[x]:z_loc[x+1]+1]
            orientation, gap_size = gap_analyzer(roi)
            buffer += F"{i},{z_loc[x]+2},{z_loc[x+1]+1},{orientation},{gap_size}\n"
            # print(scaf_arr_dict[i][z_loc[x]:z_loc[x+1]+1])
        x += 1

with open(F"{args.output}_gap.csv", "w") as ofile:
    ofile.write(buffer)
