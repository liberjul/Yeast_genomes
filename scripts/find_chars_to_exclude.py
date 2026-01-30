import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help="Alignment in FASTA format")
parser.add_argument("-t","--thresh", type=float, default = 0.5, help="Threshold for proportion of gaps at which a site with >= the prop. of gaps is excluded.")
args = parser.parse_args()

seq_dict = {}
min_len, max_len = 9999999, 0
with open(args.fasta, "r") as ifile:
    seq = ""
    lines = ifile.readlines()
    for line in lines:
        if line != "" and line[0] == ">":
            if seq != "":
                if len(seq) < min_len:
                    min_len = len(seq)
                if len(seq) > max_len:
                    max_len = len(seq)
                seq_dict[header] = seq
            header = line[1:]
            seq = ""
        elif line != "" and line != "\n":
            seq += line.replace("\n", "")
    if seq != "":
        if len(seq) < min_len:
            min_len = len(seq)
        if len(seq) > max_len:
            max_len = len(seq)
        seq_dict[header] = seq

if min_len != max_len:
    raise ValueError("Sequences in the FASTA file are not all the same length. Please align file before input.")

gaps = [0]*max_len
for i in range(max_len):
    for header in seq_dict:
        if seq_dict[header][i] == "-":
            gaps[i] += 1

seq_count = len(seq_dict.keys())
gaps_prop = [0.]*max_len
for i in range(max_len):
    gaps_prop[i] = gaps[i]/seq_count

# Sliding window exclusion ranges
range_tups = [] # inclusive range tuples
min, max = 0, 1
while min < max_len and max < max_len:
    if gaps_prop[min] >= args.thresh:
        while max < max_len and gaps_prop[max] >= args.thresh:
            max += 1
        range_tups.append((min, max))
        print(min, max)
        min = max + 1
        max = min + 1
    else:
        min += 1
        max = min + 1
range_str = "	EXSET * UNTITLED  = "
for i in range_tups:
    range_str += F"{i[0]+1} -  {i[1]} "
print(range_str.strip(",") + ";")
