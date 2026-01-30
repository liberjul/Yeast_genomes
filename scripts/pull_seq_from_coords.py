import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help="input FASTA")
parser.add_argument("-o", "--output", type=str, help="output files path and prefix")
parser.add_argument("-r", "--record", type=str, help="record/header/scaffold/contig name")
parser.add_argument("-s", "--start", type=int, help="1-indexed start position")
parser.add_argument("-p", "--stop", type=int, help="1-indexed stop position, inclusive")
args = parser.parse_args()

if ">" in args.record:
    args.record = args.record.strip(">")
if '"' in args.record:
    args.record = args.record.strip('"')

print(F"{args.record}, {args.start}, {args.stop}")
hit_dict={args.record : [args.start, args.stop]}

buffer = ""
newline = ""

with open(args.fasta, "r") as ifile:
    line = ifile.readline()
    while line != "":
        if line.strip()[1:] in hit_dict:
            header = line.strip()[1:]
            seq = ""
            line = ifile.readline()
            while line != "" and line[0] != ">":
                seq += line.strip()
                line = ifile.readline()
            start, stop = hit_dict[header]
            print(header, " Length: ", len(seq))
            if int(start) > int(stop):
                stop_new = start
                start = stop
                stop = stop_new
            buffer = F"{buffer}>{header}_[{start}..{stop}]\n{seq[int(start)-1:int(stop)]}\n"
        else:
            line = ifile.readline()

with open(args.output, "w") as ofile:
    ofile.write(buffer)
