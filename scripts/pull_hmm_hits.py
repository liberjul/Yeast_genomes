import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help="input FASTA")
parser.add_argument("-m", "--hmm", type=str, help="input hmm hit table, from nhmmer")
parser.add_argument("-o", "--output", type=str, help="output files path and prefix")
parser.add_argument("-t", "--top_hit_only", action="store_true", help="output top hit only")
args = parser.parse_args()

hit_dict={}
min_e = 9999.
with open(args.hmm, "r") as ifile:
    line = ifile.readline()
    while line != "" and line[0] == "#":
        line = ifile.readline()
    # count = 0
    while line != "" and line[0] != "#":
        # count += 1
        spl = line.strip().split(" ")
        spl = [x for x  in spl if x != ""]
        header, start, stop, e = spl[0], spl[6], spl[7], spl[12]
        # print([header, start, stop, e])
        if float(e) < min_e:
            min_e = float(e)
        if header in hit_dict:
            hit_dict[header].append([start, stop, e])
        else:
            hit_dict[header] = [[start, stop, e]]
        line = ifile.readline()
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
            for i in hit_dict[header]:
                start, stop, e = i
                if args.top_hit_only:
                    if float(e) <= min_e:
                        print(header, " Length: ", len(seq))
                        if int(start) > int(stop):
                            stop_new = start
                            start = stop
                            stop = stop_new
                        buffer = F"{buffer}>{header}_[{start}..{stop}]_e={e}\n{seq[int(start)-1:int(stop)]}\n"
                else:
                    print(header, " Length: ", len(seq))
                    if int(start) > int(stop):
                        stop_new = start
                        start = stop
                        stop = stop_new
                    buffer = F"{buffer}>{header}_[{start}..{stop}]_e={e}\n{seq[int(start)-1:int(stop)]}\n"
        else:
            line = ifile.readline()

with open(args.output, "w") as ofile:
    ofile.write(buffer)
