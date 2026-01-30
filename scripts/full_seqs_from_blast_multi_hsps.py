#!/usr/bin/env python
import argparse, subprocess

'''
cd /mnt/c/Users/jal138/Duke\ Bio_Ea\ Dropbox/He\ Lab/Julian_Liber/Microbial_traits/Aureobasidium/Genome_seqs/
blastn -query EXF150_70454.fasta -db ../GWAS/data/Ap_genomes \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" > 70454_blast.out
python ~/bin/code_utilities/full_seqs_from_blast_multi_hsps.py \
  -i 70454_blast.out -d ../GWAS/data/Ap_genomes -o 70454_aln.fasta

'''
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--db", type=str, help="BLAST database")
parser.add_argument("-i", "--input", type=str, help="Input, BLAST result, outfmt 6")
parser.add_argument("-o", "--output", type=str, help="Output file name, FASTA format")
parser.add_argument("-e", "--evalue", type=float, default = 1., help="E-value threshold, includes hits <= threshold.")
args = parser.parse_args()

hit_dict = {}
with open(args.input, "r") as ifile:
    line = ifile.readline()
    while line != "":
        spl = line.split("\t")
        sseqid, qstart, qend, sstart, send, evalue = spl[1, 6, 7, 8, 9, 10]
        if evalue <= args.evalue:
            if sseqid in hit_dict:
                hit_dict[sseqid]["sstarts"].append(sstart)
                hit_dict[sseqid]["sends"].append(send)
                hit_dict[sseqid]["evals"].append(evalue)
            else:
                hit_dict[sseqid] = {"sstarts" : [sstart], "sends" : [send], "evals" : [evalue]}
        line = ifile.readline()
seqs_to_pull = {}
for sseqid in hit_dict:
    if min(hit_dict[sseqid]["sstarts"]) < min(hit_dict[sseqid]["sends"]): #strand is +
        seqs_to_pull[sseqid] = [min(hit_dict[sseqid]["sstarts"]), max(hit_dict[sseqid]["sends"]), "+"]
    else: #strand is -
        seqs_to_pull[sseqid] = [min(hit_dict[sseqid]["sends"]), max(hit_dict[sseqid]["sstarts"]), "-"]
print(seqs_to_pull)
# subprocess.run(F"makeblastdb -in {db}.fasta -dbtype nucl -out {db}",
        # shell=True, executable="/bin/bash")