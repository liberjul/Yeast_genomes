#!/bin/bash --login

for i in SSU ITS1 5_8S ITS2 LSU RPB1 RPB2 TEF1 CYTB
do
  clipkit ../data/phylo/aligned_seqs/class.yeast_seqs.$i.mafft.fasta -m gappy -o ../data/phylo/aligned_seqs/class.yeast_seqs.$i.mafft.clipkit.fasta
done
