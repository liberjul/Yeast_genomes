#!/bin/bash --login

for i in SSU ITS1 5_8S ITS2 LSU RPB1 RPB2 TEF1 CYTB
do
  clipkit ../aligned_seqs/yeast_seqs.$i.mafft.fasta -m gappy -o ../aligned_seqs/yeast_seqs.$i.mafft.clipkit.fasta
done
