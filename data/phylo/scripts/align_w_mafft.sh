#!/bin/bash --login

for i in SSU ITS1 5_8S ITS2 LSU RPB1 RPB2 TEF1 CYTB
do
  mafft ../unaligned_seqs/yeast_seqs.$i.fasta > ../aligned_seqs/yeast_seqs.$i.mafft.fasta
done
