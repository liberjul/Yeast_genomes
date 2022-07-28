#!/bin/bash --login

ITSx -i ../unaligned_seqs/its1_58s_its2.fasta -o ../unaligned_seqs/itsx -t fungi --save_regions ITS1,5.8S,ITS2
ITSx -i ../unaligned_seqs/ssu.fasta -o ../unaligned_seqs/itsx -t fungi --save_regions SSU
ITSx -i ../unaligned_seqs/lsu.fasta -o ../unaligned_seqs/itsx -t fungi --save_regions LSU

for i in SSU ITS1 5_8S ITS2 LSU
do
  cp ../unaligned_seqs/itsx.$i.fasta ../unaligned_seqs/yeast_seqs.$i.fasta
done

cp ../unaligned_seqs/rpb1.fasta ../unaligned_seqs/yeast_seqs.RPB1.fasta
cp ../unaligned_seqs/rpb2.fasta ../unaligned_seqs/yeast_seqs.RPB2.fasta
cp ../unaligned_seqs/tef1.fasta ../unaligned_seqs/yeast_seqs.TEF1.fasta
cp ../unaligned_seqs/cytb.fasta ../unaligned_seqs/yeast_seqs.CYTB.fasta
