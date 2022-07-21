#!/bin/bash --login

ITSx -i ../unaligned_seqs/its1_58s_its2.fasta -o ../unaligned_seqs/yeast_seqs -t fungi --save_regions ITS1,5.8S,ITS2
ITSx -i ../unaligned_seqs/ssu.fasta -o ../unaligned_seqs/yeast_seqs -t fungi --save_regions SSU
ITSx -i ../unaligned_seqs/lsu.fasta -o ../unaligned_seqs/yeast_seqs -t fungi --save_regions LSU

mv ../unaligned_seqs/rpb1.fasta ../unaligned_seqs/yeast_seqs.RPB1.fasta
mv ../unaligned_seqs/rpb2.fasta ../unaligned_seqs/yeast_seqs.RPB2.fasta
mv ../unaligned_seqs/tef1.fasta ../unaligned_seqs/yeast_seqs.TEF1.fasta
mv ../unaligned_seqs/cytb.fasta ../unaligned_seqs/yeast_seqs.CYTB.fasta
