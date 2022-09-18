#!/bin/bash --login

DATADIR=../data/phylo/unaligned_seqs
for PREFIX in mbm outgroup
do
  ITSx -i $DATADIR/"$PREFIX".its1_58s_its2.fasta -o $DATADIR/"$PREFIX".itsx -t fungi --save_regions ITS1,5.8S,ITS2
  ITSx -i $DATADIR/"$PREFIX".ssu.fasta -o $DATADIR/"$PREFIX".itsx -t fungi --save_regions SSU
  ITSx -i $DATADIR/"$PREFIX".lsu.fasta -o $DATADIR/"$PREFIX".itsx -t fungi --save_regions LSU
  for i in SSU ITS1 5_8S ITS2 LSU
  do
    cat $DATADIR/"$PREFIX".itsx.$i.fasta >> $DATADIR/"$PREFIX".yeast_seqs.$i.fasta
  done
  cp $DATADIR/"$PREFIX".rpb1.fasta $DATADIR/"$PREFIX".yeast_seqs.RPB1.fasta
  cp $DATADIR/"$PREFIX".rpb2.fasta $DATADIR/"$PREFIX".yeast_seqs.RPB2.fasta
  cp $DATADIR/"$PREFIX".tef1.fasta $DATADIR/"$PREFIX".yeast_seqs.TEF1.fasta
  cp $DATADIR/"$PREFIX".cytb.fasta $DATADIR/"$PREFIX".yeast_seqs.CYTB.fasta
done
