#!/bin/bash --login

DATADIR=../data/phylo/unaligned_seqs
ITSx -i $DATADIR/its1_58s_its2.fasta -o $DATADIR/itsx -t fungi --save_regions ITS1,5.8S,ITS2
ITSx -i $DATADIR/ssu.fasta -o $DATADIR/itsx -t fungi --save_regions SSU
ITSx -i $DATADIR/lsu.fasta -o $DATADIR/itsx -t fungi --save_regions LSU

ITSx -i ../data/assembly/JL201/NODE_101.fasta -o ../data/phylo/phylo_loci_JL201/JL201_ITSx -t fungi --save_regions SSU,ITS1,5.8S,ITS2,LSU
ITSx -i ../data/assembly/JL221/NODE_120.fasta -o ../data/phylo/phylo_loci_JL221/JL221_ITSx -t fungi --save_regions SSU,ITS1,5.8S,ITS2,LSU

for strain in JL201 JL221
do
  for i in SSU ITS1 5_8S ITS2 LSU
  do
    cp $DATADIR/../phylo_loci_$strain/"$strain"_ITSx.$i.fasta $DATADIR/../phylo_loci_$strain/yeast_seqs.$i.fasta
  done
done

cp $DATADIR/rpb1.fasta $DATADIR/yeast_seqs.RPB1.fasta
cp $DATADIR/rpb2.fasta $DATADIR/yeast_seqs.RPB2.fasta
cp $DATADIR/tef1.fasta $DATADIR/yeast_seqs.TEF1.fasta
cp $DATADIR/cytb.fasta $DATADIR/yeast_seqs.CYTB.fasta
