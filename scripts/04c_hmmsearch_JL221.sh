#!/bin/bash --login

DATADIR=../data/phylo

mkdir -p $DATADIR/phylo_loci_JL221

for i in RPB1 RPB2 TEF1 CYTB
do
  nhmmer --tblout $DATADIR/phylo_loci_JL221/"$i"_hmm.out --cpu 24 $DATADIR/hmms/$i.hmm ../data/assembly/JL221/scaffolds.fasta
  python pull_hmm_hits.py -f ../data/assembly/JL221/scaffolds.fasta \
    -m $DATADIR/phylo_loci_JL221/"$i"_hmm.out -o $DATADIR/phylo_loci_JL221/"$i".fasta -t
done
