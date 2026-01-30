#!/bin/bash --login

DATADIR=../data/phylo

taxa="NB124-2"
mkdir -p $DATADIR/phylo_loci_"$taxa"

for i in RPB1 RPB2 TEF1 CYTB SSU 5_8S LSU
do
  nhmmer --tblout $DATADIR/phylo_loci_"$taxa"/"$i"_hmm.out --dna --cpu 24 $DATADIR/hmms/$i.hmm ../data/assembly/"$taxa"/scaffolds_polished.fasta
  python pull_hmm_hits.py -f ../data/assembly/"$taxa"/scaffolds_polished.fasta \
    -m $DATADIR/phylo_loci_"$taxa"/"$i"_hmm.out -o $DATADIR/phylo_loci_"$taxa"/"$i".fasta -t -s 2000 -p 2000
done
