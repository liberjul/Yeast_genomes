#!/bin/bash --login

DATADIR=../data/phylo
OUT=RS_NRRLY17275
OUTDIR=$DATADIR/phylo_loci_$OUT
GENOME=../data/assembly/RS_NRRLY17275.fna
mkdir -p $OUTDIR

for i in RPB1 RPB2 TEF1 CYTB LSU SSU 5_8S
do
  # nhmmer --tblout $OUTDIR/"$i"_hmm.out --cpu 4 $DATADIR/hmms/$i.hmm $GENOME
  python pull_hmm_hits.py -f $GENOME \
    -m $OUTDIR/"$i"_hmm.out -o $OUTDIR/"$i".fasta -t -s 1000 -p 1000
done
