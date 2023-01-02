#!/bin/bash --login

DATADIR=../data/phylo


for i in CYTB
do
  for ASSEM in $DATADIR/jgi_genomic_seqs/*ito*fasta
  do
    PREFIX=$(basename -- ${ASSEM%_*})
    echo $PREFIX
    nhmmer --tblout $DATADIR/jgi_genomic_seqs/"$i"_hmm.out --cpu 24 $DATADIR/hmms/$i.hmm $ASSEM
    python pull_hmm_hits.py -f $ASSEM \
      -m $DATADIR/jgi_genomic_seqs/"$i"_hmm.out -o ${ASSEM%fasta}"$i".fasta -t
    sed -i "s/scaffold/"$PREFIX"_scaffold/" ${ASSEM%fasta}"$i".fasta
  done
done
