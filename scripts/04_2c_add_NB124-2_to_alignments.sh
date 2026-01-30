#!/bin/bash --login

DATADIR=../data/phylo

taxa="NB124-2"

for i in RPB1 RPB2 TEF1 CYTB SSU 5_8S LSU
do
  echo ">Novel_sp_"$taxa > temp.fna
  head -n2 $DATADIR/phylo_loci_"$taxa"/"$i".fasta | tail -n1 >> temp.fna
  mafft --keeplength \
    --add $DATADIR/phylo_loci_"$taxa"/"$i".fasta \
	--maxiterate 1000 \
	--op 0.0 \
	$DATADIR/aligned_seqs/class.yeast_seqs."$i".dedup_types_no_miss.fasta > \
	$DATADIR/aligned_seqs/class.yeast_seqs."$i".dedup_types_plus_new.fasta
  rm temp.fna
done