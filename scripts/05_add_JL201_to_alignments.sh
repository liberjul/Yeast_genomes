#!/bin/bash --login


DATADIR=../data/phylo/unaligned_seqs

x=0
for i in SSU ITS1 5_8S ITS2 LSU
do
  sed "s/NODE_101_length_6643_cov_1305.952635/NS00000$x.1/" ../data/phylo/phylo_loci_JL201/JL201_ITSx."$i".fasta >> $DATADIR/yeast_seqs.$i.fasta
  x=$(( $x + 1 ))
done
echo ">NS000006.1 UNVERIFIED: Novel species strain JL201 RNA polymerase II largest subunit-like gene, partial sequence" >> $DATADIR/yeast_seqs.RPB1.fasta
tail -n+2 ../data/phylo/phylo_loci_JL201/RPB1.fasta >> $DATADIR/yeast_seqs.RPB1.fasta
echo ">NS000007.1 UNVERIFIED: Novel species strain JL201 RNA polymerase II second largest subunit-like gene, partial sequence" >> $DATADIR/yeast_seqs.RPB2.fasta
tail -n+2 ../data/phylo/phylo_loci_JL201/RPB2.fasta >> $DATADIR/yeast_seqs.RPB2.fasta
echo ">NS000008.1 UNVERIFIED: Novel species strain JL201 translation elongation factor 1-alpha-like gene, partial sequence" >> $DATADIR/yeast_seqs.TEF1.fasta
tail -n+2 ../data/phylo/phylo_loci_JL201/TEF1.fasta >> $DATADIR/yeast_seqs.TEF1.fasta
echo ">NS000009.1 UNVERIFIED: Novel species strain JL201 cytochrome b-like gene, partial sequence; mitochondrial" >> $DATADIR/yeast_seqs.CYTB.fasta
tail -n+2 ../data/phylo/phylo_loci_JL201/CYTB.fasta >> $DATADIR/yeast_seqs.CYTB.fasta
