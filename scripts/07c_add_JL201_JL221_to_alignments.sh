#!/bin/bash --login

DATADIR=../data/phylo/unaligned_seqs

python 05c_split_fastas_by_reg.py --new

ITSx -i $DATADIR/new.its1_58s_its2.fasta -o $DATADIR/new.itsx -t fungi --save_regions ITS1,5.8S,ITS2
ITSx -i $DATADIR/new.ssu.fasta -o $DATADIR/new.itsx -t fungi --save_regions SSU
ITSx -i $DATADIR/new.lsu.fasta -o $DATADIR/new.itsx -t fungi --save_regions LSU

for i in SSU ITS1 5_8S ITS2 LSU
do
  cat $DATADIR/new.itsx.$i.fasta >> $DATADIR/yeast_seqs.$i.fasta
done

cat $DATADIR/new.rpb1.fasta >> $DATADIR/yeast_seqs.RPB1.fasta
cat $DATADIR/new.rpb2.fasta >> $DATADIR/yeast_seqs.RPB2.fasta
cat $DATADIR/new.tef1.fasta >> $DATADIR/yeast_seqs.TEF1.fasta
cat $DATADIR/new.cytb.fasta >> $DATADIR/yeast_seqs.CYTB.fasta


x=0
for i in SSU ITS1 5_8S ITS2 LSU
do
  sed "s/NODE_101_length_6643_cov_1305.952635/NS00000$x.1/" ../data/phylo/phylo_loci_JL201/JL201_ITSx."$i".fasta >> $DATADIR/yeast_seqs.$i.fasta
  x=$(( $x + 1 ))
  sed "s/NODE_120_length_3795_cov_610.704949/NS00000$x.1/" ../data/phylo/phylo_loci_JL221/JL221_ITSx."$i".fasta >> $DATADIR/yeast_seqs.$i.fasta
  x=$(( $x + 1 ))
done
i=$(printf "%06g" $x); x=$(( $x + 1 ))
echo ">NS$i.1 UNVERIFIED: Novel species strain JL201 RNA polymerase II largest subunit-like gene, partial sequence" >> $DATADIR/yeast_seqs.RPB1.fasta
tail -n+2 ../data/phylo/phylo_loci_JL201/RPB1.fasta >> $DATADIR/yeast_seqs.RPB1.fasta
i=$(printf "%06g" $x); x=$(( $x + 1 ))
echo ">NS$i.1 UNVERIFIED: Novel species strain JL201 RNA polymerase II second largest subunit-like gene, partial sequence" >> $DATADIR/yeast_seqs.RPB2.fasta
tail -n+2 ../data/phylo/phylo_loci_JL201/RPB2.fasta >> $DATADIR/yeast_seqs.RPB2.fasta
i=$(printf "%06g" $x); x=$(( $x + 1 ))
echo ">NS$i.1 UNVERIFIED: Novel species strain JL201 translation elongation factor 1-alpha-like gene, partial sequence" >> $DATADIR/yeast_seqs.TEF1.fasta
tail -n+2 ../data/phylo/phylo_loci_JL201/TEF1.fasta >> $DATADIR/yeast_seqs.TEF1.fasta
i=$(printf "%06g" $x); x=$(( $x + 1 ))
echo ">NS$i.1 UNVERIFIED: Novel species strain JL201 cytochrome b-like gene, partial sequence; mitochondrial" >> $DATADIR/yeast_seqs.CYTB.fasta
tail -n+2 ../data/phylo/phylo_loci_JL201/CYTB.fasta >> $DATADIR/yeast_seqs.CYTB.fasta
i=$(printf "%06g" $x); x=$(( $x + 1 ))
echo ">NS$i.1 UNVERIFIED: Novel species strain JL221 RNA polymerase II largest subunit-like gene, partial sequence" >> $DATADIR/yeast_seqs.RPB1.fasta
tail -n+2 ../data/phylo/phylo_loci_JL221/RPB1.fasta >> $DATADIR/yeast_seqs.RPB1.fasta
i=$(printf "%06g" $x); x=$(( $x + 1 ))
echo ">NS$i.1 UNVERIFIED: Novel species strain JL221 RNA polymerase II second largest subunit-like gene, partial sequence" >> $DATADIR/yeast_seqs.RPB2.fasta
tail -n+2 ../data/phylo/phylo_loci_JL221/RPB2.fasta >> $DATADIR/yeast_seqs.RPB2.fasta
i=$(printf "%06g" $x); x=$(( $x + 1 ))
echo ">NS$i.1 UNVERIFIED: Novel species strain JL221 translation elongation factor 1-alpha-like gene, partial sequence" >> $DATADIR/yeast_seqs.TEF1.fasta
tail -n+2 ../data/phylo/phylo_loci_JL221/TEF1.fasta >> $DATADIR/yeast_seqs.TEF1.fasta
i=$(printf "%06g" $x); x=$(( $x + 1 ))
echo ">NS$i.1 UNVERIFIED: Novel species strain JL221 cytochrome b-like gene, partial sequence; mitochondrial" >> $DATADIR/yeast_seqs.CYTB.fasta
tail -n+2 ../data/phylo/phylo_loci_JL221/CYTB.fasta >> $DATADIR/yeast_seqs.CYTB.fasta
