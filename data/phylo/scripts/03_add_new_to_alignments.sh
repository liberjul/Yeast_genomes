#!/bin/bash --login

DATADIR=../unaligned_seqs

python 01_split_fastas_by_reg.py --new

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

### JL201
STRAIN=JL201
x=0
for i in SSU ITS1 5_8S ITS2 LSU
do
  sed "s/NODE_101_length_6643_cov_1305.952635/NS00000$x.1/" ../phylo_loci_$STRAIN/"$STRAIN"_ITSx."$i".fasta >> $DATADIR/yeast_seqs.$i.fasta
  x=$(( $x + 1 ))
done
echo ">NS000006.1 UNVERIFIED: Novel species strain "$STRAIN" RNA polymerase II largest subunit-like gene, partial sequence" >> $DATADIR/yeast_seqs.RPB1.fasta
tail -n+2 ../phylo_loci_"$STRAIN"/RPB1.fasta >> $DATADIR/yeast_seqs.RPB1.fasta
echo ">NS000007.1 UNVERIFIED: Novel species strain "$STRAIN" RNA polymerase II second largest subunit-like gene, partial sequence" >> $DATADIR/yeast_seqs.RPB2.fasta
tail -n+2 ../phylo_loci_"$STRAIN"/RPB2.fasta >> $DATADIR/yeast_seqs.RPB2.fasta
echo ">NS000008.1 UNVERIFIED: Novel species strain "$STRAIN" translation elongation factor 1-alpha-like gene, partial sequence" >> $DATADIR/yeast_seqs.TEF1.fasta
tail -n+2 ../phylo_loci_"$STRAIN"/TEF1.fasta >> $DATADIR/yeast_seqs.TEF1.fasta
echo ">NS000009.1 UNVERIFIED: Novel species strain "$STRAIN" cytochrome b-like gene, partial sequence; mitochondrial" >> $DATADIR/yeast_seqs.CYTB.fasta
tail -n+2 ../phylo_loci_"$STRAIN"/CYTB.fasta >> $DATADIR/yeast_seqs.CYTB.fasta
