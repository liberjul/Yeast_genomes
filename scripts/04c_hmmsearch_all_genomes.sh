#!/bin/bash --login

DATADIR=../data/phylo/jgi_genomic_seqs
OUTDIR=$DATADIR/phylo_loci_genome_tree
mkdir -p $OUTDIR $OUTDIR/combined_seqs
for i in RPB1 RPB2 TEF1 CYTB LSU SSU 5_8S
	do
	echo > $OUTDIR/combined_seqs/"$i".fna
	while read GENOME;
	do
		# mkdir -p $OUTDIR/$GENOME
		# nhmmer --tblout $OUTDIR/$GENOME/"$i"_hmm.out \
		  # --cpu 4 $DATADIR/../hmms/$i.hmm $DATADIR/full_genomes/"$GENOME"_genome.fna
		python pull_hmm_hits.py -f $DATADIR/full_genomes/"$GENOME"_genome.fna \
		  -m $OUTDIR/$GENOME/"$i"_hmm.out -o $OUTDIR/$GENOME/"$i".fasta -t -e 1e-30 -s 1000 -p 1000
		if [ -f "$OUTDIR/$GENOME/"$i".fasta" ]
		then
			sed "s|^>.*$|>$GENOME|" $OUTDIR/$GENOME/"$i".fasta >> $OUTDIR/combined_seqs/"$i".fna
		fi
	done < $DATADIR/full_genomes/new_genomes.txt
	mafft --localpair --maxiterate 1000 --thread 4 $OUTDIR/combined_seqs/"$i".fna > $OUTDIR/combined_seqs/"$i".mafft.fna
done

