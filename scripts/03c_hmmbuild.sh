#!/bin/bash --login

DATADIR=../data/phylo
# conda activate aaftf

mkdir -p $DATADIR/hmms/
for i in RPB1 RPB2 TEF1 CYTB
do
  hmmbuild $DATADIR/hmms/$i.hmm $DATADIR/aligned_seqs/yeast_seqs.$i.mafft.fasta
done

# conda deactivate

# scontrol show job $SLURM_JOB_ID     ### write job information to output file
