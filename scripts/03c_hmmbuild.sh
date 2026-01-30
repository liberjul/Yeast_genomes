#!/bin/bash --login

DATADIR=../data/phylo
# conda activate aaftf

mkdir -p $DATADIR/hmms/
for i in RPB1 RPB2 TEF1 CYTB 5_8S SSU LSU
do
  echo $i
  hmmbuild $DATADIR/hmms/$i.hmm $DATADIR/aligned_seqs/class.yeast_seqs.$i.dedup_types_no_miss.fasta
done

# conda deactivate

# scontrol show job $SLURM_JOB_ID     ### write job information to output file
