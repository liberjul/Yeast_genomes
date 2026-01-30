#!/bin/bash
#SBATCH --mem=500M
#SBATCH --partition=biodept,common,scavenger
#SBATCH --account=biodept
#SBATCH --cpus-per-task=4

cd /hpc/group/bio1/ian/envs/funannotate/EffectorP-3.0/

genome='JL221'
python EffectorP.py \
-i /hpc/group/bio1/ian/Yeast_genomes/$genome/annotations/annotate_misc/signalp/processed_entries.fasta \
-o /hpc/group/bio1/ian/Yeast_genomes/$genome/annotations/annotate_misc/effectors.txt

genome='JL201'
python EffectorP.py \
-i /hpc/group/bio1/ian/Yeast_genomes/$genome/annotations/annotate_misc/signalp/processed_entries.fasta \
-o /hpc/group/bio1/ian/Yeast_genomes/$genome/annotations/annotate_misc/effectors.txt
