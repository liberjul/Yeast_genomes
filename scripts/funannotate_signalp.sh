#!/bin/bash
#SBATCH --mem=40G
#SBATCH --partition=biodept,common,scavenger
#SBATCH --account=biodept
#SBATCH --cpus-per-task=12

# Run from project directory after activating funannotate conda environment

cd JL221
signalp6 \
--fastafile ./annotations/predict_results/*proteins.fa \
--output_dir ./annotations/annotate_misc/signalp \
--format txt \
--organism euk \
--mode fast \
--write_procs 10

cd ../JL201
signalp6 \
--fastafile ./annotations/predict_results/*proteins.fa \
--output_dir ./annotations/annotate_misc/signalp \
--format txt \
--organism euk \
--mode fast \
--write_procs 10
