#!/bin/bash
#SBATCH --mem=20G
#SBATCH --partition=biodept,common,scavenger
#SBATCH --account=biodept
#SBATCH --cpus-per-task=14

# Run from project directory after activating interproscan conda environment

cd /work/idm7/iprscan
mkdir JL201
mkdir JL221

cd /hpc/group/bio1/ian/Yeast_genomes
cd ./JL201/annotations
mkdir interproscan_results

interproscan.sh \
-i ./predict_results/*.proteins.fa --seqtype p \
--disable-precalc --excl-applications SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE \
--iprlookup --goterms --pathways \
--output-file-base ./interproscan_results/JL201 \
--tempdir /work/idm7/iprscan/JL201 \
-cpu 12 

cd /hpc/group/bio1/ian/Yeast_genomes
cd ./JL221/annotations
mkdir interproscan_results

interproscan.sh \
-i ./predict_results/*.proteins.fa --seqtype p \
--disable-precalc --excl-applications SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE \
--iprlookup --goterms --pathways \
--output-file-base ./interproscan_results/JL221 \
--tempdir /work/idm7/iprscan/JL221 \
-cpu 12 
