#!/bin/bash
#SBATCH --mem=450M
#SBATCH --partition=biodept,common
#SBATCH --account=biodept
#SBATCH --cpus-per-task=1

# Run from project directory after activating funannotate conda environment

cd JL221
funannotate remote \
-i annotations \
-m antismash \
-e ian.medeiros@duke.edu \
-o antismash_results

cd ../JL201
funannotate remote \
-i annotations \
-m antismash \
-e ian.medeiros@duke.edu \
-o antismash_results
