#!/bin/bash
#SBATCH --mem=100G
#SBATCH --partition=biodept,common,scavenger
#SBATCH --account=biodept
#SBATCH --cpus-per-task=10
#SBATCH --array=1-2

# run from project directory after activating funannotate conda environment

# This is the first annotation round. We will run "funannotate annotate" a second time to create 
# NCBI-ready outputs with the final sbt file and locus tag prefix, as well as making corrections 
# to gene names and removing EC numbers from hypothetical proteins.

export EGGNOG_DATA_DIR=/hpc/group/bio1/ian/envs/funannotate/eggnog-mapper-data
export PATH=/hpc/group/bio1/ian/envs/funannotate/share/phobius/:$PATH

# BioProject PRJNA892096

if [ ${SLURM_ARRAY_TASK_ID} == 1 ]; then
	genome='JL201'
	biosample='SAMN31367384'
	locus_tag='OIV83'
	species='Microbotryomycetes_sp_JL201'
	deep='025829'
elif [ ${SLURM_ARRAY_TASK_ID} == 2 ]; then
	genome='JL221'
	biosample='SAMN31367475'
	locus_tag='OIO90'
	species='Microbotryomycetes_sp_JL221'
	deep='032757'
fi

cd $genome

funannotate annotate \
-i annotations \
--species $species \
--strain $genome \
--iprscan ./annotations/interproscan_results/${genome}.xml \
--signalp ./annotations/annotate_misc/signalp/prediction_results.txt \
--netgpi ./annotations/annotate_misc/NetGPI/output_protein_type.txt \
--deeploc ./deeploc_out/results_20221020-${deep}.csv \
--busco_db basidiomycota \
--renumber_antismash \
--no-progress \
--cpus 10
