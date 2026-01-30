#!/bin/bash
#SBATCH --mem=60G
#SBATCH --partition=biodept,common,scavenger
#SBATCH --account=biodept
#SBATCH --cpus-per-task=10
#SBATCH --array=1-2

# Run from project directory after activating funannotate conda environment

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

# First, we remove EC number notes from genes with product "hypothetical protein" per NCBI submission requirements
cd annotations
cp ./annotate_results/${species}_${genome}.discrepency.report.txt ./${species}_${genome}_original.discrepency.report.txt

cd annotate_misc
cp eggnog.emapper.annotations eggnog.emapper.annotations.hypothetical-not-corrected

cat ../annotate_results/*.gff3 | grep "product=hypothetical protein;" | grep -E --only-matching "FUN_[0-9]{6}" | awk -F'_' '{print $2}' | sort | uniq > ./hypothetical_protein_remove_EC_number.txt
for j in $(cat hypothetical_protein_remove_EC_number.txt)
do
sed -E -i "/FUN_$j-T1/s/[0-9]\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,4}/-/g" eggnog.emapper.annotations
done
sed -E -i "s/-(,-){1,10}/-/g" eggnog.emapper.annotations

cd ../..

# Now we re-run funannotate with the cleaned eggnog results, the final .sbt file and locus tag, and any fixed gene-product combinations
funannotate annotate \
-i annotations \
--species $species \
--strain $genome \
--iprscan ./annotations/interproscan_results/${genome}.xml \
--signalp ./annotations/annotate_misc/signalp/prediction_results.txt \
--netgpi ./annotations/annotate_misc/NetGPI/output_protein_type.txt \
--effectorp ./annotations/annotate_misc/effectors.txt \
--deeploc ./deeploc_out/results_20221020-${deep}.csv \
--busco_db basidiomycota \
--rename $locus_tag \
--sbt ./template_${genome}.sbt \
--fix ./annotations/${genome}_GenesToFix.txt \
--renumber_antismash \
--no-progress \
--cpus 10
