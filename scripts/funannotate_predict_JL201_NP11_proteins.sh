#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=32G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=funannotate_predict_JL201_NP11_prot        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

## Prep assemblies

export QUARRY_PATH=~/bin/CodingQuarry_v2.0/QuarryFiles
export ZOE=~/bin/SNAP/Zoe
export AUGUSTUS_CONFIG_PATH=~/bin/anaconda3/envs/aaftf/config/
#FUNANNOTATE_DB and GENEMARK_PATH are specified in ~/.bashrc
#GlimmerHMM and SNAP are added to the path in ~/.bashrc
conda activate aaftf
module load GCC/10.2.0  OpenMPI/4.0.5 AUGUSTUS/3.4.0 DIAMOND hisat2/2.1.0 icc/2018.1.163-GCC-6.4.0-2.28  impi/2018.1.163 SAMtools/1.8 PASA/2.4.1 Trinity/2.8.4 kallisto/0.46.0 FASTA/36.3.8h_04-May-2020
export AUGUSTUS_CONFIG_PATH=~/bin/anaconda3/envs/aaftf/config/

cd /mnt/scratch/liberjul/Yeast_genomes/funannotate_working_dir
DATADIR=~/He_Lab/Yeast_genomes/data

funannotate predict -i $DATADIR/assembly/JL201/scaffolds.masked.fasta \
  -o $DATADIR/annotation/JL201 \
  -s "Microbotryomycetes sp" --strain JL201 \
  --protein_evidence $DATADIR/annotation/Rhoto1_GeneCatalog_proteins_20141208.aa.fasta \
  --busco_db basidiomycota \
  --cpus 16

conda deactivate

scontrol show job $SLURM_JOB_ID     ### write job information to output file
