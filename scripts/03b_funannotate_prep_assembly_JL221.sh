#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=00:10:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=10G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=funannotate_prep_JL221        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

## Prep assemblies
conda activate aaftf
module load GCC/10.2.0  OpenMPI/4.0.5 AUGUSTUS/3.4.0 DIAMOND hisat2/2.1.0 icc/2018.1.163-GCC-6.4.0-2.28  impi/2018.1.163 SAMtools/1.8 PASA/2.4.1 Trinity/2.8.4 kallisto/0.46.0 FASTA/36.3.8h_04-May-2020

cd /mnt/scratch/liberjul/Yeast_genomes/funannotate_working_dir
DATADIR=~/He_Lab/Yeast_genomes/data
funannotate clean -i $DATADIR/assembly/JL221/scaffolds.fasta \
  -o $DATADIR/assembly/JL221/scaffolds.clean.fasta

funannotate sort -i $DATADIR/assembly/JL221/scaffolds.clean.fasta \
  -o $DATADIR/assembly/JL221/scaffolds.sorted.fasta

funannotate mask -i $DATADIR/assembly/JL221/scaffolds.sorted.fasta \
  -o $DATADIR/assembly/JL221/scaffolds.masked.fasta

conda deactivate

scontrol show job $SLURM_JOB_ID     ### write job information to output file
