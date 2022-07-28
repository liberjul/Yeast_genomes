#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=02:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=8          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=16G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=trim_filter_reads_JL201        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

cd $SLURM_SUBMIT_DIR

MEM=16 # 16gb
CPU=8
FWD=../data/raw_data/jl201_S181_R1_001.fastq.gz
REV=../data/raw_data/jl201_S181_R2_001.fastq.gz
BASE=JL201_S181
READSDIR=../data/raw_data
TRIMREAD=../data/filt_trimmed

echo $BASE
echo $READSDIR
conda activate aaftf

AAFTF trim --method bbduk --memory $MEM -c $CPU \
 --left $FWD --right $REV\
  -o $TRIMREAD/${BASE}
conda deactivate

scontrol show job $SLURM_JOB_ID
