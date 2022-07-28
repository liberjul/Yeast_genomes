#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=02:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=24          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=32G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=assemble_JL201        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

cd $SLURM_SUBMIT_DIR

CPU=24
MEM=64
FWD=../data/filt_trimmed/JL201_S181_1P.fastq.gz
REV=../data/filt_trimmed/JL201_S181_2P.fastq.gz
WORKDIR=../data/working
OUTDIR=../data/assembly
mkdir -p $WORKDIR $OUTDIR
conda activate aaftf
BASE=$(basename ${FWD%.fastq.gz})
TRIMREAD=../data/filt_trimmed
ASMFILE=$OUTDIR/${BASE}.spades.fasta
echo $TRIMREAD/${BASE}
conda activate aaftf

spades.py --threads $CPU --mem $MEM -1 $FWD -2 $REV -o $WORKDIR/spades_$BASE
cp $WORKDIR/spades_$BASE/*.fasta $OUTDIR

conda deactivate

scontrol show job $SLURM_JOB_ID     ### write job information to output file
