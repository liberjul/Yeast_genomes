#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=06:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=24          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=64G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=assemble_JL201        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

cd $SLURM_SUBMIT_DIR

CPU=24
MEM=64
i=../data/raw_data/JL201_raw_R1_001.fastq.qz
WORKDIR=../data/working
OUTDIR=../data/assembly
mkdir -p $WORKDIR $OUTDIR
conda activate aaftf
BASE=$(basename ${i%.fastq.gz})
TRIMREAD=../data/filt_trimmed
INTERLEAVED=$TRIMREAD/${BASE}.fastq.gz
ASMFILE=$OUTDIR/${BASE}.spades.fasta
echo $TRIMREAD/${BASE}
conda activate aaftf

spades.py --threads $CPU --mem $MEM --12 $INTERLEAVED -o $WORKDIR/spades_$BASE
source deactivate
cp $WORKDIR/spades_$BASE/*.fasta $OUTDIR

conda deactivate

scontrol show job $SLURM_JOB_ID     ### write job information to output file
