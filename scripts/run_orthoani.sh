#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=2:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=4          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=10G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name OrthoANI_all          # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

cd ${SLURM_SUBMIT_DIR}

module load NCBI-BLAST/2.12.0-rhel8

DATADIR=~/Yeast_genomes/data
BLASTDIR=/hpc/dctrl/jal138/software/ncbi-blast-2.16.0+/bin/

java -jar ~/bin/OAT_cmd.jar -blastplus_dir $BLASTDIR \
  -fasta1 $DATADIR/assembly/Aimania_JL201.fsa \
  -fasta2 $DATADIR/assembly/Aimania_JL221.fsa > $DATADIR/ani/JL201-JL221_orthoani.txt
java -jar ~/bin/OAT_cmd.jar -blastplus_dir $BLASTDIR \
  -fasta1 $DATADIR/assembly/Aimania_JL201.fsa \
  -fasta2 $DATADIR/assembly/Aimania_NB124.fsa > $DATADIR/ani/JL201-NB124_orthoani.txt
java -jar ~/bin/OAT_cmd.jar -blastplus_dir $BLASTDIR \
  -fasta1 $DATADIR/assembly/Aimania_NB124.fsa \
  -fasta2 $DATADIR/assembly/Aimania_JL221.fsa > $DATADIR/ani/NB124-JL221_orthoani.txt

java -jar ~/bin/OAT_cmd.jar -blastplus_dir $BLASTDIR \
  -fasta2 $DATADIR/assembly/Aimania_JL201.fsa \
  -fasta1 $DATADIR/assembly/Aimania_JL221.fsa > $DATADIR/ani/JL221-JL201_orthoani.txt
java -jar ~/bin/OAT_cmd.jar -blastplus_dir $BLASTDIR \
  -fasta2 $DATADIR/assembly/Aimania_JL201.fsa \
  -fasta1 $DATADIR/assembly/Aimania_NB124.fsa > $DATADIR/ani/NB124-JL201_orthoani.txt
java -jar ~/bin/OAT_cmd.jar -blastplus_dir $BLASTDIR \
  -fasta2 $DATADIR/assembly/Aimania_NB124.fsa \
  -fasta1 $DATADIR/assembly/Aimania_JL221.fsa > $DATADIR/ani/JL221-NB124_orthoani.txt

scontrol show job ${SLURM_JOB_ID}
