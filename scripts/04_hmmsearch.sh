#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=00:10:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=24          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=8G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=hmmsearch        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

cd $SLURM_SUBMIT_DIR

DATADIR=../data/phylo
conda activate aaftf

for i in SSU ITS1 5_8S ITS2 LSU RPB1 RPB2 TEF1 CYTB
do
  nhmmer --tblout ../data/assembly/JL201/phylo_loci/"$i"_hmm.out --cpu 24 $DATADIR/hmms/$i.hmm ../data/assembly/JL201/scaffolds.fasta
  python pull_hmm_hits.py -f ../data/assembly/JL201/scaffolds.fasta \
    -m ../data/assembly/JL201/phylo_loci/"$i"_hmm.out -o ../data/assembly/JL201/phylo_loci/"$i".fasta -t
done

conda deactivate

scontrol show job $SLURM_JOB_ID     ### write job information to output file
