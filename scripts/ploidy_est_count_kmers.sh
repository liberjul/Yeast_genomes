#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=00:10:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=10G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=count_kmers        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########


cd ~/He_Lab/Yeast_genomes/scripts
DATADIR=~/He_Lab/Yeast_genomes/data

~/bin/bbmap/kmercountexact.sh threads=16 k=23 in=$DATADIR/raw_data/jl201_S181_R1_001.fastq.gz in2=$DATADIR/raw_data/jl201_S181_R2_001.fastq.gz khist=$DATADIR/read_mapping/JL201_kmer_hist.txt -Xmx10g peaks=$DATADIR/read_mapping/JL201_kmer_peaks.txt
~/bin/bbmap/kmercountexact.sh threads=16 k=23 in=$DATADIR/raw_data/JL221_S55_R1_001.fastq.gz in2=$DATADIR/raw_data/JL221_S55_R2_001.fastq.gz khist=$DATADIR/read_mapping/JL221_kmer_hist.txt -Xmx10g peaks=$DATADIR/read_mapping/JL221_kmer_peaks.txt

scontrol show job $SLURM_JOB_ID     ### write job information to output file
