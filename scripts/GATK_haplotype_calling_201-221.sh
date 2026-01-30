#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=01:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=4          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=10G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=count_kmers        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########


cd ~/He_Lab/Yeast_genomes/scripts
DATADIR=~/He_Lab/Yeast_genomes/data
module load Java/17.0.2

gatk --java-options "-Xmx10g" HaplotypeCaller -ERC GVCF \
  -R $DATADIR/assembly/JL201/scaffolds.masked.fasta -I $DATADIR/read_mapping/JL201.aligned.bam -O $DATADIR/read_mapping/JL201.gatkHC.vcf \
  -A DepthPerAlleleBySample -A MappingQuality -A LikelihoodRankSumTest


scontrol show job $SLURM_JOB_ID     ### write job information to output file
