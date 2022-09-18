#!/bin/bash --login


DATADIR=../data/annotation

makeblastdb -in $DATADIR/JL201/Microbotryomycetes_sp_JL201.proteins.fa \
  -input_type fasta -dbtype prot -title JL201_prot -out $DATADIR/JL201/Microbotryomycetes_sp_JL201_prot

makeblastdb -in $DATADIR/JL221/Microbotryomycetes_sp_JL221.proteins.fa \
  -input_type fasta -dbtype prot -title JL221_prot -out $DATADIR/JL221/Microbotryomycetes_sp_JL221_prot
