#!/bin/bash --login

DATADIR=../data/annotation

awk '/>FUN_003387-T1/{f=1; c=0} f; />/ && ++c==2{f=0}' $DATADIR/JL201/Microbotryomycetes_sp_JL201.proteins.fa | head -n-1 > $DATADIR/JL201/JL201_mep1_homolog1.fasta
awk '/>FUN_004464-T1/{f=1; c=0} f; />/ && ++c==2{f=0}' $DATADIR/JL201/Microbotryomycetes_sp_JL201.proteins.fa | head -n-1 > $DATADIR/JL201/JL201_mep1_homolog2.fasta
blastp -query $DATADIR/JL201/JL201_mep1_homolog1.fasta -db $DATADIR/JL221/Microbotryomycetes_sp_JL221_prot > $DATADIR/JL201/JL201_mep1_homolog1_blast.out
blastp -query $DATADIR/JL201/JL201_mep1_homolog2.fasta -db $DATADIR/JL221/Microbotryomycetes_sp_JL221_prot > $DATADIR/JL201/JL201_mep1_homolog2_blast.out

awk '/>FUN_003375-T1/{f=1; c=0} f; />/ && ++c==2{f=0}' $DATADIR/JL221/Microbotryomycetes_sp_JL221.proteins.fa | head -n-1 > $DATADIR/JL221/JL221_mep1_homolog1.fasta
awk '/>FUN_003945-T1/{f=1; c=0} f; />/ && ++c==2{f=0}' $DATADIR/JL221/Microbotryomycetes_sp_JL221.proteins.fa | head -n-1 > $DATADIR/JL221/JL221_mep1_homolog2.fasta
blastp -query $DATADIR/JL221/JL221_mep1_homolog1.fasta -db $DATADIR/JL201/Microbotryomycetes_sp_JL201_prot > $DATADIR/JL221/JL221_mep1_homolog1_blast.out
blastp -query $DATADIR/JL221/JL221_mep1_homolog2.fasta -db $DATADIR/JL201/Microbotryomycetes_sp_JL201_prot > $DATADIR/JL221/JL221_mep1_homolog2_blast.out
