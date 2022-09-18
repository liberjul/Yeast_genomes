#!/bin/bash --login


DATADIR=../data/annotation

blastp -query $DATADIR/ura3_NP11.fasta -db $DATADIR/JL221/Microbotryomycetes_sp_JL221_prot > $DATADIR/JL221/blastp_ura3.out

blastp -query $DATADIR/ura3_NP11.fasta -db $DATADIR/JL201/Microbotryomycetes_sp_JL201_prot > $DATADIR/JL201/blastp_ura3.out
