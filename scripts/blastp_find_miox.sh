#!/bin/bash --login


DATADIR=../data/annotation

blastp -query $DATADIR/UstMay_MIOX.fasta -db $DATADIR/JL221/Microbotryomycetes_sp_JL221_prot > $DATADIR/JL221/blastp_MIOX.out

blastp -query $DATADIR/UstMay_MIOX.fasta -db $DATADIR/JL201/Microbotryomycetes_sp_JL201_prot > $DATADIR/JL201/blastp_MIOX.out
