#!/bin/bash --login


DATADIR=../data/annotation

blastp -query $DATADIR/S288C_YGR121C_MEP1_protein.fsa -db $DATADIR/JL221/Microbotryomycetes_sp_JL221_prot > $DATADIR/JL221/blastp_mep1.out

blastp -query $DATADIR/S288C_YGR121C_MEP1_protein.fsa -db $DATADIR/JL201/Microbotryomycetes_sp_JL201_prot > $DATADIR/JL201/blastp_mep1.out
