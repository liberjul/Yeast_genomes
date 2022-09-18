#!/bin/bash --login
echo "usage: sh ./check_gRNAs_against_genome.sh <strain> <template_fasta> <grna_fasta> <output>"
# Script from https://github.com/liberjul/code_utilities/blob/main/check_grnas_v_db_usearch.py
python /mnt/c/Users/julia/Anaconda3/Lib/site-packages/code_utilities/check_grnas_v_db_usearch.py -d ../data/assembly/$1/scaffolds.fasta -t $2 -g $3 -o $4 -u /mnt/c/Users/julia/bin/usearch11.0.667_win32.exe

# python C:\Users\julia\Anaconda3\Lib\site-packages\code_utilities\check_grnas_v_db_usearch.py -d ..\data\assembly\JL201\scaffolds.fasta -t ..\data\knockout\JL201_URA3.fasta -g ..\data\knockout\JL201_URA3_gRNAs.fasta -o ..\data\knockout\JL201_URA3_gRNA_check_res -u C:\Users\julia\bin\usearch11.0.667_win32.exe
# sh ./check_gRNA_against_genome.sh JL201 ../data/knockout/JL201_URA3.fasta ../data/knockout/JL201_URA3_gRNAs.fasta ../data/knockout/JL201_URA3_gRNA_check_res.out

# python C:\Users\julia\Anaconda3\Lib\site-packages\code_utilities\check_grnas_v_db_usearch.py -d ..\data\assembly\JL221\scaffolds.fasta -t ..\data\knockout\JL221_URA3.fasta -g ..\data\knockout\JL221_URA3_gRNAs.fasta -o ..\data\knockout\JL221_URA3_gRNA_check_res -u C:\Users\julia\bin\usearch11.0.667_win32.exe
# sh ./check_gRNA_against_genome.sh JL221 ../data/knockout/JL221_URA3.fasta ../data/knockout/JL221_URA3_gRNAs.fasta ../data/knockout/JL221_URA3_gRNA_check_res.out
