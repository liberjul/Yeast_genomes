#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --partition=common
#SBATCH --cpus-per-task 16

top_dir=$PWD
SCRIPTSDIR=~/Yeast_genomes/scripts
echo "$top_dir"
for genome in *;
do
	if [[ "$genome" != "busco_downloads" ]]
	then
		for gene in $genome/run_microbotryomycetes_odb12/busco_sequences/single_copy_busco_sequences/*.fna;
		do
			echo $gene >> $top_dir/all_buscos.txt
		done
	fi
done
cd $top_dir
sed -i "s|^.*/run_microbotryomycetes_odb12/busco_sequences/single_copy_busco_sequences/||" all_buscos.txt
sort all_buscos.txt | uniq > all_buscos_uniq.txt

mkdir $top_dir/concatenated_buscos

for genome in *;
do
	if [[ "$genome" != "busco_downloads" ]]
	then
		mkdir -p $genome/seqs_header_fixed
		for gene in $genome/run_microbotryomycetes_odb12/busco_sequences/single_copy_busco_sequences/*.fna;
		do
			sed "s/^>.*$/>$genome/" $gene > $genome/seqs_header_fixed/$(basename -- $gene)
		done
	fi
done
# cd $top_dir

while read p;
do
	echo "concatenating $p"
	cat $top_dir/*/seqs_header_fixed/$p > $top_dir/concatenated_buscos/$p
done < all_buscos_uniq.txt

cd $top_dir/concatenated_buscos/

x=0
echo "ArrayTaskID	File" > $SCRIPTSDIR/array_config_buscos.txt
for file in $top_dir/concatenated_buscos/*.fna;
        do
		x=$((x+1))
		echo $x"	"$(basename -- $file) >> $SCRIPTSDIR/array_config_buscos.txt
        #echo "mafft --localpair --maxiterate 1000 --thread 10 ${file} > ${file//.fna}.mafft.fna" >> mafft_commands.txt;
done

# source /hpc/home/kst24/miniconda3/bin/activate iqtree

# for file in *;
        # do
        # echo "iqtree -s ${file//.fna}.mafft.fna  -alrt 1000 -bb 1000 -nt AUTO -pre ${file//.fna}_iqtree" >> iqtree_commands.txt;
# done