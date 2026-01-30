#!/bin/bash --login

ISO=NB124
DATADIR1=../data/annotation/"$ISO"_funannotate_annotate_18march25/annotate_results
DATADIR2=../data/annotation/"$ISO"_funannotate_annotate_18march25/annotate_corrected

python ./correct_gene_names_tbl.py \
  -t $DATADIR1/Aimania_"$ISO".tbl \
  -l ACM66B \
  -o $DATADIR2/Aimania_"$ISO".tbl \
  -c $DATADIR2/contig_conversion.txt


# grep -vP "CFMR|REFERENCE" $DATADIR2/Aimania_"$ISO".tbl > $DATADIR2/Aimania_"$ISO".edit.tbl 
# mv $DATADIR2/Aimania_"$ISO".edit.tbl $DATADIR2/Aimania_"$ISO".tbl
# Note: manual fix of contamination.
table2asn -indir $DATADIR2 \
  -Z -V v \
  -t $DATADIR2/template_"$ISO".sbt \
  -M n \
  -j "[Organism=Microbotryomycetes sp.] [strain=NB124-2]" \
  -euk

