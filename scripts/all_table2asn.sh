#!/bin/bash --login

ISO=JL201
DATADIR1=../data/annotation/"$ISO"_funannotate_annotate_18march25/annotate_results
DATADIR2=../data/annotation/"$ISO"_funannotate_annotate_18march25/annotate_corrected

table2asn -indir $DATADIR2 \
  -Z -V v \
  -t $DATADIR2/template_"$ISO".sbt \
  -M n \
  -j "[Organism=Microbotryomycetes sp. JL201] [strain=$ISO]" \
  -euk
  
ISO=JL221
DATADIR1=../data/annotation/"$ISO"_funannotate_annotate_18march25/annotate_results
DATADIR2=../data/annotation/"$ISO"_funannotate_annotate_18march25/annotate_corrected

table2asn -indir $DATADIR2 \
  -Z -V v \
  -t $DATADIR2/template_"$ISO".sbt \
  -M n \
  -j "[Organism=Microbotryomycetes sp. JL221] [strain=$ISO]" \
  -euk

ISO=NB124
DATADIR1=../data/annotation/"$ISO"_funannotate_annotate_18march25/annotate_results
DATADIR2=../data/annotation/"$ISO"_funannotate_annotate_18march25/annotate_corrected

table2asn -indir $DATADIR2 \
  -Z -V v \
  -t $DATADIR2/template_"$ISO".sbt \
  -M n \
  -j "[Organism=Microbotryomycetes sp. NB124-2] [strain=NB124-2]" \
  -euk

