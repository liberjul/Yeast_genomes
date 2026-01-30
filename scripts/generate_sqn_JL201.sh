DATADIR=../data/annotation

# GENOME=JL201
for GENOME in JL201 JL221
do
  echo $DATADIR/"$GENOME"/final/template_"$GENOME".sbt
  /mnt/c/Users/julia/Anaconda3/Library/bin/table2asn -indir $DATADIR/"$GENOME"/final \
    -outdir $DATADIR/"$GENOME"/final/ncbi \
    -t $DATADIR/"$GENOME"/final/template_"$GENOME".sbt \
    -euk -a s -Z -M n \
    -j "[organism=Microbotryomycetes sp.] [isolate=$GENOME]"
done
