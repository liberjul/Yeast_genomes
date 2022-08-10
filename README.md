# Yeast_genomes

### This repo is scripts and data for the assembly of novel Microbotryomycete yeast JL201 and JL221.

#### The processing pipeline involves the following steps:

1. Trim and filter reads using AAFTF trim.
2. Assembly reads using SPAdes.
3. From here, the assembly is used for three independent analyses.
  A. Align reads to genome with HISAT2, then use BCFTools to variant call.
  B. Annotate scaffolds with funannotate.
  C. Pull out loci for phylogenetic analysis using HMMER.
