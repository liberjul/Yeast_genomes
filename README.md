# Yeast_genomes

### This repo is scripts and data for the assembly of novel Microbotryomycete yeasts JL201 and JL221.

#### The processing pipeline involves the following steps:

1. Trim and filter reads using AAFTF trim. See the AAFTF [repo](https://github.com/stajichlab/AAFTF) for installation.
2. Assembly reads using SPAdes (installed as part of AAFTF).
3. From here, the assembly is used for three independent analyses.

  - A. Align reads to genome with [HISAT2](http://daehwankimlab.github.io/hisat2/), then use [BCFTools](http://samtools.github.io/bcftools/bcftools.html) to variant call.

  - B. Annotate scaffolds with [funannotate](https://github.com/nextgenusfs/funannotate/).

  - C. Pull out loci for phylogenetic analysis using [HMMER](http://hmmer.org/).
