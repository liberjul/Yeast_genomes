# Yeast_genomes

### This repo is scripts and data for the assembly of novel Microbotryomycete yeasts JL201 and JL221.

#### The processing pipeline involves the following steps:

1. Trim and filter reads using AAFTF trim. See the AAFTF [repo](https://github.com/stajichlab/AAFTF) for installation.
    - [01_aaftf_trim_filter_JL201.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/01_aaftf_trim_filter_JL201.sb)
    - [01_aaftf_trim_filter_JL221.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/01_aaftf_trim_filter_JL221.sb)
2. Assembly reads using SPAdes (installed as part of AAFTF).
    - [02_aaftf_assemble_JL201.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/02_aaftf_assemble_JL201.sb)
    - [02_aaftf_assemble_JL221.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/02_aaftf_assemble_JL221.sb)
3. From here, the assembly is used for three independent analyses.

  - A. Align reads to genome with [HISAT2](http://daehwankimlab.github.io/hisat2/), then use [BCFTools](http://samtools.github.io/bcftools/bcftools.html) to variant call.
      - [03a_hisat2_build_index_JL201.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/03a_hisat2_build_index_JL201.sb)
      - [03a_hisat2_build_index_JL221.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/03a_hisat2_build_index_JL221.sb)
      - [04a_align_reads_hisat2_JL201.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/04a_align_reads_hisat2_JL201.sb)
      - [04a_align_reads_hisat2_JL221.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/04a_align_reads_hisat2_JL221.sb)
      - [05a_variant_calling_JL201.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/05a_variant_calling_JL201.sb)
      - [05a_variant_calling_JL221.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/05a_variant_calling_JL221.sb)

  - B. Annotate scaffolds with [funannotate](https://github.com/nextgenusfs/funannotate/).
      - [03b_funannotate_prep_assembly_JL201.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/03b_funannotate_prep_assembly_JL201.sb)
      - [03b_funannotate_prep_assembly_JL221.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/03b_funannotate_prep_assembly_JL221.sb)
      - [04b_funannotate_predict_JL201_NP11_proteins.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/04b_funannotate_predict_JL201_NP11_proteins.sb)
      - [04b_funannotate_predict_JL221_NP11_proteins.sb](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/04b_funannotate_predict_JL221_NP11_proteins.sb)
      - [funannotate_signalp.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/funannotate_signalp.sh)
      - [NetGPI](https://services.healthtech.dtu.dk/service.php?NetGPI-1.1) annotations run on webserver with SignalP results
      - [EffectorP.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/EffectorP.sh) annotations run on SignalP results
      - [antiSMASH](https://fungismash.secondarymetabolites.org/#!/start) annotations run on webserver with [funannotate_remote.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/funannotate_remote.sh)
      - [interproscan.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/interproscan.sh)
      - [funannotate_annotate_round1.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/funannotate_annotate_round1.sh)
      - [funannotate_annotate_round2.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/funannotate_annotate_round2.sh)
  - C. Pull out loci for phylogenetic analysis using [HMMER](http://hmmer.org/).
      - [03c_hmmbuild.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/03c_hmmbuild.sh)
      - [04c_hmmsearch_JL201.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/04c_hmmsearch_JL201.sh)
      - [04c_hmmsearch_JL221.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/04c_hmmsearch_JL221.sh)
      - [05c_split_fastas_by_reg.py](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/05c_split_fastas_by_reg.py)
      - [06c_extract_rRNA_regions.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/06c_extract_rRNA_regions.sh)
      - [07c_add_JL201_JL221_to_alignments.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/07c_add_JL201_JL221_to_alignments.sh)
      - [08c_align_w_mafft.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/08c_align_w_mafft.sh)
      - [09c_clip_w_clipkit.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/09c_clip_w_clipkit.sh)
      - [10c_format_fasta_headers.py](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/10c_format_fasta_headers.py)
      - [11c_concat_regions.sh](https://github.com/liberjul/Yeast_genomes/blob/main/scripts/11c_concat_regions.sh)
