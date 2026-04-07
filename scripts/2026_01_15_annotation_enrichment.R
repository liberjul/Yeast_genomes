library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(stringr)
library(plotly)
library(readr)
library(DESeq2)
library(data.table)
library(purrr)
library(topGO)
library(rstatix)
# library(ALL)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
### Example data from RNASEQ
geneID2GO <- readMappings("../../RNASeq/2025_05_29_JL_PC_NW/data/metadata/AurpulNBB1_geneid2go.map")
GO2geneID <- inverseList(geneID2GO)
GO2geneID
str(head(GO2geneID))
geneNames <- names(geneID2GO)
head(geneNames)


# GOdata_up <- new("topGOdata", ontology = "MF", allGenes = geneListAp_up,
#                  annot = annFUN.gene2GO, gene2GO = geneID2GO)
# resultFisher_up <- runTest(GOdata_up, algorithm = "classic", statistic = "fisher")
# 
# geneListAp <- factor(as.integer(unique(gene_to_gene$NBB_id) %in% sig_gene_list_NBB))
# names(geneListAp) <- unique(gene_to_gene$NBB_id)

###
IPS_cols <- c("transcript", "MD5", "length", "source", "accession", "description", "start", "stop", "score", "status", "date", "IP_accession", "IP_description")

JL201_annot <- read_delim("../data/annotation/JL201/final/JL201_assembly_v2_InterProScan.tsv",
                          col_names = IPS_cols)
JL221_annot <- read_delim("../data/annotation/JL221/JL221_assembly_v2_BRAKER3_prot_InterProScan.tsv",
                          col_names = IPS_cols)
NB124_annot <- read_delim("../data/annotation/NB124-2/NB124-2_assembly_v1_BRAKER3_InterProScan.tsv",
                          col_names = IPS_cols)

JL201_annot %>%
  dplyr::filter(source == "Pfam")

JL201_annot %>%
  dplyr::group_by(source, accession) %>%
  summarize(count_JL201 = n()) -> term_counts_JL201

JL221_annot %>%
  dplyr::group_by(source, accession) %>%
  summarize(count_JL221 = n()) -> term_counts_JL221

dplyr::full_join(term_counts_JL201,
                 term_counts_JL221, 
                 by = c("source", "accession")) %>%
  dplyr::filter(source == "Pfam") %>%
  replace_na(replace = list(count_JL201 = 0,
                            count_JL221 = 0)) -> JL201_v_JL221_df

JL201_v_JL221_df[,3:4] %>%
  as.matrix() %>%
  as.table() -> xtab
dimnames(xtab) <- list(accession = JL201_v_JL221_df$accession,
                       strain = c("JL201", "JL221"))

xtab
row_wise_fisher_test(xtab) -> res
res %>%
  filter(p.adj < 0.1)

JL201_v_JL221_df %>%
  filter(accession == "PF00665")
### Let's do all of them now
aimania_annot_cts <- list()
aimania_annot_cts[["JL201"]] <- read_delim("../data/annotation/JL201/final/JL201_assembly_v2_InterProScan.tsv",
                                           col_names = IPS_cols) %>%
  dplyr::group_by(source, accession) %>%
  summarize(count = n())

aimania_annot_cts[["JL221"]] <- read_delim("../data/annotation/JL221/JL221_assembly_v2_BRAKER3_prot_InterProScan.tsv",
                                           col_names = IPS_cols) %>%
  dplyr::group_by(source, accession) %>%
  summarize(count = n())

aimania_annot_cts[["NB124"]] <- read_delim("../data/annotation/NB124-2/NB124-2_assembly_v1_BRAKER3_InterProScan.tsv",
                                           col_names = IPS_cols) %>%
  dplyr::group_by(source, accession) %>%
  summarize(count = n())
mbm_annot_cts <- list()
fp <- "../data/annotation/MBM_proteins/101_ InterProScan/"
files <- list.files(fp)
files
for (i in files){
  strain <- str_split_i(i, ".nostop", 1)
  mbm_annot_cts[[strain]] <- read_delim(paste0(fp, i), delim = "\t",
                                        col_names = IPS_cols) %>%
    dplyr::group_by(source, accession) %>%
    summarize(count = n())
}

for (i in names(mbm_annot_cts)){
  mbm_annot_cts[[i]]$strain <- i
}
list_rbind(mbm_annot_cts) %>%
  dplyr::group_by(source, accession) %>%
  summarize(total_mbm = sum(count)) -> mbm_totals

for (i in names(aimania_annot_cts)){
  aimania_annot_cts[[i]]$strain <- i
}
list_rbind(aimania_annot_cts) %>%
  dplyr::group_by(source, accession) %>%
  summarize(total_aim = sum(count)) -> aimania_totals

full_join(aimania_totals, mbm_totals,
          by = c("source", "accession")) %>%
  # dplyr::filter(source == "Pfam") %>%
  replace_na(replace = list(total_aim = 0,
                            total_mbm = 0)) -> aim_v_mbm_df

aim_v_mbm_df[,3:4] %>%
  as.matrix() %>%
  as.table() -> xtab
dimnames(xtab) <- list(accession = aim_v_mbm_df$accession,
                       group = c("Aimania", "MBM"))

xtab
row_wise_fisher_test(xtab) -> res

res %>%
  write_delim("2026_01_15_fisher_test_all_anots_aim_v_mbm.txt")

res <- read_delim("2026_01_15_fisher_test_all_anots_aim_v_mbm.txt")

res %>%
  filter(p.adj < 0.05) -> sig_res

colSums(xtab)

aim_v_mbm_df %>%
  filter(accession %in% sig_res$group) %>%
  view()

fp <- "../data/annotation/MBM_proteins/101_ InterProScan/"
all_IPS_files <- c(paste0(fp, list.files(fp)),
                   "../data/annotation/JL201/final/JL201_assembly_v2_InterProScan.tsv",
                   "../data/annotation/JL221/JL221_assembly_v2_BRAKER3_prot_InterProScan.tsv",
                   "../data/annotation/NB124-2/NB124-2_assembly_v1_BRAKER3_InterProScan.tsv")
all_IPS_files 
dscr_list <- list()
for (i in all_IPS_files){
  dscr_list[[i]] <- read_delim(i, delim = "\t",
                               col_names = IPS_cols) %>%
    dplyr::group_by(source, accession, description) %>%
    summarize(count = n())
}

dscr_list %>%
  list_rbind() %>%
  group_by(source, accession, description) %>%
  summarize(total = sum(count)) -> total_cts
total_cts

total_cts %>%
  filter(accession %in% sig_res$group) %>%
  view()

total_cts %>%
  filter(accession == "")

aim_v_mbm_df %>%
  mutate(prop_aim = total_aim/sum(total_aim),
         prop_mbm = total_mbm/sum(total_mbm),
         log2_ratio_aim_to_mbm = log2(prop_aim / prop_mbm)) -> aim_v_mbm_df_prop

aim_v_mbm_df_prop %>%
  filter(accession %in% sig_res$group) %>%
  left_join(total_cts, by = c("source", "accession")) %>%
  left_join(sig_res %>% dplyr::rename(accession = group),
            by = "accession") %>%
  ungroup() %>%
  dplyr::select(-c(total_aim:prop_mbm, source)) %>%
  dplyr::arrange(p.adj)-> enrich_df

view(enrich_df)


enrich_df %>%
  mutate(dir = case_when(log2_ratio_aim_to_mbm < 0 ~ "Depleted in Aimea",
                         .default = "Enriched in Aimea"),
         description = reorder(description, p.adj)) %>%
  filter(!description %in% c("consensus disorder prediction", "-")) %>%
  ggplot(aes(x = -log10(p.adj), y = description,
             size = total, fill = log2_ratio_aim_to_mbm)) +
  facet_wrap(~dir, scales = "free_y") +
  scale_fill_gradient2() +
  geom_point(shape = 21, color = "black") +
  labs(x = "-log10(P.adjusted)",
       y = "Annotation term",
       size = "Gene count",
       fill = "Log2 enrichment in\nAimea vs other\nMicrobotryomycetes") -> g
g
ggsave("../figs/2026_01_16_genome_annot_term_enrich_genus_class.svg", g, width = 14, height = 6)
### Lets compare within the 3 species

aimania_annot_cts %>%
  list_rbind() -> aimania_annot_df

intragen_res <- list()
for (i in c("JL201", "JL221", "NB124")){
  aimania_annot_df %>%
    filter(strain != i) %>%
    group_by(source, accession) %>%
    summarize(other = sum(count)) %>%
    left_join(aimania_annot_df %>%
                filter(strain == i),
              by = c("source", "accession")) %>%
    replace_na(replace = list(other = 0,
                              count = 0)) %>%
    mutate(prop_other = other/sum(other),
           prop_count = count/sum(count),
           log2_ratio = log2(prop_count/prop_other)) -> temp_df
  temp_df[,3:4] %>%
    as.matrix() %>%
    as.table() -> xtab
  dimnames(xtab) <- list(accession = temp_df$accession,
                         group = c(i, "other Aimea"))
  
  row_wise_fisher_test(xtab) %>%
    filter(p.adj < 0.05) -> intragen_res[[i]]
  
  intragen_res[[i]] %>%
    dplyr::rename(accession = group) %>%
    left_join(total_cts %>% dplyr::select(source, accession, description),
              by = "accession") %>%
    left_join(temp_df, by = c("source", "accession")) %>%
    mutate(strain = i) -> intragen_res[[i]]
}
for (i in names(intragen_res)){
  intragen_res[[i]]
  
}
intragen_res[[i]]
intragen_res %>%
  list_rbind() %>%
  filter(!description %in% c("consensus disorder prediction", "-")) %>%
  mutate(dir = case_when(log2_ratio < 0 ~ "Depleted",
                         .default = "Enriched"),
         strain_name = factor(case_when(strain == "JL201" ~ "Aimea erigeronia JL201",
                                        strain == "JL221" ~ "Aimea cardamina JL221",
                                        strain == "NB124" ~ "Aimea sorghi NB124-2"),
                              levels = c("Aimea erigeronia JL201",
                                         "Aimea cardamina JL221",
                                         "Aimea sorghi NB124-2")),
         description = reorder(description, p.adj)) %>%
  ggplot(aes(x = -log10(p.adj), y = description,
             size = other+count, fill = log2_ratio, shape = dir),
  ) +
  facet_grid(~strain_name, scales = "free_y") +
  scale_fill_gradient2() +
  geom_point(color = "black") +
  scale_shape_manual(breaks = c("Enriched", "Depleted"),
                     values = c(21:22)) +
  labs(x = "-log10(P.adjusted)",
       y = "Annotation term",
       size = "Gene count",
       fill = "Log2 ratio in\nisolate vs other\nAimea spp.",
       shape = "Direction") -> g
g
ggsave("../figs/2026_01_16_genome_annot_term_enrich_intragenus.svg", g, width = 14, height = 4.5)


aimania_annot_df %>%
  filter(accession == "PF00665") # Integrase

aimania_annot_df %>%
  filter(accession == "PF00078") # RT1

aimania_annot_df %>%
  filter(accession == "PF07727") # RT2

aimania_annot_df %>%
  filter(accession == "PF005001") # 

JL221_annot %>%
  filter(accession == "PF00665") %>%
  view()
### Location of retrotransposon-like sequences
retro_annots <- c("Retrotransposon gag protein",
                  "Ribonuclease H-like",
                  "Transposon Ty3-G Gag-Pol polyprotein",
                  "RNA polymerase Rpb1 C-terminal repeat",
                  "RNase_HI_RT_Ty1",
                  "Transposon Tf2-6 polyprotein-like Protein",
                  "Retrovirus-related Pol polyprotein from transposon 17.6-like protein",
                  "RNase H-like domain found in reverse transcriptase",
                  "RT_LTR",
                  "Reverse transcriptase (RNA-dependent DNA polymerase)",
                  "Integrase zinc binding domain",
                  "TRANSPOSON TY3-I GAG-POL POLYPROTEIN",
                  "RNase_HI_RT_Ty3",
                  "Retrovirus zinc finger-like domains",
                  "GAG-POL-RELATED RETROTRANSPOSON",
                  "Retrovirus-related Pol polyprotein from transposon gypsy-like Protein",
                  "Integrase catalytic domain profile.",
                  "Retrovirus-related Pol polyprotein from transposon 297-like Protein",
                  "DNA/RNA polymerases",
                  "Reverse transcriptase (RT) catalytic domain profile.",
                  "Integrase core domain")

#JL201

JL201_annot %>%
  filter(description %in% retro_annots) %>%
  pull(var = "transcript") %>%
  unique() -> JL201_TE_assoc

JL201_gtf <- read_delim("../data/annotation/JL201/final/JL201_assembly_v2_BRAKER3_RNA_products.gtf",
                        col_names = c("chrom_old", "source", "locus_type", "start", "end", "score", "strand", "frame", "id")) %>%
  filter(locus_type == "transcript",
         id %in% JL201_TE_assoc) %>%
  mutate(center = ((end - start)/2) + start,
         chrom_old = str_remove(chrom_old, "_polypolish"))

contigs_JL201 <- read_delim("../data/annotation/JL201_funannotate_annotate_18march25/annotate_corrected/telomere_pos_JL201.txt",
                            col_select = 1:2)

contigs_convert_JL201 <- read_delim("../data/annotation/JL201_funannotate_annotate_18march25/annotate_corrected/contig_conversion.txt",
                                    col_names = c("chrom_old", "chrom_new"))
telomeres_JL201 <- read_delim("../data/annotation/JL201_funannotate_annotate_18march25/annotate_corrected/telomere_pos_JL201.txt") %>%
  pivot_longer(cols = `5' telomere`:`3' telomere`,
               names_to = "side",
               values_to = "present") %>%
  filter(present == "Yes") %>%
  mutate(position = case_when(side == "3' telomere" ~ Length,
                              .default = 0))

contigs_JL201 %>%
  arrange("length") %>%
  ggplot() +
  geom_linerange(aes(y = Contig, xmin = 0, xmax = Length)) +
  geom_point(data = JL201_gtf %>%
               left_join(contigs_convert_JL201, by = "chrom_old"),
             aes(x = center, y = chrom_new),
             alpha = 0.2) +
  geom_point(data = telomeres_JL201,
             aes(x = position, y = Contig), color = "blue") +
  scale_y_discrete(limits = contigs_JL201$Contig[order(contigs_JL201$Length)]) +
  labs(x = "Position (bp)",
       y = "Contig",
       title = "A. erigeronia JL201") +
  theme_minimal() -> g_JL201
g_JL201
ggsave("../figs/retrotransposons_positions_JL201.svg", g_JL201, width = 6, height = 4)

#JL221

JL221_annot %>%
  filter(description %in% retro_annots) %>%
  pull(var = "transcript") %>%
  unique() -> JL221_TE_assoc

JL221_gtf <- read_delim("../data/annotation/JL221/JL221_assembly_v2_BRAKER3_prot_products.gtf",
                        col_names = c("chrom_old", "source", "locus_type", "start", "end", "score", "strand", "frame", "id")) %>%
  filter(locus_type == "transcript",
         id %in% JL221_TE_assoc) %>%
  mutate(center = ((end - start)/2) + start,
         chrom_old = str_remove(chrom_old, "_polypolish"))

contigs_JL221 <- read_delim("../data/annotation/JL221_funannotate_annotate_18march25/annotate_corrected/telomere_pos_JL221.txt",
                            col_select = 1:2)

contigs_convert_JL221 <- read_delim("../data/annotation/JL221_funannotate_annotate_18march25/annotate_corrected/contig_conversion.txt",
                                    col_names = c("chrom_old", "chrom_new"))
telomeres_JL221 <- read_delim("../data/annotation/JL221_funannotate_annotate_18march25/annotate_corrected/telomere_pos_JL221.txt") %>%
  pivot_longer(cols = `5' telomere`:`3' telomere`,
               names_to = "side",
               values_to = "present") %>%
  filter(present == "Yes") %>%
  mutate(position = case_when(side == "3' telomere" ~ Length,
                              .default = 0))

contigs_JL221 %>%
  arrange("length") %>%
  ggplot() +
  geom_linerange(aes(y = Contig, xmin = 0, xmax = Length)) +
  geom_point(data = JL221_gtf %>%
               left_join(contigs_convert_JL221, by = "chrom_old"),
             aes(x = center, y = chrom_new),
             alpha = 0.2) +
  geom_point(data = telomeres_JL221,
             aes(x = position, y = Contig), color = "blue") +
  scale_y_discrete(limits = contigs_JL221$Contig[order(contigs_JL221$Length)]) +
  labs(x = "Position (bp)",
       y = "Contig",
       title = "A. cardamina JL221") +
  theme_minimal() -> g_JL221
g_JL221
ggsave("../figs/retrotransposons_positions_JL221.svg", g_JL221, width = 6, height = 4)

# NB124
NB124_annot %>%
  filter(description %in% retro_annots) %>%
  pull(var = "transcript") %>%
  unique() -> NB124_TE_assoc
NB124_gtf <- read_delim("../data/annotation/NB124-2/NB124-2_assembly_v1_BRAKER3_prot_products.gtf",
                        col_names = c("chrom_old", "source", "locus_type", "start", "end", "score", "strand", "frame", "id")) %>%
  filter(locus_type == "transcript",
         id %in% NB124_TE_assoc) %>%
  mutate(center = ((end - start)/2) + start,
         chrom_old = str_remove(chrom_old, "_polypolish"))

contigs_NB124 <- read_delim("../data/annotation/NB124_funannotate_annotate_18march25/annotate_corrected/Aimania_NB124.fasta.fai",
                              col_names = c("chrom_new", "length"))

contigs_convert_NB124 <- read_delim("../data/annotation/NB124_funannotate_annotate_18march25/annotate_corrected/contig_conversion.txt",
                              col_names = c("chrom_old", "chrom_new"))
telomeres_NB124 <- read_delim("../data/annotation/NB124_funannotate_annotate_18march25/annotate_corrected/telomere_pos_NB124.txt") %>%
  pivot_longer(cols = `5' telomere`:`3' telomere`,
               names_to = "side",
               values_to = "present") %>%
  filter(present == "Yes") %>%
  mutate(position = case_when(side == "3' telomere" ~ Length,
                   .default = 0))
telomeres_NB124
contigs_NB124 %>%
  arrange("length") %>%
  ggplot() +
  geom_linerange(aes(y = chrom_new, xmin = 0, xmax = length)) +
  geom_point(data = NB124_gtf %>%
               left_join(contigs_convert_NB124, by = "chrom_old"),
             aes(x = center, y = chrom_new),
             alpha = 0.2) +
  geom_point(data = telomeres_NB124,
             aes(x = position, y = Contig), color = "blue") +
  scale_y_discrete(limits = contigs_NB124$chrom_new[order(contigs_NB124$length)]) +
  labs(x = "Position (bp)",
       y = "Contig",
       title = "A. sorghi NB124-2") +
  theme_minimal() -> g_NB124
g_NB124
ggsave("../figs/retrotransposons_positions_NB124.svg", g_NB124, width = 6, height = 4)

