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


GOdata_up <- new("topGOdata", ontology = "MF", allGenes = geneListAp_up,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher_up <- runTest(GOdata_up, algorithm = "classic", statistic = "fisher")

geneListAp <- factor(as.integer(unique(gene_to_gene$NBB_id) %in% sig_gene_list_NBB))
names(geneListAp) <- unique(gene_to_gene$NBB_id)

###
IPS_cols <- c("transcript", "MD5", "length", "source", "accession", "description", "start", "stop", "score", "status", "date", "IP_accession", "IP_description")

JL201_annot <- read_delim("../data/annotation/JL201/final/JL201_assembly_v2_InterProScan.tsv",
                          col_names = IPS_cols)
JL221_annot <- read_delim("../data/annotation/JL221/JL221_assembly_v2_BRAKER3_prot_InterProScan.tsv",
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
  mutate(dir = case_when(log2_ratio_aim_to_mbm < 0 ~ "Depleted in Aimania",
                         .default = "Enriched in Aimania"),
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
       fill = "Log2 enrichment in\nAimania vs other\nMicrobotryomycetes") -> g
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
                         group = c(i, "other Aimania"))

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
         strain_name = factor(case_when(strain == "JL201" ~ "Aimania erigeronia JL201",
                                 strain == "JL221" ~ "Aimania cardamina JL221",
                                 strain == "NB124" ~ "Aimania sorghi NB124-2"),
                              levels = c("Aimania erigeronia JL201",
                                         "Aimania cardamina JL221",
                                         "Aimania sorghi NB124-2")),
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
       fill = "Log2 ratio in\nisolate vs other\nAimania spp.",
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
