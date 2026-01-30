library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(stringr)
library(purrr)
library(readr)
# library(lme4)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

gene_to_product <- read_delim("../data/annotation/ncbi_cleaned_gene_products.txt", skip = 3, col_names = c("prod_name", "prod_name_clean"))
JL201_annot <- read_delim("../data/annotation/JL201/final/JL201_assembly_v2_InterProScan.tsv",
                          col_names = c("protID", "uniqueID", "length",
                                        "source", "annotation_term", "annotation_description",
                                        "start", "stop", "evalue"))
JL201_annot %>%
  filter(source %in% c("TIGRFAM", "PANTHER", "Pfam", "SUPERFAMILY")) %>%
  arrange(evalue) %>%
  pivot_wider(id_cols = protID,
              names_from = source, 
              values_from = annotation_description) -> source_info_wide_JL201

make_product_names <- function(source_info_wide){
  prod_name <- vector(length = dim(source_info_wide)[1])
  for (i in 1:dim(source_info_wide)[1]){
    if (source_info_wide$TIGRFAM[i] != "NULL"){
      prod_name[i] = str_split_i(source_info_wide$TIGRFAM[[i]][1], ": ", 2)
    } else if (! source_info_wide$PANTHER[i] %in% c("NULL", "-", "UNCHARACTERIZED")){
      prod_name[i] = tolower(source_info_wide$PANTHER[[i]][1])%>%
        str_replace_all("dna ", "DNA ") %>%
        str_replace_all("dna-", "DNA-") %>%
        str_replace_all("rna", "RNA") %>%
        str_replace_all("nad+", "NAD+") %>%
        str_replace_all("nad ", "NAD ") %>%
        str_replace_all("nadh", "NADH") %>%
        str_replace_all("atp", "ATP") %>%
        str_replace_all("adp", "ADP") %>%
        str_replace_all("amp", "AMP") %>%
        str_replace_all("gtp", "GTP") %>%
        str_replace_all("gdp", "GDP") %>%
        str_replace_all("gmp", "GMP") %>%
        str_replace_all("ctp", "CTP") %>%
        str_replace_all("cdp", "CDP") %>%
        str_replace_all("cdp", "CMP") %>%
        str_replace_all("utp", "UTP") %>%
        str_replace_all("udp", "UDP") %>%
        str_replace_all("udp", "UMP") %>%
        str_replace_all("abc transporter", "ABC transporter") %>%
        str_replace_all("ca2", "Ca2") %>%
        str_replace_all("na\\(", "Na(") %>%
        str_replace_all("h\\(", "H(") %>%
        str_replace_all("fe\\(", "Fe(") %>%
        str_replace_all("mn\\(", "Mn(") %>%
        str_replace_all("mg2", "Mg2") %>%
        str_replace_all("zn-", "Zn-")
    }
    else if (source_info_wide$SUPERFAMILY[i] != "NULL"){
      prod_name[i] = source_info_wide$SUPERFAMILY[[i]][1]
      if (str_detect(prod_name[i], "ases$")){
        prod_name[i] <- paste0(prod_name[i], " superfamily")
        }
    }
    else if (source_info_wide$Pfam[i] != "NULL"){
      prod_name[i] = source_info_wide$Pfam[[i]][1]
    } else {
      prod_name[i] = "hypothetical protein"
    }
  }
  for (i in 1:dim(source_info_wide)[1]){
    if (str_detect(prod_name[i], "domain$")){
      prod_name[i] <- paste0(prod_name[i], "-containing protein")
    }
    if (str_detect(prod_name[i], "unnamed product")){
      prod_name[i] <- "hypothetical protein"
    }
    if (str_detect(prod_name[i], "wd40")){
      prod_name[i] <- str_replace(prod_name[i], "wd40", "WD40")
    }
    if (str_detect(prod_name[i], "  ")){
      prod_name[i] <- str_replace(prod_name[i], "  ", " ")
    }
    if (str_detect(prod_name[i], "Domain of [u|U]nknown function")){
      prod_name[i] <- paste0(prod_name[i], "-containing protein")
    }
    if (str_detect(prod_name[i], ",")){
      prod_name[i] <- str_replace(prod_name[i], ",", "%2C")
    }
  }
  return(prod_name)
}

source_info_wide_JL201$prod_name <- make_product_names(source_info_wide_JL201)

gtf_JL201 <- read_delim("../data/annotation/JL201/final/JL201_assembly_v2_BRAKER3_RNA.gtf", 
                        col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

gtf_JL201%>%
  mutate(protID = str_extract(attribute, "g[0-9]*\\.t[1-9]")) %>%
  left_join(source_info_wide_JL201 %>% dplyr::select(protID, prod_name),
            by = "protID") %>%
  mutate(prod_name = case_when(!is.na(prod_name) ~ prod_name,
                               .default = "hypothetical protein"),
         attribute = case_when(feature == "CDS" & !is.na(prod_name) ~ paste0(attribute, " product \"", prod_name, "\";"),
                             .default =  attribute))-> gtf_JL201_w_prod
  
gtf_JL201_w_prod %>%
  dplyr::select(seqname:attribute) %>%
  write_delim("../data/annotation/JL201/final/JL201_assembly_v2_BRAKER3_RNA_products.gtf",
              na = "", col_names = F, delim = "\t")
###

JL221_annot <- read_delim("../data/annotation/JL221/JL221_assembly_v2_BRAKER3_prot_InterProScan.tsv",
                          col_names = c("protID", "uniqueID", "length",
                                        "source", "annotation_term", "annotation_description",
                                        "start", "stop", "evalue"))
JL221_annot %>%
  filter(source %in% c("TIGRFAM", "PANTHER", "Pfam", "SUPERFAMILY")) %>%
  arrange(evalue) %>%
  pivot_wider(id_cols = protID,
              names_from = source, 
              values_from = annotation_description) -> source_info_wide_JL221
source_info_wide_JL221$prod_name <- make_product_names(source_info_wide_JL221)

gtf_JL221 <- read_delim("../data/annotation/JL221/JL221_assembly_v2_BRAKER3_prot.gtf", 
                        col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

gtf_JL221%>%
  mutate(protID = str_extract(attribute, "g[0-9]*\\.t[1-9]")) %>%
  left_join(source_info_wide_JL221 %>% dplyr::select(protID, prod_name),
            by = "protID") %>%
  mutate(prod_name = case_when(!is.na(prod_name) ~ prod_name,
                               .default = "hypothetical protein"),
         attribute = case_when(feature == "CDS" & !is.na(prod_name) ~ paste0(attribute, " product \"", prod_name, "\";"),
                               .default =  attribute))-> gtf_JL221_w_prod

gtf_JL221_w_prod %>%
  dplyr::select(seqname:attribute) %>%
  write_delim("../data/annotation/JL221/JL221_assembly_v2_BRAKER3_prot_products.gtf",
              na = "", col_names = F, delim = "\t")

###

NB1242_annot <- read_delim("../data/annotation/NB124-2/NB124-2_assembly_v1_BRAKER3_InterProScan.tsv",
                          col_names = c("protID", "uniqueID", "length",
                                        "source", "annotation_term", "annotation_description",
                                        "start", "stop", "evalue"))
NB1242_annot %>%
  filter(source %in% c("TIGRFAM", "PANTHER", "Pfam", "SUPERFAMILY")) %>%
  arrange(evalue) %>%
  pivot_wider(id_cols = protID,
              names_from = source, 
              values_from = annotation_description) -> source_info_wide_NB1242
source_info_wide_NB1242$prod_name <- make_product_names(source_info_wide_NB1242)

gtf_NB1242 <- read_delim("../data/annotation/NB124-2/NB124-2_assembly_v1_BRAKER3_prot.gtf", 
                        col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

gtf_NB1242%>%
  mutate(protID = str_extract(attribute, "g[0-9]*\\.t[1-9]")) %>%
  left_join(source_info_wide_NB1242 %>% dplyr::select(protID, prod_name),
            by = "protID") %>%
  mutate(prod_name = case_when(!is.na(prod_name) ~ prod_name,
                               .default = "hypothetical protein"),
         attribute = case_when(feature == "CDS" & !is.na(prod_name) ~ paste0(attribute, " product \"", prod_name, "\";"),
                               .default =  attribute))-> gtf_NB1242_w_prod

gtf_NB1242_w_prod %>%
  dplyr::select(seqname:attribute) %>%
  write_delim("../data/annotation/NB124-2/NB124-2_assembly_v1_BRAKER3_prot_products.gtf",
              na = "", col_names = F, delim = "\t")

NB1242_annot %>%
  filter(!protID %in% source_info_wide_NB1242$protID) %>%
  pull(var = "protID") %>%
  unique() -> other_protID_NB1242

NB1242_annot %>%
  filter(protID %in% other_protID_NB1242) %>%
  filter(source == "ProSiteProfiles") %>%
  view()
  group_by(source) %>%
  summarize(count = n())
