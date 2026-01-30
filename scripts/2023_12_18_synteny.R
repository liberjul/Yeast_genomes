library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

gap_n <- 800000
contigs_JL201 <- read.csv("../data/annotation/JL201/final/JL201_assembly_v2_contigs.csv") %>%
  dplyr::rename(chrom.pos = length..bp.) %>%
  arrange(chrom.pos) %>%
  mutate(contig = str_extract(contig, "contig_[0-9]*"),
         chrom.pos = as.numeric(str_remove_all(chrom.pos, ","))) %>%
  arrange(-chrom.pos) %>%
  mutate(#query.acc = contig,
    subj.acc = contig,
    x.pos = cumsum(chrom.pos) + row_number()*gap_n)
contigs_JL221_to <- read.csv("../data/2023_12_16_Withers_Q5F_results_JL221_ONT/Withers_Q5F_results/Withers_Q5F_1/annotation/2023_12_18_JL221_ONT_only_contigs.csv") %>%
  dplyr::rename(chrom.pos = length..bp.) %>%
  arrange(chrom.pos) %>%
  mutate(chrom.pos = as.numeric(str_remove_all(chrom.pos, ","))) %>%
  arrange(-chrom.pos) %>%
  mutate(#query.acc = contig,
         subj.acc = contig,
         x.pos = cumsum(chrom.pos) + row_number()*gap_n)
contigs_JL221_from <- read.csv("../data/2023_12_16_Withers_Q5F_results_JL221_ONT/Withers_Q5F_results/Withers_Q5F_1/annotation/2023_12_18_JL221_ONT_only_contigs.csv") %>%
  dplyr::rename(chrom.pos = length..bp.) %>%
  arrange(chrom.pos) %>%
  mutate(chrom.pos = as.numeric(str_remove_all(chrom.pos, ","))) %>%
  arrange(-chrom.pos) %>%
  mutate(query.acc = contig,
    #subj.acc = contig,
    x.pos = cumsum(chrom.pos) + row_number()*gap_n)

contigs_NB124_2 <- read.csv("../data/raw_data/2023_12_16_Withers_B5M_results_NB124-2_ONT/Withers_B5M_results/Withers_B5M_1/annotation/2023_12_18_NB124-2_ONT_only_contigs.csv") %>%
  dplyr::rename(chrom.pos = length..bp.) %>%
  mutate(chrom.pos = as.numeric(str_remove_all(chrom.pos, ","))) %>%
  arrange(-chrom.pos) %>%
  mutate(query.acc = contig,
         # subj.acc = contig,
         x.pos = cumsum(chrom.pos) + row_number()*gap_n)


#Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, 
#q. start, q. end, s. start, s. end, evalue, bit score
### NB124-2 to JL201
blast_NB124_2_to_JL201 <- read.delim("../data/annotation/JL201/final/NB124-2_to_JL201_all_v_all_blastn.txt",
                                     header = F,
                                     col.names = c("query.acc", "subj.acc", "perc.id", "aln.len", "mismatches",
                                                   "gap.opens", "q.start", "q.end", "s.start", "s.end", "e.value", "bit.score")) %>%
  left_join(contigs_NB124_2, by = "query.acc") %>% #these have the suffix x
  left_join(contigs_JL201, by = "subj.acc") %>% #these have the suffix y
  mutate(x1.start = q.start + x.pos.x,
         x1.end = q.end + x.pos.x, 
         x2.start = s.start + x.pos.y,
         x2.end = s.end + x.pos.y,
         x1.mid =  (q.start + q.end)/2 + x.pos.x,
         x2.mid =  (s.start + s.end)/2 + x.pos.y) %>%
  rowwise() %>%
  # mutate(hit_id = as.integer(runif(n = 1, min = 1, max = 2^16))) %>%
  pivot_longer(cols = x1.start:x2.end, names_to = "x.coord", values_to = "x.pos") %>%
  mutate(y.pos = case_when(x.coord %in% c("x1.start", "x1.end", "x1.mid") ~ 0,
                           x.coord %in% c("x2.start", "x2.end", "x1.mid") ~ 1))

blast_NB124_2_to_JL201 %>%
  ggplot() +
  # geom_polygon(aes(x = x.pos, y = y.pos, fill = perc.id)) +
  geom_segment(aes(x = x1.mid, xend = x2.mid, color = query.acc), y = 1, yend = 0) +
  geom_text(data = contigs_JL201,
            aes(x = x.pos + 0.5*chrom.pos, y = 0, label = contig), 
            hjust = 1, angle = 90, vjust = 0.5) +
  geom_text(data = contigs_NB124_2,
            aes(x = x.pos + 0.5*chrom.pos, y = 1, label = contig), 
            hjust = 0, angle = 90, vjust = 0.5) +
  geom_segment(data = contigs_JL201,
            aes(x = x.pos,
                xend = x.pos + chrom.pos,
                y = 0,
                yend = 0)) +
  geom_segment(data = contigs_NB124_2,
            aes(x = x.pos,
                xend = x.pos + chrom.pos,
                y = 1,
                yend = 1)) +
  annotate(geom = "text", x = 2.5e6, y = 1, label = "NB124-2", hjust = 1, vjust = 0.5) +
  annotate(geom = "text", x = 2.5e6, y = 0, label = "JL201", hjust = 1, vjust = 0.5) +
  ylim(-0.2, 1.2) + xlim(-1e6, 4e7) +
  theme_void() +
  theme(legend.position = "none") -> g
g
ggsave("2023_12_19_NB124-2_to_JL201_blastn_synteny.svg", g, width = 8, height = 4)

blast_NB124_2_to_JL201 %>%
  ggplot() +
  # geom_polygon(aes(x = x.pos, y = y.pos, fill = perc.id)) +
  geom_segment(aes(x = x1.mid, xend = x2.mid, color = perc.id), y = 1, yend = 0) +
  geom_text(data = contigs_JL201,
            aes(x = x.pos + 0.5*chrom.pos, y = 0, label = contig), 
            hjust = 1, angle = 90, vjust = 0.5) +
  geom_text(data = contigs_NB124_2,
            aes(x = x.pos + 0.5*chrom.pos, y = 1, label = contig), 
            hjust = 0, angle = 90, vjust = 0.5) +
  geom_segment(data = contigs_JL201,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 0,
                   yend = 0)) +
  geom_segment(data = contigs_NB124_2,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 1,
                   yend = 1)) +
  annotate(geom = "text", x = 2.5e6, y = 1, label = "NB124-2", hjust = 1, vjust = 0.5) +
  annotate(geom = "text", x = 2.5e6, y = 0, label = "JL221", hjust = 1, vjust = 0.5) +
  ylim(-0.2, 1.2) + xlim(-1e6, 4e7) +
  theme_void()-> g
g
ggsave("2023_12_19_NB124-2_to_JL201_blastn_synteny_percid.svg", g, width = 8, height = 4)

###NB124-2 to JL221

blast_NB124_2_to_JL221 <- read.delim("../data/2023_12_16_Withers_Q5F_results_JL221_ONT/Withers_Q5F_results/Withers_Q5F_1/annotation/NB124-2_to_JL221_all_v_all_blastn.txt",
                                     header = F,
                                     col.names = c("query.acc", "subj.acc", "perc.id", "aln.len", "mismatches",
                                                   "gap.opens", "q.start", "q.end", "s.start", "s.end", "e.value", "bit.score")) %>%
  left_join(contigs_NB124_2, by = "query.acc") %>% #these have the suffix x
  left_join(contigs_JL221_to, by = "subj.acc") %>% #these have the suffix y
  mutate(x1.start = q.start + x.pos.x,
         x1.end = q.end + x.pos.x, 
         x2.start = s.start + x.pos.y,
         x2.end = s.end + x.pos.y,
         x1.mid =  (q.start + q.end)/2 + x.pos.x,
         x2.mid =  (s.start + s.end)/2 + x.pos.y) %>%
  rowwise() %>%
  # mutate(hit_id = as.integer(runif(n = 1, min = 1, max = 2^16))) %>%
  pivot_longer(cols = x1.start:x2.end, names_to = "x.coord", values_to = "x.pos") %>%
  mutate(y.pos = case_when(x.coord %in% c("x1.start", "x1.end", "x1.mid") ~ 0,
                           x.coord %in% c("x2.start", "x2.end", "x1.mid") ~ 1))

blast_NB124_2_to_JL221 %>%
  ggplot() +
  # geom_polygon(aes(x = x.pos, y = y.pos, fill = perc.id)) +
  geom_segment(aes(x = x1.mid, xend = x2.mid, color = query.acc), y = 1, yend = 0) +
  geom_text(data = contigs_JL221_to,
            aes(x = x.pos + 0.5*chrom.pos, y = 0, label = contig), 
            hjust = 1, angle = 90, vjust = 0.5) +
  geom_text(data = contigs_NB124_2,
            aes(x = x.pos + 0.5*chrom.pos, y = 1, label = contig), 
            hjust = 0, angle = 90, vjust = 0.5) +
  geom_segment(data = contigs_JL221_to,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 0,
                   yend = 0)) +
  geom_segment(data = contigs_NB124_2,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 1,
                   yend = 1)) +
  annotate(geom = "text", x = 2.5e6, y = 1, label = "NB124-2", hjust = 1, vjust = 0.5) +
  annotate(geom = "text", x = 2.5e6, y = 0, label = "JL221", hjust = 1, vjust = 0.5) +
  ylim(-0.2, 1.2) + xlim(-1e6, 4e7) +
  theme_void() +
  theme(legend.position = "none")-> g
g
ggsave("2023_12_19_NB124-2_to_JL221_blastn_synteny.svg", g, width = 8, height = 4)

blast_NB124_2_to_JL221 %>%
  ggplot() +
  # geom_polygon(aes(x = x.pos, y = y.pos, fill = perc.id)) +
  geom_segment(aes(x = x1.mid, xend = x2.mid, color = perc.id), y = 1, yend = 0) +
  geom_text(data = contigs_JL221_to,
            aes(x = x.pos + 0.5*chrom.pos, y = 0, label = contig), 
            hjust = 1, angle = 90, vjust = 0.5) +
  geom_text(data = contigs_NB124_2,
            aes(x = x.pos + 0.5*chrom.pos, y = 1, label = contig), 
            hjust = 0, angle = 90, vjust = 0.5) +
  geom_segment(data = contigs_JL221_to,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 0,
                   yend = 0)) +
  geom_segment(data = contigs_NB124_2,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 1,
                   yend = 1)) +
  annotate(geom = "text", x = 2.5e6, y = 1, label = "NB124-2", hjust = 1, vjust = 0.5) +
  annotate(geom = "text", x = 2.5e6, y = 0, label = "JL221", hjust = 1, vjust = 0.5) +
  ylim(-0.2, 1.2) + xlim(-1e6, 4e7) +
  theme_void()-> g
g
ggsave("2023_12_19_NB124-2_to_JL221_blastn_synteny_percid.svg", g, width = 8, height = 4)

### JL221 to JL201
blast_JL221_to_JL201 <- read.delim("../data/annotation/JL201/final/JL221_to_JL201_all_v_all_blastn.txt",
                                     header = F,
                                     col.names = c("query.acc", "subj.acc", "perc.id", "aln.len", "mismatches",
                                                   "gap.opens", "q.start", "q.end", "s.start", "s.end", "e.value", "bit.score")) %>%
  left_join(contigs_JL221_from, by = "query.acc") %>% #these have the suffix x
  left_join(contigs_JL201, by = "subj.acc") %>% #these have the suffix y
  mutate(x1.start = q.start + x.pos.x,
         x1.end = q.end + x.pos.x, 
         x2.start = s.start + x.pos.y,
         x2.end = s.end + x.pos.y,
         x1.mid =  (q.start + q.end)/2 + x.pos.x,
         x2.mid =  (s.start + s.end)/2 + x.pos.y) %>%
  rowwise() %>%
  # mutate(hit_id = as.integer(runif(n = 1, min = 1, max = 2^16))) %>%
  pivot_longer(cols = x1.start:x2.end, names_to = "x.coord", values_to = "x.pos") %>%
  mutate(y.pos = case_when(x.coord %in% c("x1.start", "x1.end", "x1.mid") ~ 0,
                           x.coord %in% c("x2.start", "x2.end", "x1.mid") ~ 1))

blast_JL221_to_JL201 %>%
  ggplot() +
  # geom_polygon(aes(x = x.pos, y = y.pos, fill = perc.id)) +
  geom_segment(aes(x = x1.mid, xend = x2.mid, color = query.acc), y = 1, yend = 0) +
  geom_text(data = contigs_JL201,
            aes(x = x.pos + 0.5*chrom.pos, y = 0, label = contig), 
            hjust = 1, angle = 90, vjust = 0.5) +
  geom_text(data = contigs_JL221_from,
            aes(x = x.pos + 0.5*chrom.pos, y = 1, label = contig), 
            hjust = 0, angle = 90, vjust = 0.5) +
  geom_segment(data = contigs_JL201,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 0,
                   yend = 0)) +
  geom_segment(data = contigs_JL221_from,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 1,
                   yend = 1)) +
  annotate(geom = "text", x = 2.5e6, y = 1, label = "JL221", hjust = 1, vjust = 0.5) +
  annotate(geom = "text", x = 2.5e6, y = 0, label = "JL201", hjust = 1, vjust = 0.5) +
  ylim(-0.2, 1.2) + xlim(-1e6, 4e7) +
  theme_void() +
  theme(legend.position = "none")-> g
g
ggsave("2023_12_19_JL221_JL201_blastn_synteny.svg", g, width = 8, height = 4)

blast_JL221_to_JL201 %>%
  ggplot() +
  # geom_polygon(aes(x = x.pos, y = y.pos, fill = perc.id)) +
  geom_segment(aes(x = x1.mid, xend = x2.mid, color = perc.id), y = 1, yend = 0) +
  geom_text(data = contigs_JL201,
            aes(x = x.pos + 0.5*chrom.pos, y = 0, label = contig), 
            hjust = 1, angle = 90, vjust = 0.5) +
  geom_text(data = contigs_JL221_from,
            aes(x = x.pos + 0.5*chrom.pos, y = 1, label = contig), 
            hjust = 0, angle = 90, vjust = 0.5) +
  geom_segment(data = contigs_JL201,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 0,
                   yend = 0)) +
  geom_segment(data = contigs_JL221_from,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 1,
                   yend = 1)) +
  annotate(geom = "text", x = 2.5e6, y = 1, label = "JL221", hjust = 1, vjust = 0.5) +
  annotate(geom = "text", x = 2.5e6, y = 0, label = "JL201", hjust = 1, vjust = 0.5) +
  ylim(-0.2, 1.2) + xlim(-1e6, 4e7) +
  theme_void()-> g
g
ggsave("2023_12_19_JL201_JL221_blastn_synteny_percid.svg", g, width = 8, height = 4)
# Going to try to make circos plots
# As of 1/26/26, not functional

bt_gen_gap <- 5e6

diff_124 <- max(contigs_NB124_2$chrom.pos + contigs_NB124_2$x.pos + bt_gen_gap)
###
blast_NB124_2_to_JL201 <- read.delim("../data/annotation/JL201/final/NB124-2_to_JL201_all_v_all_blastn.txt",
                                     header = F,
                                     col.names = c("query.acc", "subj.acc", "perc.id", "aln.len", "mismatches",
                                                   "gap.opens", "q.start", "q.end", "s.start", "s.end", "e.value", "bit.score")) %>%
  left_join(contigs_NB124_2, by = "query.acc") %>% #these have the suffix x
  left_join(contigs_JL201, by = "subj.acc") %>% #these have the suffix y
  mutate(x1.start = q.start + x.pos.x,
         x1.end = q.end + x.pos.x, 
         x2.start = s.start + x.pos.y + diff_124,
         x2.end = s.end + x.pos.y + diff_124,
         x1.mid =  (q.start + q.end)/2 + x.pos.x,
         x2.mid =  (s.start + s.end)/2 + x.pos.y + diff_124) %>%
  rowwise() %>%
  # mutate(hit_id = as.integer(runif(n = 1, min = 1, max = 2^16))) %>%
  pivot_longer(cols = x1.start:x2.end, names_to = "x.coord", values_to = "x.pos") %>%
  mutate(y.pos = case_when(x.coord %in% c("x1.start", "x1.end", "x1.mid") ~ 0,
                           x.coord %in% c("x2.start", "x2.end", "x1.mid") ~ 1),
         
         )
blast_NB124_2_to_JL201 %>%
  ggplot() +
  # geom_polygon(aes(x = x.pos, y = y.pos, fill = perc.id)) +
  # geom_segment(aes(x = x1.mid, xend = x2.mid, color = query.acc),
  #              y = 0, yend = 0) +
  # geom_text(data = contigs_JL201,
  #           aes(x = x.pos + 0.5*chrom.pos + diff_124, y = 0, label = contig),
  #           hjust = 1, angle = 0, vjust = 0.5) +
  # geom_text(data = contigs_NB124_2,
  #           aes(x = x.pos + 0.5*chrom.pos, y = 0, label = contig),
  #           hjust = 0, angle = 0, vjust = 0.5) +
  geom_segment(data = contigs_JL201,
               aes(x = x.pos + diff_124,
                   xend = x.pos + chrom.pos + diff_124,
                   y = 1,
                   yend = 1), color ="red") +
  geom_segment(data = contigs_NB124_2,
               aes(x = x.pos,
                   xend = x.pos + chrom.pos,
                   y = 1,
                   yend = 1), color = "blue") +
  geom_segment(x = 1e7, xend = 7e7, color = "black", y = 1, yend = 0) +
  annotate(geom = "text", x = 7.5e8, y = 0, label = "NB124-2", hjust = 1, vjust = 0.5) +
  annotate(geom = "text", x = 2.5e8, y = 0, label = "JL201", hjust = 1, vjust = 0.5) +
  # ylim(-0.2, 1.2) + xlim(-1e6, 4e7) +
  xlim(0, diff_124 + max(contigs_JL201$chrom.pos + contigs_JL201$x.pos)) +
  theme_void() +
  coord_polar() -> g
g
g

hist(blast_NB124_2_to_JL201$x1.mid)

hist(blast_NB124_2_to_JL201$x2.mid)

hist(contigs_JL201$x.pos + diff_124)

diff_124 + max(contigs_JL201$chrom.pos + contigs_JL201$x.pos) # 8.6e7
