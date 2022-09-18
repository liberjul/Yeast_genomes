library(ggplot2)
library(tidyr)
library(dplyr)

setwd("C:/Users/julia/Dropbox (Duke Bio_Ea)/He Lab/Julian_Liber/Yeast_genomes/")
JL201_gaps <- read.csv("./data/annotation/JL201/JL201_gap.csv") %>%
  mutate(strain = "JL201")
JL221_gaps <- read.csv("./data/annotation/JL221/JL221_gap.csv") %>%
  mutate(strain = "JL221")

comb_dat <- rbind(JL201_gaps, JL221_gaps)

comb_dat %>%
  ggplot(aes(x = gap_size, fill = orientation)) +
  facet_wrap(~strain) +
  geom_histogram(alpha = 0.6) +
  scale_x_log10() +
  labs(x = "Intergenic gap size (bp)", y = "Count") -> g
ggsave("data/annotation/gap_size_histograms.svg", g, width = 7, height = 4)

comb_dat %>%
  ggplot(aes(x = orientation, y = gap_size)) +
  facet_wrap(~strain) +
  geom_boxplot() +
  scale_y_log10() +
  labs(y = "Intergenic gap size (bp)", x = "Orientation") -> g
ggsave("data/annotation/gap_size_boxplot.svg", g, width = 7, height = 4)

JL201_gaps %>%
  filter(orientation == "convergent") %>%
  dplyr::arrange(desc(gap_size)) %>%
  .[1:10,] %>%
  write.csv(file="data/annotation/JL201/JL201_largest_conv_gaps.csv")


JL221_gaps %>%
  filter(orientation == "convergent") %>%
  dplyr::arrange(desc(gap_size)) %>%
  .[1:10,] %>%
  write.csv(file="data/annotation/JL221/JL221_largest_conv_gaps.csv")
