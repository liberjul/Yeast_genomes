library(phangorn)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(readr)
library(ggtree)
library(treeio)
library(ggrepel)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tree1 <- read.beast("../trees/20251024_mlst_aimania_iqtree.nexus")


g <- ggtree(tree1)
g$data %>%
  dplyr::rename(support = label.y,
                label = label.x) -> tree_data1

g$data <- tree_data1
g + geom_tiplab(aes(label = str_replace_all(label, "_", " ")), size = 3.1,
               fontface = "italic") -> g2
ggsave("../trees/2025_10_24_tree_nodes_labeled.svg",
       g2 + geom_text_repel(aes(label = node, x = x, y=y)) +xlim(-0.2, 0.7),
       height = 24, width = 10) 

clade_labs <- read.csv("../trees/2025_10_24_clade_labels.csv")
g2
tip_vec <- tree_data1 %>% filter(isTip) %>% pull(var = "node")
g -> g3
# g3$data %>% view()

clade_labs$min_x <- NA
for (node in clade_labs$Node){
  node_off <- offspring(tree1, node)
  tip_nodes <- node_off[node_off %in% tip_vec]
  tip_count <- length(tip_nodes)
  tree_data1 %>%
    filter(node %in% tip_nodes) %>%
    pull(var = "x") %>% min() -> min_pos
  clade_labs$min_x[clade_labs$Node == node] <- min_pos
  print(node)
  print(tip_count)
  scaleClade(g3, node = node, scale = 2/tip_count) -> g3
}
for (node in clade_labs$Node){
  ggtree::collapse(g3, node, 'min') -> g3
}
g3 +
  geom_tiplab(data = g3$data %>% filter(!str_detect(label, "Novel_sp")),
              aes(label = case_when(str_detect(label, "_Rhod") ~ str_replace(label, "_Rhodotorula__", "\"Rhodotorula\" "),
                                               .default = str_replace_all(label, "_", " "))),
                                    size = 3.1,
              fontface = "italic") +
  geom_tiplab(data = g3$data %>% filter(str_detect(label, "Novel_sp")),
              aes(label = case_when(label == "Novel_sp_JL201" ~ "Aimania erigeronia sp. nov. JL201",
                                    label == "Novel_sp_JL221" ~ "Aimania cardamina sp. nov. JL221",
                                    label == "Novel_sp_NB124-2" ~ "Aimania sorghi sp. nov. NB124-2")),
              size = 2.9,
              fontface = "bold.italic") +
  geom_text_repel(data = g3$data %>% filter(!isTip, !support %in% c("100", "1/100")),
                  aes(label = support), hjust = 1.5, min.segment.length = 1,
                  vjust = -0.1, segment.alpha = 0.2,
                  ) +
  geom_nodepoint(data = g3$data %>% filter(!isTip, support %in% c("100", "1/100")),
                 size = 1.5) +
  geom_cladelab(data = clade_labs,
                mapping = aes(node = Node, label = clade_label),
                offset.text = 0.068,
                # align = TRUE,
                vjust = -0.1, size = 1.5,
                fontface = "italic", fill = NA) +
  xlim(0, 0.7)-> g4

g4
ggsave("../trees/2025_10_24_collapsed_clades_support_no_AR_PF.svg", g4, height = 10, width = 7.5)

offspring(tree1, 193)
child(tree1, 193)

g4$data %>% view()
