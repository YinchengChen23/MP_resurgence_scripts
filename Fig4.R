library(ape)
library(ggh4x)
library(metR)
library(aplot)
library(dplyr)
library(ggsci)
library(vegan)
library(ggtree)
library(treeio)
library(ggplot2)
library(tidytree)
library(tidyverse)
library(patchwork)
library(ggnewscale)

setwd('~/GCH_Hsiesh/MP3')
meta <- read.csv('data/meta.csv',sep=',')
rownames(meta) <- meta$Sample_ID
meta$label <- meta$Sample_ID
meta$Country[meta$Country %in% c("Egypt","Kenya","Tunisia")] <- 'Africa'
meta$Country[meta$Country %in% c("Denmark","Germany","Spain","France","UK")] <- 'Europe'
meta$Country[meta$Country %in% c("USA","Guatemala")] <- 'America'
meta$Country <- factor(meta$Country, levels = c('China','Taiwan','Japan','Korea','America','Europe','Africa'))
meta <- meta[meta$ST == '3',]

tree <- read.tree("data/P1_type1/ref_M129.node_labelled.final_tree.tre")
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(meta)))
root_dist <- dist.nodes(tree)
tip_idx <- 1:length(tree$tip.label)
root_dist <- root_dist[tip_idx, Ntip(tree) + 1]
names(root_dist) <- tree$tip.label
meta$root_tip <- root_dist[rownames(meta)]
meta$time_cood <- meta$Year + ifelse(is.na(meta$Month/12), 0,meta$Month/12)

print(summary(lm(root_tip~time_cood, data=meta)))

print(cor.test(meta$root_tip, meta$time_cood, method='pearson'))

bp3 <- ggplot(data = meta, aes(time_cood, root_tip, color=Country)) + geom_point(alpha=0.5) +
  scale_color_brewer(palette = 'Set1') + xlab('') + ylab('Root-to-tip Distance') + theme_bw()+ 
  geom_smooth(method = "lm", se = FALSE, color = "grey20", linewidth = 0.8, alpha=0.5, linetype = "longdash") +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_line(color = "grey80")) +
  annotate("text", label = 'Slope (rate) = 0.918', x = 1980, y =180, size = 3) +
  annotate("text", label = 'R squared = 0.327', x = 1979, y =169.5, size = 3) +
  annotate("text", label = 'Correlation Coefficient = 0.573', x = 1989.1, y =160, size = 3)



meta <- read.csv('data/meta.csv',sep=',')
submeta <- meta[meta$Country == 'Taiwan', ]
submeta$Time[submeta$Year == 2018] <- '2018'
submeta$count <- 1
summary_meta <- submeta %>%
  group_by(Subclude, Time) %>%
  summarize(total_count = sum(count, na.rm = TRUE), .groups = 'drop') 
summary_meta$Subclude <- factor(summary_meta$Subclude, levels = c('subclade_pre','subclade_post_a','subclade_post_b','subclade_post_c','subclade_post_d','other ST3','other STs'))
summary_meta$Time <- factor(summary_meta$Time, levels = c('2018','2019','2020','2023','2024 Jan','2024 Apr','2024 Jul','2024 Oct'))

bp4 <- ggplot() + 
  geom_bar(data=summary_meta,stat="identity",position="stack",aes(x=Time , y=total_count, fill=Subclude)) +
  scale_fill_manual(values = c('subclade_post_a'='#8DD3C7','subclade_post_b'='#FFFFB3','subclade_post_c'='#BEBADA',
                               'subclade_post_d'='#FB8072','subclade_pre'='#80B1D3','other ST3'='#AAAAAA','other STs'='#666666'), name='Mono') +
  theme_bw() + xlab('') + ylab('Count') +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

ppn <- bp3|bp4
ggsave(filename = "fig4.svg", plot = ppn, width = 9, height = 3, units = "in", dpi = 300)
