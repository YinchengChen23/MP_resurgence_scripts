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
meta <- meta[meta$AdhesinType == 'type1',]

tree <- read.tree("data/P1_type1/ref_M129.node_labelled.final_tree.tre")
tree <- drop.tip(tree, setdiff(tree$tip.label,rownames(meta)))
y <- full_join(as_tibble(tree),as_tibble(meta), by = 'label')
p11 <- ggtree(tree, right = FALSE) %<+% y + 
  geom_tippoint(aes(color=Clade), size=1) +
  scale_color_manual(values = c("#8DD3C7","#FFFFB3","#BEBADA"))
# "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462"

d2 <- data.frame(node=1:Nnode2(tree), lty = 1); d2[c(1233,1234), 2] <- 2
p11 <- p11 %<+% d2 + aes(linetype=I(lty)) + layout_dendrogram()

tip_df <- p11$data[p11$data$isTip, ]
tip_df$x <- tip_df$x*0.975
p11 <- p11 + geom_segment(data = tip_df, aes(x = x, xend = max(x), y = y, yend = y, linetype=Clade),
                          linewidth = 0.2, color = "grey40") +
  theme(plot.margin = margin(0, 0, -5, 0), legend.position = 'top') +
  guides(linetype = "none")



get_taxa_name(p11) -> strain_name
supplementary <- meta[strain_name,]
supplementary$Sample_ID <-factor(supplementary$Sample_ID, levels = rev(strain_name), ordered = T)

supplementary$ST_ <- paste0('ST',supplementary$ST)
supplementary$ST_[!supplementary$ST_ %in% c('ST20','ST7','ST1','ST2','ST17','ST14','ST3')] <- 'others'
suppl1 <- supplementary[,c('Sample_ID','Country','ST_','Time','Usage')]
suppl1$country <- 'Country'
suppl1$st <- 'ST'
suppl1$time <- 'Time'
suppl1$usage <- 'Source'
suppl1$ST_ <- factor(suppl1$ST_, levels = c('ST3','ST14','ST17','ST2','ST1','ST7','ST20','others'))
suppl1$Source <- factor(suppl1$Usage, levels = c('This study','Public data'))
suppl1$Time <- factor(suppl1$Time, levels = c('< 2010','< 2019','2019','2020','2023','2024 Jan','2024 Apr','2024 Jul','2024 Oct'))
suppl1$era <- ifelse(suppl1$Time %in% c('< 2010','< 2019','2019','2020'),'per','post')


p12 <- ggplot(data = suppl1) + 
  geom_tile(aes(x = Sample_ID, y = usage, fill = Source), color = NA, alpha = 0.7, linewidth = 0.5) + 
  scale_fill_manual(values = c('This study'='#377EB8','Public data'='#FFFFFF'), name='Source', guide = 'none') +
  new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = time, fill = Time), color = NA, alpha = 0.7, linewidth = 0.5) + 
  scale_fill_manual(values = c('< 2010'='#AAAAFF','< 2019'='#5555FF','2019'='#0000FF','2020'='#FFDADA',
                               '2023'='#FFB6B6','2024 Jan'='#FF9191','2024 Apr'='#FF6D6D','2024 Jul'='#FF2424','2024 Oct'='#FF0000'),
                    name='Time') +
  new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = st, fill = ST_),color = NA, alpha = 0.7, linewidth = 0.5) +
  scale_fill_manual(values = c('ST3'="#66C2A5",'ST14'="#FC8D62",'ST17'="#8DA0CB",'ST2'="#E78AC3",'ST1="#A6D854"','ST7'="#FFD92F",'ST20'="#E5C494",'others'="#B3B3B3"),name='ST') +
  new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = country, fill = Country), color = NA, alpha = 0.7, linewidth = 0.5) + 
  scale_fill_manual(values = c("China"="#E41A1C", "Taiwan"="#377EB8", "Japan"="#4DAF4A", "Korea"="#984EA3",
                               "America"="#FF7F00","Europe"="#FFFF33","Africa"="#A65628"),name='Country') + 
  
  labs(x = "", y = "")  + scale_y_discrete(limits = c("Time", "Country", "Source","ST")) + #theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.background = element_rect(fill = "white"),
        legend.position = 'None',
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0))


ppn <- p11/p12 + plot_layout(height=c(6,1))
ppn
ggsave(filename = "fig3_a1.svg", plot = ppn, width = 9, height = 4, units = "in", dpi = 300)



meta <- read.csv('data/meta.csv',sep=',')
rownames(meta) <- meta$Sample_ID
meta$label <- meta$Sample_ID
meta$Country[meta$Country %in% c("Egypt","Kenya","Tunisia")] <- 'Africa'
meta$Country[meta$Country %in% c("Denmark","Germany","Spain","France","UK")] <- 'Europe'
meta$Country[meta$Country %in% c("USA","Guatemala")] <- 'America'
meta$Country <- factor(meta$Country, levels = c('China','Taiwan','Japan','Korea','America','Europe','Africa'))
meta <- meta[meta$AdhesinType == 'type2',]

tree <- read.tree("data/P1_type2/p2FH.node_labelled.final_tree.tre")
tree <- drop.tip(tree, setdiff(tree$tip.label,rownames(meta)))
y <- full_join(as_tibble(tree),as_tibble(meta), by = 'label')
p11 <- ggtree(tree, right = FALSE) %<+% y + 
  geom_tippoint(aes(color=Clade), size=1) +
  scale_color_manual(values = c("#FB8072","#80B1D3"))

d2 <- data.frame(node=1:Nnode2(tree), lty = 1)#; d2[c(1233,1234), 2] <- 2
p11 <- p11 %<+% d2 + aes(linetype=I(lty)) + layout_dendrogram()

tip_df <- p11$data[p11$data$isTip, ]
tip_df$x <- tip_df$x*0.975
p11 <- p11 + geom_segment(data = tip_df, aes(x = x, xend = max(x), y = y, yend = y, linetype=Clade),
                          linewidth = 0.2, color = "grey40") +
  theme(plot.margin = margin(0, 0, -5, 0), legend.position = 'top') +
  guides(linetype = "none")

get_taxa_name(p11) -> strain_name
supplementary <- meta[strain_name,]
supplementary$Sample_ID <-factor(supplementary$Sample_ID, levels = rev(strain_name), ordered = T)

supplementary$ST_ <- paste0('ST',supplementary$ST)
supplementary$ST_[!supplementary$ST_ %in% c('ST20','ST7','ST1','ST2','ST17','ST14','ST3')] <- 'others'
suppl1 <- supplementary[,c('Sample_ID','Country','ST_','Time','Usage')]
suppl1$country <- 'Country'
suppl1$st <- 'ST'
suppl1$time <- 'Time'
suppl1$usage <- 'Source'
suppl1$ST_ <- factor(suppl1$ST_, levels = c('ST3','ST14','ST17','ST2','ST1','ST7','ST20','others'))
suppl1$Source <- factor(suppl1$Usage, levels = c('This study','Public data'))
suppl1$Time <- factor(suppl1$Time, levels = c('< 2010','< 2019','2019','2020','2023','2024 Jan','2024 Apr','2024 Jul','2024 Oct'))
suppl1$era <- ifelse(suppl1$Time %in% c('< 2010','< 2019','2019','2020'),'per','post')

p12 <- ggplot(data = suppl1) + 
  geom_tile(aes(x = Sample_ID, y = usage, fill = Source), color = NA, alpha = 0.7, linewidth = 0.5) + 
  scale_fill_manual(values = c('This study'='#377EB8','Public data'='#FFFFFF'), name='Source', guide = 'none') +
  new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = time, fill = Time), color = NA, alpha = 0.7, linewidth = 0.5) + 
  scale_fill_manual(values = c('< 2010'='#AAAAFF','< 2019'='#5555FF','2019'='#0000FF','2020'='#FFDADA',
                               '2023'='#FFB6B6','2024 Jan'='#FF9191','2024 Apr'='#FF6D6D','2024 Jul'='#FF2424','2024 Oct'='#FF0000'),
                    name='Time') +
  new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = st, fill = ST_),color = NA, alpha = 0.7, linewidth = 0.5) +
  scale_fill_manual(values = c('ST3'="#66C2A5",'ST14'="#FC8D62",'ST17'="#8DA0CB",'ST2'="#E78AC3",'ST1="#A6D854"','ST7'="#FFD92F",'ST20'="#E5C494",'others'="#B3B3B3"),name='ST') +
  new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = country, fill = Country), color = NA, alpha = 0.7, linewidth = 0.5) + 
  scale_fill_manual(values = c("China"="#E41A1C", "Taiwan"="#377EB8", "Japan"="#4DAF4A", "Korea"="#984EA3",
                               "America"="#FF7F00","Europe"="#FFFF33","Africa"="#A65628"),name='Country') + 
  
  
  labs(x = "", y = "")  + scale_y_discrete(limits = c("Time", "Country", "Source","ST")) + #theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.background = element_rect(fill = "white"),
        legend.position = 'None',
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0))


ppn <- p11/p12 + plot_layout(height=c(6,1))
ppn
ggsave(filename = "fig3_a2.svg", plot = ppn, width = 4, height = 4, units = "in", dpi = 300)



#=======================================================================================================
meta <- read.csv('data/meta.csv',sep=',')
rownames(meta) <- meta$Sample_ID
meta$label <- meta$Sample_ID
meta$Country[meta$Country %in% c("Egypt","Kenya","Tunisia")] <- 'Africa'
meta$Country[meta$Country %in% c("Denmark","Germany","Spain","France","UK")] <- 'Europe'
meta$Country[meta$Country %in% c("USA","Guatemala")] <- 'America'
meta$Country <- factor(meta$Country, levels = c('China','Taiwan','Japan','Korea','America','Europe','Africa'))

meta <- meta[meta$ST == '3',]
tree <- read.tree("data/P1_type1/ref_M129.node_labelled.final_tree.tre")
tree <- drop.tip(tree, setdiff(tree$tip.label,rownames(meta)))

y <- full_join(as_tibble(tree),as_tibble(meta), by = 'label')
p21 <- ggtree(tree, right = F) %<+% y + geom_tippoint(aes(color=Country), size=1) +
  scale_color_brewer(palette = 'Set1') + layout_dendrogram()

tip_df <- p21$data[p21$data$isTip, ]
tip_df$x <- tip_df$x*0.975
tip_df$lt <- tip_df$Subclude
tip_df$lt[is.na(tip_df$lt)] <- "other ST3"
tip_df$lt[tip_df$lt != 'other ST3'] <- 'hited'
p21 <- p21 + geom_segment(data = tip_df, aes(x = x, xend = max(x), y = y, yend = y, linetype=lt), linewidth = 0.2, color = "grey40") +
  theme(plot.margin = margin(0, 0, -5, 0), legend.position = 'top') + guides(linetype = "none")



get_taxa_name(p21) -> strain_name
supplementary <- meta[strain_name,]
supplementary$Sample_ID <-factor(supplementary$Sample_ID, levels = rev(strain_name), ordered = T)

supplementary$ICU_ <- ifelse(is.na(supplementary$ICU), 'na.', supplementary$ICU)
supplementary$ICU_[supplementary$ICU_ == '1'] <- 'ICU'
supplementary$ICU_[supplementary$ICU_ == '0'] <- 'no'
supplementary$ICU_ <- factor(supplementary$ICU_, levels = c('ICU','no','na.'))

supplementary$CNS_ <- ifelse(is.na(supplementary$CNS),'na.',supplementary$CNS)
supplementary$CNS_[supplementary$CNS_ == '1'] <- 'CNS'
supplementary$CNS_[supplementary$CNS_ == '0'] <- 'no'
supplementary$CNS_ <- factor(supplementary$CNS_, levels = c('CNS','no','na.'))

supplementary$Subclude[is.na(supplementary$Subclude)] <- 'other ST3'
supplementary$Subclude <- factor(supplementary$Subclude, levels = c('subclade_pre','subclade_post_a','subclade_post_b','subclade_post_c','subclade_post_d','other ST3'))

suppl1 <- supplementary[,c('Sample_ID','ICU_','CNS_','Subclude','Time','Country')]
suppl1$Mono <- suppl1$Subclude
suppl1$icu <- 'ICU'
suppl1$cns <- 'CNS'
suppl1$subclude <- 'Mono'
suppl1$time <- 'Time'
suppl1$Time <- factor(suppl1$Time, levels = c('< 2010','< 2019','2019','2020','2023','2024 Jan','2024 Apr','2024 Jul','2024 Oct'))

p22 <- ggplot(data = suppl1) +
  geom_tile(aes(x = Sample_ID, y = cns, fill = CNS_), color = NA, alpha = 1, linewidth = 0.5) +
  scale_fill_manual(values = c('CNS'='#FF0000','no'='#33CCFF','na.'='#AAAAAA'),
                    name='CNS', guide = guide_legend(order = 4, nrow = 2, byrow = TRUE)) + new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = icu, fill = ICU_), color = NA, alpha = 1, linewidth = 0.5) +
  scale_fill_manual(values = c('ICU'='#FFAA33','no'='#33CCFF','na.'='#AAAAAA'),
                    name='ICU', guide = guide_legend(order = 3, nrow = 2, byrow = TRUE)) + new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = time, fill = Time), color = NA, alpha = 1, linewidth = 0.5) +
  scale_fill_manual(values = c('< 2010'='#AAAAFF','< 2019'='#5555FF','2019'='#0000FF','2020'='#FFDADA',
                               '2023'='#FFB6B6','2024 Jan'='#FF9191','2024 Apr'='#FF6D6D','2024 Jul'='#FF2424','2024 Oct'='#FF0000'),
                    name='Time', guide = guide_legend(order = 2)) + new_scale_fill() +
  
  geom_tile(aes(x = Sample_ID, y = subclude, fill = Subclude), color = NA, alpha = 1, linewidth = 0.5) +
  scale_fill_manual(values = c('subclade_post_a'='#8DD3C7','subclade_post_b'='#FFFFB3','subclade_post_c'='#BEBADA',
                               'subclade_post_d'='#FB8072','subclade_pre'='#80B1D3','other ST3'='#AAAAAA'),
                    name='Mono', guide = guide_legend(order = 1)) +
  
  labs(x = "", y = "")  + scale_y_discrete(limits = c("CNS", "ICU", "Time","Mono")) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

ppn <- p21/p22 + plot_layout(height=c(6,1))
ggsave(filename = "fig3_b.svg", plot = ppn, width = 13, height = 5.5, units = "in", dpi = 300)
