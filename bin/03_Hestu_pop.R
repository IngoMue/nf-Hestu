#!/usr/bin/env Rscript
#Plot Heterozygosity with population information
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(ggpubr)

#Individual Het
sfslist <- args[2:length(args)]

het_df <- data.frame(matrix(nrow = length(sfslist), ncol = 2))
colnames(het_df) <- c("ind", "het")

for (i in 1:length(sfslist)) {
  a <- scan(sfslist[i])
  het <- a[2]/sum(a)
  
  het_df$ind[i] <- substr(sfslist[i], 1, nchar(sfslist[i])-4)
  het_df$het[i] <- het
}

write.table(het_df, file = "Ind_het.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

het_df$ind <- factor(het_df$ind, levels = het_df$ind[order(het_df$het)])

ind_het <- ggplot(het_df, aes(x = ind, y = het)) +
  geom_bar(stat = "identity", fill = "forestgreen", colour = "black") +
  labs(x = "Individual", y = "Heterozygosity") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("Ind_het.png", ind_het, device = "png", width = 15, height = 10, dpi = 300)


#Heterozygosity against coverage
DoC_df <- read.table(args[1], header = T)

Het_DoC_df <- left_join(het_df, DoC_df, by = "ind")

Het_DoC <- ggplot(Het_DoC_df, aes(x = DoC, y = het)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", colour = "red", fill = "firebrick", linewidth = 0.5) +
  stat_cor(method = "kendall") +
  labs(x = "Depth of Coverage [X]", y = "Heterozygosity")

ggsave("Het_DoC.png", Het_DoC, device = "png", width = 15, height = 10, dpi = 300)

#the same plot but with population information
pop_het <- ggplot(Het_DoC_df, aes(x = pop, y = het)) +
  geom_boxplot(aes(fill = pop)) +
  labs(x = "Population", y = "Heterozygosity") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "null")

ggsave("Pop_het.png", pop_het, device = "png", width = 15, height = 10, dpi = 300)

Het_DoC_pop <- ggplot(Het_DoC_df, aes(x = DoC, y = het, fill = pop)) +
  geom_point(size = 2, pch = 21) +
  geom_smooth(method = "lm", aes(colour = pop), linewidth = 0.5) +
  stat_cor(method = "kendall", aes(colour = pop)) +
  labs(x = "Depth of Coverage [X]", y = "Heterozygosity") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d()

ggsave("Het_DoC_pop.png", Het_DoC_pop, device = "png", width = 15, height = 10, dpi = 300)
