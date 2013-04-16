#!/usr/bin/env Rscript
library(ggplot2)

theme_set(theme_bw(16))

cons <- read.csv('consensus_to_source.csv', as.is=TRUE)
cons <- transform(cons, tree=basename(tree),
                  type=ifelse(grepl('output/[^/]+/[^/]+.sum.tre', query), 'MrBayes', 'sts-online'),
                  n_taxa_label=paste(n_taxa, 'taxa'),
                  trim_count_label=sprintf('%02d trimmed', trim_count))

p <- ggplot(cons, aes(x=type, y=euclidean_distance)) +
    facet_grid(n_taxa_label~trim_count_label, scales='free_y') +
    geom_boxplot() +
    xlab("Consensus Tree Source") +
    ylab("Branch Length Distance (L2)")
svg('consensus_to_source.svg', width=10, height=7)
print(p)
dev.off()
