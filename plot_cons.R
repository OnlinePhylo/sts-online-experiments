#!/usr/bin/env Rscript
library(ggplot2)

theme_set(theme_bw(16))

cons <- read.csv('consensus_to_source.csv', as.is=TRUE)
cons <- transform(cons, tree=basename(tree),
                  type=ifelse(grepl('output/[^/]+/[^/]+.sum.tre', query), 'MrBayes', 'sts-online'),
                  n_taxa_label=paste(n_taxa, 'taxa'),
                  trim_count_label=sprintf('%02d trimmed', trim_count))
cons <- subset(cons, type != 'MrBayes' | particle_factor == 1)
cons <- transform(cons, particle_factor=ifelse(type == 'MrBayes', 0, particle_factor))

p <- ggplot(cons, aes(x=ordered(particle_factor), y=euclidean_distance, fill=type)) +
    facet_grid(n_taxa_label~trim_count_label, scales='free_y') +
    geom_boxplot() +
    xlab("Particle Factor (x number of trees in posterior)") +
    ylab("Branch Length Distance (L2)")
svg('consensus_to_source.svg', width=10, height=7)
print(p)
dev.off()

p <- ggplot(cons, aes(x=ordered(particle_factor),fill=type, y=rf_distance)) +
    facet_grid(n_taxa_label~trim_count_label, scales='free_y') +
    geom_boxplot() +
    xlab("Particle Factor (x number of trees in posterior)") +
    ylab("Robinson Foulds Distance")
svg('consensus_to_source_rf.svg', width=10, height=7)
print(p)
dev.off()
