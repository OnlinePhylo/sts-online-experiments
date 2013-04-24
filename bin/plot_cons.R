#!/usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(TRUE)
stopifnot(length(args) == 3)
input <- args[1]
out_l2 <- args[2]
out_rf <- args[3]

theme_set(theme_bw(16))
theme_update(legend.position='bottom')

cons <- read.csv(input, as.is=TRUE)
cons <- transform(cons, tree=basename(tree),
                  type=ifelse(grepl('output/[^/]+/[^/]+.sum.tre', query), 'MrBayes', 'sts-online'),
                  n_taxa_label=paste(n_taxa, 'taxa'),
                  tree_num=as.integer(gsub('.*taxon-(\\d+).*', '\\1', tree)),
                  trim_count_label=sprintf('%02d trimmed', trim_count))
cons <- subset(cons, type != 'MrBayes' | particle_factor == 1)
cons <- transform(cons, particle_factor=ifelse(type == 'MrBayes', 0, particle_factor))

p <- ggplot(cons, aes(x=ordered(particle_factor), y=euclidean_distance, color=type)) +
    facet_grid(n_taxa_label~trim_count_label, scales='free_y') +
    geom_point(aes(shape=factor(tree_num))) +
    xlab("Particle Factor (x number of trees in posterior)") +
    ylab("Branch Length Distance (L2)")
svg(out_l2, width=10, height=7)
print(p)
dev.off()

p <- ggplot(cons, aes(x=ordered(particle_factor), color=type, y=rf_distance)) +
    facet_grid(n_taxa_label~trim_count_label, scales='free_y') +
    geom_point(aes(shape=factor(tree_num))) +
    xlab("Particle Factor (x number of trees in posterior)") +
    ylab("Robinson Foulds Distance")
svg(out_rf, width=10, height=7)
print(p)
dev.off()
