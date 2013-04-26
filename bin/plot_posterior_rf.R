#!/usr/bin/env Rscript
library(ggplot2)
theme_set(theme_bw(16))

args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
pc <- read.csv(args[1], as.is=TRUE)
pc <- transform(pc, tree=sub('\\.nwk$', '', basename(tree)))

p <- ggplot(pc, aes(x=trim_taxon, y=rf_distance, fill=type)) +
    geom_boxplot() +
    facet_wrap(~tree, ncol=1) +
    theme(legend.position='bottom') +
    xlab('Trimmed taxon') +
    ylab('RF Distance (unweighted)')

ggsave(args[2], height=11, width=5)
