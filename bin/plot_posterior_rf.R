#!/usr/bin/env Rscript
library(ggplot2)
theme_set(theme_bw(16))

args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
pc <- read.csv(args[1], as.is=TRUE)
pc <- transform(pc, tree=sub('\\.nwk$', '', basename(tree)))
pc <- transform(pc, tree_number=as.integer(sub('\\d+taxon-(\\d+)', '\\1', tree)))

#print(aggregate(rf_distance~file, pc, length))

message('RF')
p <- ggplot(pc, aes(x=trim_taxon, y=rf_distance, fill=type, weight=exp(log_weight))) +
    geom_boxplot() +
    facet_grid(tree_number~n_taxa, scales='free_x') +
    theme(legend.position='bottom') +
    xlab('Trimmed taxon') +
    ylab('RF Distance')

ggsave(args[2], height=11, width=8)

message('Weighted RF')
p <- ggplot(pc, aes(x=trim_taxon, y=weighted_rf, fill=type, weight=exp(log_weight))) +
    geom_boxplot() +
    facet_grid(tree_number~n_taxa, scales='free_x') +
    theme(legend.position='bottom') +
    xlab('Trimmed taxon') +
    ylab('RF distance (with branch lengths)')
ggsave(paste(args[2], '.weightedrf.svg', sep=''), height=11, width=8)

message('Euclidean')
p <- ggplot(pc, aes(x=trim_taxon, y=euclidean, fill=type, weight=exp(log_weight))) +
    geom_boxplot() +
    facet_grid(tree_number~n_taxa, scales='free_x') +
    theme(legend.position='bottom') +
    xlab('Trimmed taxon') +
    ylab('L2')
ggsave(paste(args[2], '.euclidean.svg', sep=''), height=11, width=8)
