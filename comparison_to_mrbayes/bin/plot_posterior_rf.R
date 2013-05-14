#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(plyr)
theme_set(theme_bw(16))
theme_update(axis.text.x=element_text(angle=90))

args <- commandArgs(TRUE)
stopifnot(length(args) == 3)
pc <- read.csv(args[1], as.is=TRUE)
pc <- transform(pc, tree=sub('\\.nwk$', '', basename(tree)))
pc <- transform(pc, tree_number=as.integer(sub('\\d+taxon-(\\d+)', '\\1', tree)))
pc <- transform(pc, tree_label=paste('tree', tree_number),
                n_taxa_label=paste(n_taxa, 'taxa'), n_trimmed=ifelse(trim_taxon=='', 0, nchar(gsub('[^-]', '', trim_taxon)) + 1))

tt <- unique(pc[, c('trim_taxon', 'n_trimmed')])
tt <- tt[order(tt$n_trimmed, tt$trim_taxon),]
pc <- transform(pc, trim_taxon=factor(trim_taxon, levels=tt$trim_taxon))

#print(aggregate(rf_distance~file, pc, length))

m <- melt(pc, measure.vars=c('rf_distance', 'weighted_rf', 'euclidean'))

measure_names <- data.frame(variable=c('rf_distance', 'weighted_rf', 'euclidean'),
                            measure=c('RF Distance', 'L2', 'L1'))

m <- transform(m, measure=measure_names$measure[match(as.character(m$variable), measure_names$variable)])

pdf(args[2], width=11, height=11)
# One plot per measure
d_ply(m, .(measure), function(piece) {
  measure <- piece$measure[1]
  message(measure)
  p <- ggplot(piece, aes(x=trim_taxon, y=value, fill=type, weight=exp(log_weight))) +
      geom_boxplot() +
      facet_grid(tree_label~n_taxa_label, scales='free_x') +
      theme(legend.position='bottom') +
      xlab('Trimmed taxon') +
      ylab(measure)
  print(p)
})
dev.off()

pdf(args[3], width=7, height=7)
# One plot per measure, tree size
d_ply(m, .(measure, n_taxa), function(piece) {
  measure <- piece$measure[1]
  n_taxa_label <- piece$n_taxa_label[1]
  message(paste(n_taxa_label, measure))
  p <- ggplot(piece, aes(x=trim_taxon, y=value, fill=type, weight=exp(log_weight))) +
      geom_boxplot() +
      facet_wrap(~tree_label, scales='free_x') +
      theme(legend.position='bottom') +
      xlab('Trimmed taxon') +
      ylab(measure) +
      ggtitle(n_taxa_label)
  print(p)
})

dev.off()
