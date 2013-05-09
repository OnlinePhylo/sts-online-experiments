#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)
library(reshape2)

args <- commandArgs(TRUE)
stopifnot(length(args) == 3)
input <- args[1]
out_l2 <- args[2]
out_rf <- args[3]

theme_set(theme_bw(base_size=16))
theme_update(legend.position='bottom')

cons <- read.csv(input, as.is=TRUE)
cons <- transform(cons, tree=basename(tree),
                  type=ifelse(grepl('output/[^/]+/[^/]+.sum.tre', query), 'MrBayes', 'sts-online'),
                  n_taxa_label=paste(n_taxa, 'taxa'),
                  tree_num=as.integer(gsub('.*taxon-(\\d+).*', '\\1', tree)),
                  trim_count_label=sprintf('%02d trimmed', trim_count))
cons <- subset(cons, type != 'MrBayes' | particle_factor == 1)
cons <- transform(cons, particle_factor=ifelse(type == 'MrBayes', 0, particle_factor))

m <- melt(cons, measure.vars=c('euclidean_distance', 'rf_distance'))

y_labels <- data.frame(variable=c('euclidean_distance', 'rf_distance'),
                       label=c('Branch Length Distance (L2)', 'RF Distance'),
                       path=c(out_l2, out_rf))

d_ply(m, .(variable), function(piece) {
   s <- ddply(piece, .(particle_factor, tree_num, n_taxa_label, type, trim_taxon), function(p) {
      with(p, data.frame(q25=quantile(value, 0.25), min_value=min(value), max_value=max(value), median_value=median(value), q75=quantile(value, 0.75)))
   })

   m <- y_labels[match(piece$variable[1], y_labels$variable),]
   p <- ggplot(s, aes(x=ordered(particle_factor), y=median_value, ymin=q25, ymax=q75, color=type)) +
     facet_grid(tree_num~n_taxa_label, scales='free_y') +
     geom_pointrange(position=position_jitter()) +
     xlab("Particle Factor (x number of trees in posterior)") +
     ylab(m$label)
   message(m$path)
   svg(as.character(m$path), width=6, height=11)
   print(p)
   dev.off()
})
