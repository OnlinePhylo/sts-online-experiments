#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(TRUE)

stopifnot(length(args) == 2)

infile <- args[1]
outfile <- args[2]

pp <- read.csv(infile, as.is=TRUE)
pp <- transform(pp, tree=sub('\\.nwk$', '', basename(tree)))

theme_set(theme_bw(16))

p <- ggplot(pp, aes(x=pp1, y=pp2, color=type, shape=type)) +
  geom_abline(slope=1, color='grey') +
  geom_point(alpha=0.7, size=1.4) +
  facet_grid(move_type~tree) +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Split posterior probability 1 (MrBayes)") +
  ylab("Split posterior probability 2") +
  theme(legend.position='bottom',
        axis.text.x=element_text(angle=35, hjust=1)) +
  ggtitle("Split posterior probability comparison")
ggsave(outfile, plot=p, width=11, height=6)
