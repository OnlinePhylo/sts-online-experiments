#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_bw(16) + theme(strip.background = element_blank()))

args <- commandArgs(TRUE)

stopifnot(length(args) == 2)

input_path <- args[1]
output_path <- args[2]

comp <- read.csv(input_path, as.is=TRUE)

ml_bl <- group_by(comp, proposal_method, seed) %>%
  filter(type == 'empirical') %>%
  arrange(-weight) %>%
  do(head(., 1)) %>%
  ungroup

pdf(output_path, width = 30, height = 10, useDingbats = FALSE)
tryCatch({
  p <- ggplot(comp, aes(x = length, weight = weight)) +
    geom_density(aes(color = type)) +
    geom_vline(aes(xintercept = length), linetype = 'dashed', data = ml_bl) +
    facet_grid(proposal_method ~ seed, scales= 'free_y') +
    theme(legend.position = 'bottom') +
    xlab("Branch length") +
    ggtitle(paste("branch length ", unique(comp$branch_length)))
  if(all(comp$branch_length < 0.05))
    p <- p + xlim(0, 0.1)
    #ggtitle(sprintf("sim_bl=%.2f seed=%d proposal=%s", control$branch_length, control$seed, control$proposal_method)) +
  else
    p <- p + xlim(0, 0.4)
  print(p)
}, finally = dev.off())
