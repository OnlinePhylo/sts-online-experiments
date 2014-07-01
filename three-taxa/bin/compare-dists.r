#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(jsonlite))
#theme_set(theme_bw(16))

args <- commandArgs(TRUE)

stopifnot(length(args) == 3)

sts_output_path <- args[1]
empirical_output_path <- args[2]
outfile <- args[3]

sts_doc <- fromJSON(sts_output_path)
control <- fromJSON(file.path(dirname(sts_output_path), 'control.json'))
sts_length <- data.frame(length=sts_doc$trees[, 'treeLength'],
                         weight=rep(1, length(sts_doc$trees)),
                         type='sts')
empirical_length <- read.csv(empirical_output_path, as.is = TRUE)
empirical_length <- transform(empirical_length,
                              length = branch_length,
                              weight = exp(posterior - max(posterior)),
                              type = 'empirical')
sel <- c('length', 'weight', 'type')
both <- rbind(sts_length[, sel], empirical_length[, sel])
both <- group_by(both, type) %>%
  mutate(weight = weight / sum(weight)) %>%
  ungroup()

write.csv(both, outfile, row.names=FALSE)
#pdf(outfile, width = 7, height = 5, useDingbats = FALSE)
#tryCatch({
  #p <- ggplot(both, aes(x = length, weight = weight, color=type)) +
    #geom_density() +
    #geom_vline(xintercept = control$branch_length, linetype = 'dashed') +
    #ggtitle(sprintf("sim_bl=%.2f seed=%d proposal=%s", control$branch_length, control$seed, control$proposal_method)) +
    #xlab("Branch length")
  #print(p)
#}, finally = dev.off())
