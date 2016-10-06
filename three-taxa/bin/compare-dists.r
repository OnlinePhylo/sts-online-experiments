#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(jsonlite))
#theme_set(theme_bw(16))

args <- commandArgs(TRUE)

stopifnot(length(args) == 4)

sts_output_path <- args[1]
mb_output_path <- args[2]
empirical_output_path <- args[3]
outfile <- args[4]

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
mb_length <- read.delim(mb_output_path, as.is = TRUE, skip = 1, sep = '\t')
mb_length <- transform(mb_length,
                       length = TL,
                       weight = rep(1, nrow(mb_length)),
                       type = 'mb')
sel <- c('length', 'weight', 'type')
both <- rbind(sts_length[, sel], empirical_length[, sel], mb_length[, sel])
both <- group_by(both, type) %>%
  mutate(weight = weight / sum(weight)) %>%
  ungroup()

write.csv(both, outfile, row.names=FALSE)
