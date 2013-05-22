#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

theme_set(theme_bw(16))

files <- list.files(pattern='50tax_trim_t1.*.comp.csv')
df <- do.call(rbind, lapply(as.list(files), read.csv, as.is=TRUE))
df <- transform(df, n_trees=particle_factor*750)

m <- melt(subset(df, particle_factor==1), measure.vars=c('rf_distance', 'weighted_rf', 'euclidean'))

names <- data.frame(variable=c('rf_distance', 'weighted_rf', 'euclidean'),
                    name=c('RF Distance', 'L2 distance', 'L1 distance'))
guided_labels <- data.frame(move_type=c('guided', 'noguided'),
                            move_label=c('Empirical Bayes', 'Simple'))

m <- transform(m, name=names$name[match(variable, names$variable)],
               move_label=guided_labels$move_label[match(move_type, guided_labels$move_type)])




p <- ggplot(m, aes(x=move_label, y=value, fill=move_label)) +
    geom_boxplot() +
    facet_wrap(~name, scales='free_y') +theme(legend.position='none')  +
    xlab("Proposal Mechanism")

ggsave('proposal_comparison.svg', width=8, height=4)
