#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)

comparisons <- read.csv('comparisons.csv')
comparisons <- transform(comparisons,
                         is_sts1=grepl('_online', tree1),
                         is_sts2=grepl('_online', tree2),
                         tree=sub('output/(tree\\d{3})/.*$', '\\1', tree1))
comparisons <- transform(comparisons, type=ifelse(is_sts1 & is_sts2, 'within-online',
                                        ifelse(!is_sts1 & !is_sts2, 'within-mrbayes', 'between')))
head(comparisons)

print(ddply(comparisons, .(type), function(x) summary(x$d12_mean)))
