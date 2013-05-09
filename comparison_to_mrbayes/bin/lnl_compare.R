#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(TRUE)

stopifnot(length(args) == 5)

cols <- c('LnL', 'source', 'file_name')

read_sts <- function(s) {
  sts <- read.table(s, col.names=c('LnL', 'tree'))
  transform(sts, source='STS-online', file_name=basename(s))[,cols]
}

read_mb <- function(s) {
  mb <- read.table(s, header=TRUE, skip=1)
  # Remove 25% as burn-in
  mb <- mb[as.integer(nrow(mb) / 4):nrow(mb),]
  transform(mb, source='MrBayes', file_name=basename(s))[,cols]
}

stacked <- rbind(read_sts(args[1]), read_sts(args[2]),
                 read_mb(args[3]), read_mb(args[4]))


p1 <- ggplot(stacked, aes(x=LnL, linetype=source, color=file_name)) +
  geom_density()
p2 <- ggplot(stacked, aes(x=file_name, y=LnL, fill=source)) +
  geom_boxplot()

pdf(args[5])
print(p1)
print(p2)
dev.off()
