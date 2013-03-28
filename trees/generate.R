#!/usr/bin/env Rscript

library(ape)

set.seed(1)

n_trees <- 10
n_taxa <- 10

idx <- 1:n_trees

r <- lapply(as.list(idx), function(i) {
  path <- sprintf('tree%03d.nwk', i)
  write.tree(rtree(n_taxa, br=rexp, rate=10),
             path)
  message(path)
  path
})
