#!/usr/bin/env Rscript

library(ape)

set.seed(1)

n_trees_per_taxon_count <- 5
n_taxa <- c(10, 50)

r <- lapply(as.list(n_taxa), function(i) {
  lapply(as.list(1:n_trees_per_taxon_count), function(j) {
    path <- sprintf('%02dtaxon-%02d.nwk', i, j)
    write.tree(rtree(i, br=rexp, rate=10),
               path)
    message(path)
    path
  })
})
