#!/usr/bin/env Rscript

library(assertthat)
library(TreeSim)

set.seed(1)

n_trees_per_taxon_count <- 5
n_taxa <- c(10, 50, 100)

r <- lapply(as.list(n_taxa), function(n) {
  trees <- sim.bd.taxa(n = n,
                       numbsim = n_trees_per_taxon_count,
                       mu = 2.0,
                       lambda = 6.0,
                       frac = 0.7,
                       complete = FALSE,
                       stochsampling = TRUE)[[1]]
  assert_that(is.list(trees))
  assert_that(length(trees) == n_trees_per_taxon_count)
  lapply(seq_along(trees), function(i) {
    path <- sprintf('%02dtaxon-%02d.nwk', n, i)
    message(path)
    write.tree(trees[[i]], path)
    path
  })
})
