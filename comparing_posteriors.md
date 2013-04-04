# For each tree:

1. Simulate an alignment
2. Run MrBayes on full alignment (multiple runs, long chains?)
3. Vary:

    - Number of particles
    - Number of taxa trimmed
    - Number of MCMC moves (?)

For each combination:

1. Trim the tree
2. Run MrBayes
3. Run `sts-online`
4. Compare sample from `sts-online` to reference run on full set of taxa using

     - RF distance
     - Branch lengths?
     - Log-likelihood?
