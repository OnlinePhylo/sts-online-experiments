Before I kick off the 50 taxon runs, RFC on process:

1. Simulate 5 50-taxon trees with branch lengths drawn from exponential with mean 0.1
2. Simulate 1000 site alignment for each tree
3. Re-fit branch lengths using PhyML
4. Run MrBayes on each alignment (4 runs)
5. Prune 1 taxon from each alignment, run MrBayes followed by STS (repeated for 5 different taxa)
6. Compare
