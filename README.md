# Compare-EvAM
Compare outputs from evolutionary accumulation models

Illustration of the ordering-matrix-based comparison approach for EvAM outputs from [1].

This code requires the `hyperinf` package for various EvAM approaches. Install with

`remotes::install_github("StochasticBiology/hyperinf")`

`comparison.R` generates synthetic test data -- both cross-sectional and simulated on a tree (using simulation code from HyperDAGs [2]) -- and visualises comparisons from the outputs. `cancer.R` uses data from [3] to compare inferred EvAM dynamics of chromosomal aberrations between tumour type pairs, and inference outputs from [4] to compare EvAM dynamics of drug resistance features in Klebsiella between countries.

[1] https://arxiv.org/abs/2608.02781

[2] https://github.com/StochasticBiology/hyperdags

[3] https://www.nature.com/articles/s41586-019-1907-7/

[4] https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3003848
