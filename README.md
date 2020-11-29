# characterization_ERA-5_daily_precipitation
This repository contains the R-codes used for the comparison of intensity precipitation from ERA-5, and EOBS and CMORPH.
They are based on the Extended Generalized Pareto distribution (EGPD), Tencaliec et al. (2019).


quantile_confidence_interval.R   contains a function to extract bootstrapped EGPD parameters and quantiles.
KB_test_functions.R              contains functions to perform the Kullback-Leibler test.
