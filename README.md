# SPGL1: A solver for large-scale sparse least squares

[![GitHub license](https://img.shields.io/github/license/mpf/spgl1)](https://github.com/mpf/spgl1/blob/master/COPYING)
[![DOI:10.1137/080714488](https://zenodo.org/badge/DOI/10.1137/080714488.svg)](https://doi.org/10.1137/080714488)

* **Documentation**: https://friedlander.io/spgl1

## Introduction

SPGL1 is a Matlab solver
for large-scale one-norm regularized least squares.  It is designed to
solve any of the following three problems:

1. Basis pursuit denoise (BPDN):
   minimize  ||x||_1  subject to  ||Ax - b||_2 <= sigma,

2. Basis pursuit (BP):
   minimize   ||x||_1  subject to  Ax = b
 
3. Lasso:
   minimize  ||Ax - b||_2  subject to  ||x||_1 <= tau,

The matrix A can be defined explicily, or as an operator (i.e., a
function) that return both both Ax and A'y.  SPGL1 can solve these
three problems in both the real and complex domains.


Home page: https://friedlander.io/spgl1


## References :notebook:

The algorithm implemented by SPGL1 is described in the paper

- E. van den Berg and M. P. Friedlander, "Probing the Pareto frontier
  for basis pursuit solutions", SIAM J. on Scientific Computing,
  31(2):890-912, November 2008
