# Information structure

```text
===========  ================================================
Field        Description
===========  ================================================
tau          Value of :math:`\tau` in final Lasso subproblem
rNorm        Norm of the final residual
gNorm        Norm of the final gradient vector
rGap         Relative duality gap
stat         Exit status (see below)
success      Exit status (see below)
statusStr    Exit status (see below)
iter         Number of iterations
nProdA       Number of products with A
nProdAt      Number of products with A transpose
nNewton      Number of Newton root-finding iterations
timeProject  Time spent in the projection function
timeMatProd  Time spent in matrix-vector products
timeTotal    Total runtime
options      Options provided to the solver
xNorm1       History of the primal norm of x (optional)
rNorm2       History of the residual norm (optional)
lambda       History of the gradient norm (optional)
===========  ================================================
```
The history fields are updated every iteration if ``options.history`` is set to true.
