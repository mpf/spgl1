# User's guide

## Lasso

The standard Lasso problem, given by

$$\mathop{\mathrm{minimize}}_{x}\quad {\textstyle\frac{1}{2}}\Vert Ax-b\Vert_2^2\quad\mathrm{subject\ to}\quad \Vert x\Vert_1 \leq \tau,$$

can be solved by SPGL1 using the ``spg_lasso`` function:
```matlab
[x,r,g,info] = spg_lasso(A,b,tau,options)
```

The ``options`` parameter controls the stopping criteria and behavior of the solver and consists of keyword-value pairs or structures (see [solver options](options.md) for more details). Default options are used when no options are specified.

As an example, let's create random matrix $A$ and vector $b$, and solve the
Lasso problem for $\tau=0.5$ using:

```matlab 
s = RandStream('mt19937ar','Seed',0);
A = randn(s,5,10);
b = randn(s,5,1);
[x,r,g,info] = spg_lasso(A,b,0.5);
```

The SPGL1 solver gives an output similar to this:

```text
 ===================================================================
 SPGL1  v2.0 (28 Nov 2019)
 ===================================================================
 No. rows            :        5       Initial tau         : 5.00e-01
 No. columns         :       10       2-norm of b         : 1.86e+00
 Optimality tol      : 1.00e-04       bound on 1-norm of x: 5.00e-01
 Basis pursuit tol   : 1.00e-04       Maximum iterations  :       50

  Iter     ||Ax-b||_2   Relative Gap  ||A'r||_d
     0  1.8616905e+00  9.7024931e-01   3.36e+00
     1  1.4495790e+00  9.5092862e-01   2.74e+00
     2  1.3023687e+00  6.1655207e-01   2.30e+00
     3  1.2009939e+00  1.7815336e-01   1.43e+00
     4  1.1775752e+00  1.5030186e-01   1.55e+00
     5  1.1558273e+00  1.2492861e-01   1.34e+00
     6  1.1513362e+00  7.1502516e-02   1.28e+00
     7  1.1492592e+00  4.1721960e-02   1.22e+00
     8  1.1465284e+00  3.8571207e-02   1.28e+00
     9  1.1491912e+00  4.1627659e-02   1.47e+00
    10  1.1486180e+00  4.0969188e-02   1.44e+00
    11  1.1446930e+00  4.3987941e-03   1.25e+00
    12  1.1446779e+00  4.0502125e-03   1.25e+00
    13  1.1446237e+00  1.5756440e-04   1.25e+00
    14  1.1446236e+00  2.5627973e-05   1.25e+00

 EXIT -- Optimal solution found

 Products with A     :      14        Total time   (secs) :      0.1
 Products with A'    :      15        Project time (secs) :      0.0
 Newton iterations   :       0        Mat-vec time (secs) :      0.0
 Line search its     :       0
```

The exit status indicates that an optimal solution was found. The output of the solver is given by the tuple ``[x,r,g,info]`` with the following fields: 

* the computed solution ``x``
* the residual ``r = b-Ax``
* the gradient ``g := -A'r`` of the objective at the solution
* an ``info`` structure that contains information such as the exit status, the number of iterations, and the total run time; see [Information Structure](infostruct.md) for details.

Continuing with our example, we find that if we solve Lasso with $\tau=1$, the required number of iterations exceed the default of 10 times the number of rows in $A$. We can increase the maximum number of iterations to 100 by calling
```matlab
[x,r,g,info] = spg_lasso(A,b,1,'iterations',100);
```

In this case around 90 iterations suffice to find an optimal solution.


## Basis-pursuit denoise

The standard basis-pursuit denoise problem is given by

$$\mathop{\mathrm{minimize}}_{x}\quad  \Vert x\Vert_1\quad \mathrm{subject\ to}\quad \Vert Ax-b\Vert_2 \leq \sigma.$$

Solving the problem in SPGL1 can be done by calling the ``spg_bpdn`` function

```matlab
[x,r,g,info] = spg_bpdn(A,b,sigma,options)
```

In the special case where $\sigma=0$, the problem reduces to the basis pursuit problem

$$\mathop{\mathrm{minimize}}_{x}\quad  \Vert x\Vert_1\quad \mathrm{subject\ to}\quad Ax=b,$$

which can be solved using the basis-pursuit solver:

```matlab
    [x,r,g,info] = spg_bp(A,b,options)
```

The basis-pursuit solution is obtained by solving a sequence of Lasso problems with suitable values of $\tau$. The gradient output ``g`` thus refers to the objective value of the Lasso subproblem corresponding to the latest value of $\tau$. As an example, consider the same example used for the Lasso problem, now solved with basis-pursuit denoise for $\sigma=1.0$:

```matlab 
s = RandStream('mt19937ar','Seed',0);
A = randn(s,5,10);
b = randn(s,5,1);
[x,r,g,info] = spg_bpdn(A,b,1.0);
```

The solver gives an output similar to

```text
 ===================================================================
 SPGL1  v2.0 (28 Nov 2019)
 ===================================================================
 No. rows            :        5       Initial tau         : 0.00e+00
 No. columns         :       10       2-norm of b         : 1.86e+00
 Optimality tol      : 1.00e-04       Target ||Ax-b||_2   : 1.00e+00
 Basis pursuit tol   : 1.00e-04       Maximum iterations  :       50

  Iter     ||Ax-b||_2   Relative Gap  Rel Error  ||A'r||_d   ||x||_p
     0  1.8616905e+00  0.0000000e+00   4.63e-01   3.36e+00  4.77e-01
     1  1.4973610e+00  5.2616060e-01   3.32e-01   2.12e+00
     2  1.2974269e+00  1.9995584e-01   2.29e-01   1.62e+00
     3  1.2030579e+00  8.1971765e-02   1.69e-01   1.77e+00
     4  1.1995095e+00  7.7709077e-02   1.66e-01   1.91e+00  6.25e-01
     5  1.0272341e+00  8.8274756e-02   2.65e-02   1.11e+00
     6  1.0212485e+00  8.2144096e-02   2.08e-02   1.09e+00
     7  1.0187911e+00  2.8385712e-02   1.84e-02   9.66e-01
     8  1.0183832e+00  2.2393424e-02   1.81e-02   9.64e-01  6.45e-01
     9  9.9870605e-01  1.6574628e-02   1.29e-03   9.45e-01
    10  1.0008367e+00  1.8704785e-02   8.36e-04   1.11e+00
    11  1.0070820e+00  2.4974860e-02   7.03e-03   1.30e+00
    12  9.9794627e-01  1.5816119e-02   2.05e-03   1.01e+00
    13  9.9785653e-01  9.7287655e-03   2.14e-03   9.87e-01  6.42e-01
    14  9.9993939e-01  9.6116065e-03   6.06e-05   9.92e-01

 EXIT -- Found a root

 Products with A     :      18        Total time   (secs) :      0.1
 Products with A'    :      16        Project time (secs) :      0.0
 Newton iterations   :       4        Mat-vec time (secs) :      0.0
 Line search its     :       3
```
  
The last column of the output shows that four Lasso subproblems with different values of $\tau$ are solved to obtain the basis-pursuit denoise solution.
 
## Group norms

SPGL1 supports group-sparse versions of the three main problem classes (Lasso, basis-pursuit denoise, and basis pursuit). The *group norm* is defined by given subsets of the vector elements, such that each element occurs in exactly one subset. The norm is then define as the sum of the Euclidean norms of the subvectors in each set. To solve the group-norm basis-pursuit denoise formulation we can use
 
    [x,r,g,info] = spg_group(A,b,groups,sigma,options)


The ``groups`` parameter is a vector that contains the group number for each of the elements in `x`. The group numbers themselves can be chosen arbitrary, as long as elements in the same group have the same number. Groups do not have to consist of contiguous elements, although in practice they often are.

In the following example we use three groups, labeled `1`, `2`, and `3`:

```matlab
s = RandStream('mt19937ar','Seed',0);
A = randn(s,5,10);
b = randn(s,5,1);
sigma = 1.2
x = spg_group(A,b,[1,1,1,2,2,2,2,3,3,3],sigma);
```

In this case we get a result that is group-sparse: the elements in some groups are all zero, whereas those in other groups are all non-zero:

```text
  x =

         0
         0
         0
   -0.1804
   -0.2038
   -0.1107
    0.0037
         0
         0
         0
```

## Multiple measurement vectors (MMV)

A special case of group sparsity is the multiple-measurement vectors (MMV) problem. In this problem we are given a matrix of measurements `B`, and assume that the unknown matrix `X` is such that all columns have the same support. This is often achieved by minimizing the sum of Euclidean norms of the rows in `X`. This can be reformulated as a group-norm problem by appropriately vectorizing `B` and `X`, and suitably redefining `A`. 

For convenience SPGL1 provides the function

    [X,R,G,info] = spg_mmv(A,B,sigma,options);

As an example, we solve

```matlab
s = RandStream('mt19937ar','Seed',0);
A = randn(s,5,7);
B = randn(s,5,6);
[X,R,G,info] = spg_mmv(A,B,3.5);
```

which gives

```text
 X =
 
    0.1102   -0.0296   -0.0177    0.0365   -0.0298   -0.0218
    0.0123    0.0042    0.1348   -0.0577   -0.0208    0.1842
   -0.0020   -0.0001   -0.0008    0.0009   -0.0087    0.0027
         0         0         0         0         0         0
    0.0806    0.0252    0.0703   -0.1165    0.0985    0.0176
         0         0         0         0         0         0
    0.2815    0.1036   -0.1306    0.0197   -0.0596   -0.2179
```

## Generic interface

The generic interface to SPGL1 is given by

```matlab
[x,r,g,info] = spgl1(A,b,tau,sigma,x0,options)
[x,r,g,info] = spgl1(A,b,tau,sigma,options)
[x,r,g,info] = spgl1(A,b,tau,options)
```

* The ``options`` parameters are optional and can be a mixture of structure objects and key-value pairs. The first option parameter can be a string identifying a predefined parameter set, which can then be modified by the parameters that follow. See the examples in [options](options.md).

* The ``x0`` parameter can be provided to initialize $x$. If set, the parameters ``tau`` and ``sigma`` must also be provided.

* To solve the Lasso problem formulation with an initial value for ``x``, set ``sigma`` to the empty vector ``[]``. In order to solve basis pursuit denoise we can set ``tau`` to ``0`` or an empty vector ``[]`` (strongly recommended), or provide an initial value for ``tau`` (generally not recommended). In case an initial value of ``tau`` is specified, it is important that this value be smaller than the value for tau at the solution; reduction of ``tau`` from a value that is too large can be very time consuming or result in a suboptimal basis-pursuit solution. When the ``x0`` parameter is set to the empty vector ``[]``, a default initial value for $x$ is used.



## Custom norms

SPGL1 can be extended to solve Lasso and basis-pursuit denoise problems with
custom primal norms. For this we need to provide three functions

1. ``options.primal_norm``, which evaluates the primal norm $\|x\|_p$ of a given vector $x$;
2. ``options.dual_norm``, which evaluates the dual norm $\|y\|_d$ corresponding to the primal norm;
3. ``options.project``, which evaluates orthogonal projection onto a scaled unit-norm ball corresponding to the primal norm, i.e., onto the set $\mathcal{B}_p:=\{x \mid \|x\|_p\le\tau\}$ for some positive value $\tau$.

Good examples of the definition of these norms can be found in the ``spg_group`` solver.


