Introduction
==================

SPGL1 solves two problem formulations. The first formulation, the Lasso problem, is given by

.. math::

   \mathop{\mathrm{minimize}}_{x}\ {\textstyle\frac{1}{2}}\Vert Ax-b\Vert_2^2\quad\mathrm{subject\ to}\quad \Vert x\Vert_p \leq \tau

where :math:`\tau \geq 0` is given, and :math:`\Vert\cdot\Vert_p` denotes a suitable primal norm with corresponding dual :math:`\Vert\cdot\Vert_d`. The second formulation, the basis-pursuit denoise problem, is given by

.. math::

   \mathop{\mathrm{minimize}}_{x}\ \Vert x\Vert_p\quad \mathrm{subject\ to}\quad \Vert Ax-b\Vert_2 \leq \sigma

where :math:`\sigma` is the maximum two-norm of the residual.  The default primal norm is the :math:`\ell_1` norm given by sum of absolute values of the input. The corresponding `lasso`_ and `basis-pursuit denoise`_ solvers are described below. Several other norms included in SPGL1 are the group :math:`\ell_{1,2}` norm and the special multiple-measurement vector (MMV) case. For convenience, SPGL1 also supports :math:`\ell_2`  regularization of x.


.. _lasso:

Lasso
----------------

The standard Lasso problem, given by

.. math::

   \mathop{\mathrm{minimize}}_{x}\ {\textstyle\frac{1}{2}}\Vert Ax-b\Vert_2^2\quad\mathrm{subject\ to}\quad \Vert x\Vert_1 \leq \tau

can be solved by SPGL1 using the ``spg_lasso`` function:

``[x,r,g,info] = spg_lasso(A,b,tau,options)``

The ``options`` parameter controls the stopping criteria and behavior of the solver and consists of keyword-value pairs or structures (see the `options <options.html>`_ section for more details). Default options are used when no options are specified.

As an example, let's create random matrix A and vector b, and solve the Lasso problem for :math:`\tau=0.5` using:

::
 
 s = RandStream('mt19937ar','Seed',0);
 A = randn(s,5,10);
 b = randn(s,5,1);
 [x,r,g,info] = spg_lasso(A,b,0.5);


The Lasso solver gives an output similar to

::

 ================================================================================
 SPGL1  v2.0 (28 Nov 2019)
 ================================================================================
 No. rows              :        5      No. columns           :       10
 Initial tau           : 5.00e-01      Two-norm of b         : 1.86e+00
 Optimality tol        : 1.00e-04      Target one-norm of x  : 5.00e-01
 Basis pursuit tol     : 1.00e-04      Maximum iterations    :       50

  Iter      Objective   Relative Gap      gNorm
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

 Products with A     :      14        Total time   (secs) :     0.1
 Products with A'    :      15        Project time (secs) :     0.0
 Newton iterations   :       0        Mat-vec time (secs) :     0.0
 Line search its     :       0
 

The exit status indicates that an optimal solution was found. The output of the solver is given by the ``[x,r,g,info]`` tuple consisting of the solution ``x``, the residual ``r=b-Ax``, the gradient ``g`` of the objective at the solution, and an ``info`` structure that contains information such as the exit status, the number of iterations, and the total run time.

Continuing with our example, we find that if we solve Lasso with :math:`\tau=1`, the required number of iterations exceed the default of 10 times the number of rows in A. We can increase the maximum number of iterations to 100 by calling

::
 
 [x,r,g,info] = spg_lasso(A,b,1,'iterations',100);

In this case around 90 iterations suffice to find an optimal solution.


.. _basis-pursuit denoise:


Basis-pursuit denoise
-----------------------------

The standard basis-pursuit denoise problem is given by

.. math::

   \mathop{\mathrm{minimize}}_{x}\ \Vert x\Vert_1\quad \mathrm{subject\ to}\quad \Vert Ax-b\Vert_2 \leq \sigma


Solving the problem in SPGL1 can be done by calling the ``spg_bpdn`` function

::

 [x,r,g,info] = spg_bpdn(A,b,sigma,options)


In the special case where :math:`\sigma=0`, we can use the basis-pursuit solver:

::
  
  [x,r,g,info] = spg_bp(A,b,options)


The basis-pursuit solution is obtained by solving a sequence of Lasso problems with suitable values of :math:`\tau`. The gradient output ``g`` thus refers to the objective value of the Lasso subproblem corresponding to the latest value of :math:`\tau`. As an example, consider the same example used for the Lasso problem, now solved with basis-pursuit denoise for :math:`\sigma=1.0`:

::
 
 s = RandStream('mt19937ar','Seed',0);
 A = randn(s,5,10);
 b = randn(s,5,1);
 [x,r,g,info] = spg_bpdn(A,b,1.0);


The solver gives an output similar to

::

 ================================================================================
 SPGL1  v2.0 (28 Nov 2019)
 ================================================================================
 No. rows              :        5      No. columns           :       10
 Initial tau           : 0.00e+00      Two-norm of b         : 1.86e+00
 Optimality tol        : 1.00e-04      Target objective      : 1.00e+00
 Basis pursuit tol     : 1.00e-04      Maximum iterations    :       50

  Iter      Objective   Relative Gap  Rel Error      gNorm            tau
     0  1.8616905e+00  0.0000000e+00   4.63e-01  3.363e+00  4.7704625e-01
     1  1.4973610e+00  5.2616060e-01   3.32e-01  2.119e+00
     2  1.2974269e+00  1.9995584e-01   2.29e-01  1.616e+00
     3  1.2030579e+00  8.1971765e-02   1.69e-01  1.771e+00
     4  1.1995095e+00  7.7709077e-02   1.66e-01  1.905e+00  6.2512430e-01
     5  1.0272341e+00  8.8274756e-02   2.65e-02  1.114e+00
     6  1.0212485e+00  8.2144096e-02   2.08e-02  1.094e+00
     7  1.0187911e+00  2.8385712e-02   1.84e-02  9.659e-01
     8  1.0183832e+00  2.2393424e-02   1.81e-02  9.642e-01  6.4454119e-01
     9  9.9870605e-01  1.6574628e-02   1.29e-03  9.452e-01
    10  1.0008367e+00  1.8704785e-02   8.36e-04  1.114e+00
    11  1.0070820e+00  2.4974860e-02   7.03e-03  1.303e+00
    12  9.9794627e-01  1.5816119e-02   2.05e-03  1.010e+00
    13  9.9785653e-01  9.7287655e-03   2.14e-03  9.869e-01  6.4237393e-01
    14  9.9993939e-01  9.6116065e-03   6.06e-05  9.921e-01

 EXIT -- Found a root

 Products with A     :      18        Total time   (secs) :     0.0
 Products with A'    :      16        Project time (secs) :     0.0
 Newton iterations   :       4        Mat-vec time (secs) :     0.0
 Line search its     :       3


The last column of the output shows that four Lasso subproblems with different :math:`\tau` are solved to obtain the basis-pursuit denoise solution.
 


Group norm and MMV
-------------------------

By default, SPGL1 comes with two additional norms. For the group norm we are given subsets of the vector elements, such that each element occurs in exactly one subset. The norm is then defines as the sum of the Euclidean norm of the elements in each set. To solve the group-norm basis-pursuit denoise formulation we can use

::
 
 [x,r,g,info] = spg_group(A,b,groups,sigma,options)


The ``groups`` parameter is a vector that contains the group number for each of the elements in x. The group numbers themselves can be chosen arbitrary, as long as elements in the same group have the same number. Groups do not have to consist of contiguous elements, although in practice they often are. In the following example we use three groups

::

 s = RandStream('mt19937ar','Seed',0);
 A = randn(s,5,10);
 b = randn(s,5,1);
 x = spg_group(A,b,[1,1,1,2,2,2,2,3,3,3], 1.2);

In this case we get a result that is group-sparse: the elements in some groups are all zero, whereas those in other groups are all non-zero:

::

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


A special case of group sparsity is the multiple-measurement vectors (MMV) problem. In this problem we are given a matrix of measurements B, and assume that the unknown matrix X is such that all columns have the same support. This is often achieved by minimizing the sum of Euclidean norms of the rows in X. This can be reformulated as a group-norm problem by appropriately vectorizing B and X, and redefining A. For convenience SPGL1 provides the function

::

 [X,R,G,info] = spg_mmv(A,B,sigma,options);


As an example, we solve

::

 s = RandStream('mt19937ar','Seed',0);
 A = randn(s,5,7);
 B = randn(s,5,6);
 [X,R,G,info] = spg_mmv(A,B,3.5);
 
which gives

::

 X =
 
    0.1102   -0.0296   -0.0177    0.0365   -0.0298   -0.0218
    0.0123    0.0042    0.1348   -0.0577   -0.0208    0.1842
   -0.0020   -0.0001   -0.0008    0.0009   -0.0087    0.0027
         0         0         0         0         0         0
    0.0806    0.0252    0.0703   -0.1165    0.0985    0.0176
         0         0         0         0         0         0
    0.2815    0.1036   -0.1306    0.0197   -0.0596   -0.2179


Generic interface
--------------------

The generic interface to SPGL1 is given by

::

 [x,r,g,info] = spgl1(A,b,tau,sigma,x0,options)
 [x,r,g,info] = spgl1(A,b,tau,sigma,options)
 [x,r,g,info] = spgl1(A,b,tau,options)

The options parameters are optional and can be a mixture of structure objects and key-value pairs. The first option parameter can be a string identifying a predefined parameter set, which can then be modified by the parameters that follow (for examples, see the examples in the `options<options.html>`_ section).

The ``x0`` parameter can be provided to initialize x. If set we need to provide both tau and sigma. To solve the Lasso problem formulation with an initial value for x, set sigma to the empty vector ``[]``. In order to solve basis pursuit denoise we can set tau to 0 or an empty vector ``[]`` (strongly recommended), or provide an initial value for tau (generally not recommended). In case an initial value of tau is specified it is important that this value be smaller than the value for tau at the solution; reduction of tau from a value that is too large can be very time consuming or result in a suboptimal basis-pursuit solution. When the x0 parameter is set to the empty vector ``[]``, a default initial value for x is used.



Custom norms
----------------------

SPGL1 can be extended to solve Lasso and basis-pursuit denoise problems with
custom primal norms. For this we need to provide three functions

1. ``options.primal_norm``, which evaluates the primal norm of a given vector;
2. ``options.dual_norm``, which evaluates the dual norm corresponding to the primal norm;
3. ``options.project``, which evaluates orthogonal projection onto a scaled unit-norm ball corresponding to the primal norm.

Good examples of the definition of these norms can be found in the ``spg_group`` solver.


.. _regularized versions:

Regularization
--------------------

Regularize version of Lasso and basis-pursuit denoise can be obtained by augmenting A with the weighted identity matrix, and augmenting b with a vector of all zero. For convenience, SPGL1 also supports direct regularization, which changes the Lasso formulation to

.. math::

   \mathop{\mathrm{minimize}}_{x}\ {\textstyle\frac{1}{2}}\Vert Ax-b\Vert_2^2 + {\textstyle\frac{\mu}{2}}\Vert x\Vert_2^2\quad\mathrm{subject\ to}\quad \Vert x\Vert_p \leq \tau

and basis-pursuit denoise to:

.. math::

   \mathop{\mathrm{minimize}}_{x}\ \Vert x\Vert_p\quad \mathrm{subject\ to}\quad \left\Vert \left[\begin{array}{c}A\\ \sqrt{\mu}I\end{array}\right]x-\left[\begin{array}{c}b\\0\end{array}\right]\right\Vert_2\leq \sigma

The :math:`\mu` parameter can be specified in the options as ``options.mu`` or as parameter ``'mu'``.



Information structure
----------------------

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

The history fields are updated every iteration if ``options.history`` is set to true.


Exit status
-----------------------

The information structure returned by SPGL1 contains the exit code ``info.stat``, along with a boolean flag that indicates whether the solve was successful or nor (``info.success``) and the error string ``info.statusStr``.

======  ======================  =========  ==========================================
Value   Status code             Success    Status string
======  ======================  =========  ==========================================
1       ``EXIT_ROOT_FOUND``     Yes        Found a root
2       ``EXIT_BPSOL_FOUND``    Yes        Found a BP solution
3       ``EXIT_LEAST_SQUARES``  Yes        Found a least-squares solution
4       ``EXIT_OPTIMAL``        Yes        Optimal solution found
5       ``EXIT_ITERATIONS``     No         Too many iterations
6       ``EXIT_LINE_ERROR``     No         Linesearch error
7       ``EXIT_SUBOPTIMAL_BP``  No         Found a suboptimal BP solution
8       ``EXIT_MATVEC_LIMIT``   No         Maximum matrix-vector operations reached
9       ``EXIT_RUNTIME``        No         Maximum runtime reached
10      ``EXIT_PROJECTION``     No         Inaccurate projection
======  ======================  =========  ==========================================

