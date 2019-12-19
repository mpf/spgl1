# Exit status

The information structure returned by SPGL1 contains the exit code
``info.stat``, along with a boolean flag that indicates whether the solve was
successful or nor (``info.success``) and the error string ``info.statusStr``.

```text
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
```
