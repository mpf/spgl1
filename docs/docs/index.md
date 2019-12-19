# SPGL1: A solver for sparse least squares

SPGL1 is an open-source Matlab solver for sparse least-squares. It is designed to solve any one of these three problem formulations:

**Lasso problem**

$$\mathop{\mathrm{minimize}}_{x}\quad {\textstyle\frac{1}{2}}\Vert Ax-b\Vert_2^2\quad\mathrm{subject\ to}\quad \Vert x\Vert_p \leq \tau$$


**Basis pursuit denoise**

$$\mathop{\mathrm{minimize}}_{x}\quad \Vert x\Vert_p\quad \mathrm{subject\ to}\quad \Vert Ax-b\Vert_2 \leq \sigma$$


**Basis pursuit**

$$\mathop{\mathrm{minimize}}_{x}\quad \Vert x\Vert_p\quad \mathrm{subject\ to}\quad Ax=b$$

## Features

* **Linear operator A.** SPGL1 relies on matrix-vector operations `A*x` and `A'*y`, and accepts both explicit matrices (dense or sparse) and functions that evaluate these products.

* **Sparsity regularizer.** Any sparsity-inducing norm $\Vert\cdot\Vert_p$ can be used if the user provides functions evaluating the primal norm, the corresponding dual norm, and a procedure for projecting onto a level-set of the primal norm. The default is the 1-norm.  Several other norms included in SPGL1 are the group (1,2)-norm and the special multiple-measurement vector (MMV) case.

* **Real and complex domains.** SPGL1 is suitable for problems that live in either the real or complex domains. In the complex domain, the correct corresponding 1-norm (sum of magnitudes) is used.

## Feedback

We are glad to hear from you if you find SPGL1 useful, or if you have any suggestions, contributions, or bug reports. Please send these to the SPGL1 authors

* Ewout van den Berg (Email ``<vandenberg.ewout@gmail.com>``)
* [Michael P. Friedlander](https://friedlander.io) (Email: ``<mpf@cs.ubc.ca>``)

## Credits

This research is supported in part by an NSERC Discovery Grant and by a grant for the Office of Naval Research (N00014-17-1-2009).

