# Regularization

Regularize version of Lasso and basis-pursuit denoise can be obtained by
augmenting $A$ with the weighted identity matrix, and augmenting $b$ with a
vector of all zero. For convenience, SPGL1 also supports direct regularization,
which changes the Lasso formulation to

$$\mathop{\mathrm{minimize}}_{x}\quad {\textstyle\frac{1}{2}}\Vert Ax-b\Vert_2^2 + {\textstyle\frac{\mu}{2}}\Vert x\Vert_2^2\quad\mathrm{subject\ to}\quad \Vert x\Vert_p \leq \tau$$

and basis-pursuit denoise to

$$\mathop{\mathrm{minimize}}_{x}\quad \Vert x\Vert_p\quad \mathrm{subject\ to}\quad \left\Vert \left[\begin{array}{c}A\\ \sqrt{\mu}I\end{array}\right]x-\left[\begin{array}{c}b\\0\end{array}\right]\right\Vert_2\leq \sigma.$$

The $\mu$ parameter can be specified in the options as ``options.mu`` or as parameter ``mu``.


