Pareto curve
=======================================

The Pareto curve plays a central role in SPGL1. It is given as the lowest residual norm, as a function of :math:`\tau`, the maximum allowable primal norm of :math:`x`.

 .. figure:: fig/ParetoSigma.png
   :scale: 50

   Example Pareto curve (blue) along with a desired misfit level :math:`\sigma` and Lasso parameter :math:`\tau` that attains the desired solution.

Solving the Lasso problem for a given value of :math:`\tau` amounts to evaluating a point on the curve. The basis-pursuit denoise problem is more difficult: we are given :math:`\sigma`  and need to  find the value of :math:`\tau` where the values of the Pareto curve attains the desired value. 


Root finding
---------------------

SPGL1 solves the basis-pursuit problem by solving a sequence of Lasso problems with different values of :math:`\tau`. Given the solution of a Lasso problem, it determines the gradient of the Pareto curve at that point, which allows us to do root-finding. As an illustration, suppose we start at :math:`\tau=0` and solve the Lasso problem for this value (this always gives :math:`x=0`). We then evaluate the gradient, which leaves us at the situation shown in the left plot below. We then update :math:`\tau` to the point where thelinear approximation to the Pareto curve intersect with the desired value of :math:`\sigma`, indicated by the red dot. In this case this gives, say :math:`\tau= 1.2`.

.. image:: fig/RootFind1.png
   :width: 48%
.. image:: fig/RootFind2.png
   :width: 48%

We then solve the Lasso problem with the new value of :math:`\tau`, as shown in the right plot, and repeate the procedure until the solution has a residual norm that is sufficiently close to :math:`\sigma`.


Implementation
------------------------

SPGL1 provides options that control the root-finding process. The most important option, ``options.rootfindMode``, specified whether root finding is done using the primal or the dual mode, with parameter values 0 and 1 respectively.

Primal mode
~~~~~~~~~~~~~~~~~~

When using the primal objective, SPGL1 uses relaxed conditions to determine whether the Lasso subproblem is solved and approximates the gradient based on the primal objective. This root-finding mode can be very fast, but may result in value of :math:`\tau` that exceeds the minimum. As a result, this mode is useful when a quick solution is needed, or for large-scale problems where highly accurate solves of the subproblem are prohibitive.
An exaggerated illustration of root-finding on the Pareto curve in the primal mode is given in the image below. The subproblem is solved approximately for :math:`\tau=1.2`. The solution (indicated by the gray point) may be slighly above the Pareto curve, or the gradient may be too flat. This makes it possible to 'overshoot' the target and obtain a solution with a value :math:`\tau` that exceeds the minimum possible value. 


.. image:: fig/RootFind5.png


Dual mode
~~~~~~~~~~~~~~~~~~

In the dual root-finding mode we apply root finding based on linear under-approximation of the Pareto curve. Provided that the initial :math:`\tau` is below the final value (a value of zero is always highly recommended), there is never any overestimation. As illustrated below, we use the dual objective value (indicated by the gray dot) along with a gradient to determine the next value of :math:`\tau`. A new value for the intersection with the desired misfit level :math:`\sigma` is computed at every iteration and the maximum value is maintained for the next root-finding step. Root finding is initiated when the duality gap is a fraction (``options.rootfindTol``) of the difference between :math:`\sigma` and the primal objective. As the overall gap shrinks towards the solution, increasingly accurate solves are automatically enforced.

.. image:: fig/RootFind6.png

