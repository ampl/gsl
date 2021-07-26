.. index::
   single: nonlinear least squares
   single: least squares, nonlinear

*******************************
Nonlinear Least-Squares Fitting
*******************************

.. include:: include.rst

This chapter describes functions for multidimensional nonlinear
least-squares fitting.  There are generally two classes of
algorithms for solving nonlinear least squares problems, which
fall under line search methods and trust region methods.
GSL currently implements only trust region methods and
provides the user with
full access to intermediate steps of the iteration. The user
also has the ability to tune a number of parameters which affect
low-level aspects of the algorithm which can help to accelerate
convergence for the specific problem at hand. GSL provides
two separate interfaces for nonlinear least squares fitting. The
first is designed for small to moderate sized problems, and the
second is designed for very large problems, which may or may not
have significant sparse structure.

The header file :file:`gsl_multifit_nlinear.h` contains prototypes for the
multidimensional nonlinear fitting functions and related declarations
relating to the small to moderate sized systems.

The header file :file:`gsl_multilarge_nlinear.h` contains prototypes for the
multidimensional nonlinear fitting functions and related declarations
relating to large systems.

.. index::
   single: nonlinear least squares, overview

Overview
========

The problem of multidimensional nonlinear least-squares fitting requires
the minimization of the squared residuals of :math:`n` functions,
:math:`f_i`, in :math:`p` parameters, :math:`x_i`,

.. only:: not texinfo

   .. math::

      \Phi(x) &= {1 \over 2} || f(x) ||^2 \\
              &= {1 \over 2} \sum_{i=1}^{n} f_i (x_1, \dots, x_p)^2

.. only:: texinfo

   ::

      \Phi(x) = (1/2) || f(x) ||^2
              = (1/2) \sum_{i=1}^{n} f_i(x_1, ..., x_p)^2 

In trust region methods, the objective (or cost) function :math:`\Phi(x)` is approximated
by a model function :math:`m_k(\delta)` in the vicinity of some point :math:`x_k`. The
model function is often simply a second order Taylor series expansion around the
point :math:`x_k`, ie:

.. only:: not texinfo

   .. math:: \Phi(x_k + \delta) \approx m_k(\delta) = \Phi(x_k) + g_k^T \delta + {1 \over 2} \delta^T B_k \delta

.. only:: texinfo

   ::

      \Phi(x_k + \delta) ~=~ m_k(\delta) = \Phi(x_k) + g_k^T \delta + 1/2 \delta^T B_k \delta

where :math:`g_k = \nabla \Phi(x_k) = J^T f` is the gradient vector at the point :math:`x_k`,
:math:`B_k = \nabla^2 \Phi(x_k)` is the Hessian matrix at :math:`x_k`, or some
approximation to it, and :math:`J` is the :math:`n`-by-:math:`p` Jacobian matrix

.. only:: not texinfo

   .. math:: J_{ij} = \partial f_i / \partial x_j

.. only:: texinfo

   ::

      J_{ij} = d f_i / d x_j

In order to find the next step :math:`\delta`, we minimize the model function
:math:`m_k(\delta)`, but search for solutions only within a region where
we trust that :math:`m_k(\delta)` is a good approximation to the objective
function :math:`\Phi(x_k + \delta)`. In other words,
we seek a solution of the trust region subproblem (TRS)

.. only:: not texinfo

   .. math:: \min_{\delta \in R^p} m_k(\delta) = \Phi(x_k) + g_k^T \delta + {1 \over 2} \delta^T B_k \delta, \qquad\hbox{s.t.}\quad || D_k \delta || \le \Delta_k

.. only:: texinfo

   ::

      \min_(\delta \in R^p) m_k(\delta), s.t. || D_k \delta || <= \Delta_k

where :math:`\Delta_k > 0` is the trust region radius and :math:`D_k` is
a scaling matrix. If :math:`D_k = I`, then the trust region is a ball
of radius :math:`\Delta_k` centered at :math:`x_k`. In some applications,
the parameter vector :math:`x` may have widely different scales. For
example, one parameter might be a temperature on the order of
:math:`10^3` K, while another might be a length on the order of
:math:`10^{-6}` m. In such cases, a spherical trust region may not
be the best choice, since if :math:`\Phi` changes rapidly along
directions with one scale, and more slowly along directions with
a different scale, the model function :math:`m_k` may be a poor
approximation to :math:`\Phi` along the rapidly changing directions.
In such problems, it may be best to use an elliptical trust region,
by setting :math:`D_k` to a diagonal matrix whose entries are designed
so that the scaled step :math:`D_k \delta` has entries of approximately the same
order of magnitude.

The trust region subproblem above normally amounts to solving a
linear least squares system (or multiple systems) for the step
:math:`\delta`. Once :math:`\delta` is computed, it is checked whether
or not it reduces the objective function :math:`\Phi(x)`. A useful
statistic for this is to look at the ratio

.. only:: not texinfo

   .. math:: \rho_k = { \Phi(x_k) - \Phi(x_k + \delta_k) \over m_k(0) - m_k(\delta_k) }

.. only:: texinfo

   ::

      \rho_k = ( \Phi(x_k) - \Phi(x_k + \delta_k) / ( m_k(0) - m_k(\delta_k) )

where the numerator is the actual reduction of the objective function
due to the step :math:`\delta_k`, and the denominator is the predicted
reduction due to the model :math:`m_k`. If :math:`\rho_k` is negative,
it means that the step :math:`\delta_k` increased the objective function
and so it is rejected. If :math:`\rho_k` is positive,
then we have found a step which reduced the objective function and
it is accepted. Furthermore, if :math:`\rho_k` is close to 1,
then this indicates that the model function is a good approximation
to the objective function in the trust region, and so on the next
iteration the trust region is enlarged in order to take more ambitious
steps. When a step is rejected, the trust region is made smaller and
the TRS is solved again. An outline for the general trust region method
used by GSL can now be given.

**Trust Region Algorithm**

1. Initialize: given :math:`x_0`, construct :math:`m_0(\delta)`, :math:`D_0` and :math:`\Delta_0 > 0`

2. For k = 0, 1, 2, ...

   a. If converged, then stop
   b. Solve TRS for trial step :math:`\delta_k`
   c. Evaluate trial step by computing :math:`\rho_k`

      1). if step is accepted, set :math:`x_{k+1} = x_k + \delta_k` and increase radius,
          :math:`\Delta_{k+1} = \alpha \Delta_k`
      2). if step is rejected, set :math:`x_{k+1} = x_k` and decrease radius,
          :math:`\Delta_{k+1} = {\Delta_k \over \beta}`; goto 2(b)

   d. Construct :math:`m_{k+1}(\delta)` and :math:`D_{k+1}`

GSL offers the user a number of different algorithms for solving the trust
region subproblem in 2(b), as well as different choices of scaling matrices
:math:`D_k` and different methods of updating the trust region radius
:math:`\Delta_k`. Therefore, while reasonable default methods are provided,
the user has a lot of control to fine-tune the various steps of the
algorithm for their specific problem.

Solving the Trust Region Subproblem (TRS)
=========================================

Below we describe the methods available for solving the trust region
subproblem. The methods available provide either exact or approximate
solutions to the trust region subproblem. In all algorithms below,
the Hessian matrix :math:`B_k` is approximated as :math:`B_k \approx J_k^T J_k`,
where :math:`J_k = J(x_k)`. In all methods, the solution of the TRS
involves solving a linear least squares system involving the Jacobian
matrix. For small to moderate sized problems (:code:`gsl_multifit_nlinear` interface),
this is accomplished by factoring the full Jacobian matrix, which is provided
by the user, with the Cholesky, QR, or SVD decompositions. For large systems
(:code:`gsl_multilarge_nlinear` interface), the user has two choices. One
is to solve the system iteratively, without needing to store the full
Jacobian matrix in memory. With this method, the user must provide a routine
to calculate the matrix-vector products :math:`J u` or :math:`J^T u` for a given vector :math:`u`.
This iterative method is particularly useful for systems where the Jacobian has
sparse structure, since forming matrix-vector products can be done cheaply. The
second option for large systems involves forming the normal equations matrix
:math:`J^T J` and then factoring it using a Cholesky decomposition. The normal
equations matrix is :math:`p`-by-:math:`p`, typically much smaller than the full
:math:`n`-by-:math:`p` Jacobian, and can usually be stored in memory even if the full
Jacobian matrix cannot. This option is useful for large, dense systems, or if the
iterative method has difficulty converging.

.. index::
   single: Levenberg-Marquardt algorithm
   single: nonlinear least squares, levenberg-marquardt

Levenberg-Marquardt
-------------------

There is a theorem which states that if :math:`\delta_k` is a solution
to the trust region subproblem given above, then there exists
:math:`\mu_k \ge 0` such that

.. only:: not texinfo

   .. math:: \left( B_k + \mu_k D_k^T D_k \right) \delta_k = -g_k

.. only:: texinfo

   ::

      ( B_k + \mu_k D_k^T D_k ) \delta_k = -g_k

with :math:`\mu_k (\Delta_k - ||D_k \delta_k||) = 0`. This
forms the basis of the Levenberg-Marquardt algorithm, which controls
the trust region size by adjusting the parameter :math:`\mu_k`
rather than the radius :math:`\Delta_k` directly. For each radius
:math:`\Delta_k`, there is a unique parameter :math:`\mu_k` which
solves the TRS, and they have an inverse relationship, so that large values of
:math:`\mu_k` correspond to smaller trust regions, while small
values of :math:`\mu_k` correspond to larger trust regions.

With the approximation :math:`B_k \approx J_k^T J_k`, on each iteration,
in order to calculate the step :math:`\delta_k`,
the following linear least squares problem is solved:

.. only:: not texinfo

   .. math::

      \left[
        \begin{array}{c}
          J_k \\
          \sqrt{\mu_k} D_k
        \end{array}
      \right]
      \delta_k =
      -
      \left[
        \begin{array}{c}
          f_k \cr
          0
        \end{array}
      \right]

.. only:: texinfo

   ::

      [J_k; sqrt(mu_k) D_k] \delta_k = - [f_k; 0]

If the step :math:`\delta_k` is accepted, then
:math:`\mu_k` is decreased on the next iteration in order
to take a larger step, otherwise it is increased to take
a smaller step. The Levenberg-Marquardt algorithm provides
an exact solution of the trust region subproblem, but
typically has a higher computational cost per iteration
than the approximate methods discussed below, since it
may need to solve the least squares system above several
times for different values of :math:`\mu_k`.

.. index::
   single: Levenberg-Marquardt algorithm, geodesic acceleration
   single: nonlinear least squares, levenberg-marquardt, geodesic acceleration

Levenberg-Marquardt with Geodesic Acceleration
----------------------------------------------

This method applies a so-called geodesic acceleration correction to
the standard Levenberg-Marquardt step :math:`\delta_k` (Transtrum et al, 2011).
By interpreting :math:`\delta_k` as a first order step along a geodesic in the
model parameter space (ie: a velocity :math:`\delta_k = v_k`), the geodesic
acceleration :math:`a_k` is a second order correction along the
geodesic which is determined by solving the linear least squares system

.. only:: not texinfo

   .. math::

      \left[
        \begin{array}{c}
          J_k \\
          \sqrt{\mu_k} D_k
        \end{array}
      \right]
      a_k =
      -
      \left[
        \begin{array}{c}
          f_{vv}(x_k) \\
          0
        \end{array}
      \right]

.. only:: texinfo

   ::

      [J_k; sqrt(mu_k) D_k] a_k = - [f_vv(x_k); 0]

where :math:`f_{vv}` is the second directional derivative of
the residual vector in the velocity direction :math:`v`,
:math:`f_{vv}(x) = D_v^2 f = \sum_{\alpha\beta} v_{\alpha} v_{\beta} \partial_{\alpha} \partial_{\beta} f(x)`,
where :math:`\alpha` and :math:`\beta` are summed over the :math:`p`
parameters. The new total step is then :math:`\delta_k' = v_k + {1 \over 2}a_k`.
The second order correction :math:`a_k` can be calculated with a modest additional
cost, and has been shown to dramatically reduce the number of iterations
(and expensive Jacobian evaluations) required to reach convergence on a variety
of different problems. In order to utilize the geodesic acceleration, the user must supply a
function which provides the second directional derivative vector
:math:`f_{vv}(x)`, or alternatively the library can use a finite
difference method to estimate this vector with one additional function
evaluation of :math:`f(x + h v)` where :math:`h` is a tunable step size
(see the :code:`h_fvv` parameter description).

.. index::
   single: Dogleg algorithm
   single: nonlinear least squares, dogleg

Dogleg
------

This is Powell's dogleg method, which finds an approximate
solution to the trust region subproblem, by restricting
its search to a piecewise linear "dogleg" path,
composed of the origin, the Cauchy point which represents
the model minimizer along the steepest descent direction,
and the Gauss-Newton point, which is the overall minimizer
of the unconstrained model. The Gauss-Newton step is calculated by
solving

.. math:: J_k \delta_{gn} = -f_k

which is the main computational task for each iteration,
but only needs to be performed once per iteration. If
the Gauss-Newton point is inside the trust region, it is
selected as the step. If it is outside, the method then
calculates the Cauchy point, which is located along the
gradient direction. If the Cauchy point is also outside
the trust region, the method assumes that it is still far
from the minimum and so proceeds along the gradient
direction, truncating the step at the trust region
boundary. If the Cauchy point is inside the trust region,
with the Gauss-Newton point outside, the method
uses a dogleg step, which is a linear combination of the
gradient direction and the Gauss-Newton direction, stopping at the trust
region boundary.

.. index::
   single: double Dogleg algorithm
   single: Dogleg algorithm, double
   single: nonlinear least squares, double dogleg

Double Dogleg
-------------

This method is an improvement over the classical dogleg
algorithm, which attempts to include information about
the Gauss-Newton step while the iteration is still far from
the minimum. When the Cauchy point is inside the trust region
and the Gauss-Newton point is outside, the method computes
a scaled Gauss-Newton point and then takes a dogleg step
between the Cauchy point and the scaled Gauss-Newton point.
The scaling is calculated to ensure that the reduction
in the model :math:`m_k` is about the same as the reduction
provided by the Cauchy point.

Two Dimensional Subspace
------------------------

The dogleg methods restrict the search for the TRS solution
to a 1D curve defined by the Cauchy and Gauss-Newton points.
An improvement to this is to search for a solution using
the full two dimensional subspace spanned by the Cauchy
and Gauss-Newton directions. The dogleg path is of course
inside this subspace, and so this method solves the TRS
at least as accurately as the dogleg methods. Since this
method searches a larger subspace for a solution, it can
converge more quickly than dogleg on some problems. Because
the subspace is only two dimensional, this method is
very efficient and the main computation per iteration is
to determine the Gauss-Newton point.

Steihaug-Toint Conjugate Gradient
---------------------------------

One difficulty of the dogleg methods is calculating the
Gauss-Newton step when the Jacobian matrix is singular. The
Steihaug-Toint method also computes a generalized dogleg
step, but avoids solving for the Gauss-Newton step directly,
instead using an iterative conjugate gradient algorithm. This
method performs well at points where the Jacobian is singular,
and is also suitable for large-scale problems where factoring
the Jacobian matrix could be prohibitively expensive.

Weighted Nonlinear Least-Squares
================================

Weighted nonlinear least-squares fitting minimizes the function

.. only:: not texinfo

   .. math::

      \Phi(x) &= {1 \over 2} || f ||_W^2 \\
              &= {1 \over 2} \sum_{i=1}^{n} w_i f_i (x_1, \dots, x_p)^2

.. only:: texinfo

   ::

      \Phi(x) = (1/2) || f(x) ||_W^2
              = (1/2) \sum_{i=1}^{n} f_i(x_1, ..., x_p)^2 

where :math:`W = \diag(w_1,w_2,...,w_n)` is the weighting matrix,
and :math:`||f||_W^2 = f^T W f`.
The weights :math:`w_i` are commonly defined as :math:`w_i = 1/\sigma_i^2`,
where :math:`\sigma_i` is the error in the :math:`i`-th measurement.
A simple change of variables :math:`\tilde{f} = W^{1 \over 2} f` yields
:math:`\Phi(x) = {1 \over 2} ||\tilde{f}||^2`, which is in the
same form as the unweighted case. The user can either perform this
transform directly on their function residuals and Jacobian, or use
the :func:`gsl_multifit_nlinear_winit` interface which automatically
performs the correct scaling. To manually perform this transformation,
the residuals and Jacobian should be modified according to

.. only:: not texinfo

   .. math::

      \tilde{f}_i & = \sqrt{w_i} f_i = {f_i \over \sigma_i} \\
      \tilde{J}_{ij} & = \sqrt{w_i} { \partial f_i \over \partial x_j } = { 1 \over \sigma_i} { \partial f_i \over \partial x_j }

.. only:: texinfo

   ::

      f~_i = f_i / \sigma_i
      J~_ij = 1 / \sigma_i df_i/dx_j

For large systems, the user must perform their own weighting.

.. _sec_tunable-parameters:

Tunable Parameters
==================

The user can tune nearly all aspects of the iteration at allocation
time. For the :code:`gsl_multifit_nlinear` interface, the user may
modify the :type:`gsl_multifit_nlinear_parameters` structure, which is
defined as follows:

.. type:: gsl_multifit_nlinear_parameters

   ::

      typedef struct
      {
        const gsl_multifit_nlinear_trs *trs;        /* trust region subproblem method */
        const gsl_multifit_nlinear_scale *scale;    /* scaling method */
        const gsl_multifit_nlinear_solver *solver;  /* solver method */
        gsl_multifit_nlinear_fdtype fdtype;         /* finite difference method */
        double factor_up;                           /* factor for increasing trust radius */
        double factor_down;                         /* factor for decreasing trust radius */
        double avmax;                               /* max allowed |a|/|v| */
        double h_df;                                /* step size for finite difference Jacobian */
        double h_fvv;                               /* step size for finite difference fvv */
      } gsl_multifit_nlinear_parameters;

For the :code:`gsl_multilarge_nlinear` interface, the user may
modify the :type:`gsl_multilarge_nlinear_parameters` structure, which is
defined as follows:

.. type:: gsl_multilarge_nlinear_parameters

   ::

      typedef struct
      {
        const gsl_multilarge_nlinear_trs *trs;       /* trust region subproblem method */
        const gsl_multilarge_nlinear_scale *scale;   /* scaling method */
        const gsl_multilarge_nlinear_solver *solver; /* solver method */
        gsl_multilarge_nlinear_fdtype fdtype;        /* finite difference method */
        double factor_up;                            /* factor for increasing trust radius */
        double factor_down;                          /* factor for decreasing trust radius */
        double avmax;                                /* max allowed |a|/|v| */
        double h_df;                                 /* step size for finite difference Jacobian */
        double h_fvv;                                /* step size for finite difference fvv */
        size_t max_iter;                             /* maximum iterations for trs method */
        double tol;                                  /* tolerance for solving trs */
      } gsl_multilarge_nlinear_parameters;

Each of these parameters is discussed in further detail below.

.. type:: gsl_multifit_nlinear_trs
          gsl_multilarge_nlinear_trs

   The parameter :data:`trs` determines the method used to solve the trust region
   subproblem, and may be selected from the following choices,

   .. var:: gsl_multifit_nlinear_trs * gsl_multifit_nlinear_trs_lm
            gsl_multilarge_nlinear_trs * gsl_multilarge_nlinear_trs_lm

      This selects the Levenberg-Marquardt algorithm.

   .. var:: gsl_multifit_nlinear_trs * gsl_multifit_nlinear_trs_lmaccel
            gsl_multilarge_nlinear_trs * gsl_multilarge_nlinear_trs_lmaccel

      This selects the Levenberg-Marquardt algorithm with geodesic
      acceleration.

   .. var:: gsl_multifit_nlinear_trs * gsl_multifit_nlinear_trs_dogleg
            gsl_multilarge_nlinear_trs * gsl_multilarge_nlinear_trs_dogleg

      This selects the dogleg algorithm.

   .. var:: gsl_multifit_nlinear_trs * gsl_multifit_nlinear_trs_ddogleg
            gsl_multilarge_nlinear_trs * gsl_multilarge_nlinear_trs_ddogleg

      This selects the double dogleg algorithm.

   .. var:: gsl_multifit_nlinear_trs * gsl_multifit_nlinear_trs_subspace2D
            gsl_multilarge_nlinear_trs * gsl_multilarge_nlinear_trs_subspace2D

      This selects the 2D subspace algorithm.

   .. var:: gsl_multilarge_nlinear_trs * gsl_multilarge_nlinear_trs_cgst

      This selects the Steihaug-Toint conjugate gradient algorithm. This
      method is available only for large systems.

.. type:: gsl_multifit_nlinear_scale
          gsl_multilarge_nlinear_scale

   The parameter :data:`scale` determines the diagonal scaling matrix :math:`D` and
   may be selected from the following choices,

   .. var:: gsl_multifit_nlinear_scale * gsl_multifit_nlinear_scale_more
            gsl_multilarge_nlinear_scale * gsl_multilarge_nlinear_scale_more

      This damping strategy was suggested by |More|, and
      corresponds to :math:`D^T D = \max(\diag(J^T J))`,
      in other words the maximum elements of
      :math:`\diag(J^T J)` encountered thus far in the iteration.
      This choice of :math:`D` makes the problem scale-invariant,
      so that if the model parameters :math:`x_i` are each scaled
      by an arbitrary constant, :math:`\tilde{x}_i = a_i x_i`, then
      the sequence of iterates produced by the algorithm would
      be unchanged. This method can work very well in cases
      where the model parameters have widely different scales
      (ie: if some parameters are measured in nanometers, while others
      are measured in degrees Kelvin). This strategy has been proven
      effective on a large class of problems and so it is the library
      default, but it may not be the best choice for all problems.

   .. var:: gsl_multifit_nlinear_scale * gsl_multifit_nlinear_scale_levenberg
            gsl_multilarge_nlinear_scale * gsl_multilarge_nlinear_scale_levenberg

      This damping strategy was originally suggested by Levenberg, and
      corresponds to :math:`D^T D = I`. This method has also proven
      effective on a large class of problems, but is not scale-invariant.
      However, some authors (e.g. Transtrum and Sethna 2012) argue
      that this choice is better for problems which are susceptible
      to parameter evaporation (ie: parameters go to infinity)

   .. var:: gsl_multifit_nlinear_scale * gsl_multifit_nlinear_scale_marquardt
            gsl_multilarge_nlinear_scale * gsl_multilarge_nlinear_scale_marquardt

      This damping strategy was suggested by Marquardt, and
      corresponds to :math:`D^T D = \diag(J^T J)`. This
      method is scale-invariant, but it is generally considered
      inferior to both the Levenberg and |More| strategies, though
      may work well on certain classes of problems.

.. type:: gsl_multifit_nlinear_solver
          gsl_multilarge_nlinear_solver

   Solving the trust region subproblem on each iteration almost always
   requires the solution of the following linear least squares system

   .. only:: not texinfo

      .. math::

         \left[
           \begin{array}{c}
             J \\
             \sqrt{\mu} D
           \end{array}
         \right]
         \delta =
         -
         \left[
           \begin{array}{c}
             f \\
             0
           \end{array}
         \right]

   .. only:: texinfo

      ::

         [J; sqrt(mu) D] \delta = - [f; 0]

   The :data:`solver` parameter determines how the system is
   solved and can be selected from the following choices:

   .. var:: gsl_multifit_nlinear_solver * gsl_multifit_nlinear_solver_qr

      This method solves the system using a rank revealing QR
      decomposition of the Jacobian :math:`J`. This method will
      produce reliable solutions in cases where the Jacobian
      is rank deficient or near-singular but does require about
      twice as many operations as the Cholesky method discussed
      below.

   .. var:: gsl_multifit_nlinear_solver * gsl_multifit_nlinear_solver_cholesky
            gsl_multilarge_nlinear_solver * gsl_multilarge_nlinear_solver_cholesky

      This method solves the alternate normal equations problem

      .. only:: not texinfo

         .. math:: \left( J^T J + \mu D^T D \right) \delta = -J^T f

      .. only:: texinfo

         ::

            ( J^T J + \mu D^T D ) \delta = -J^T f

      by using a Cholesky decomposition of the matrix
      :math:`J^T J + \mu D^T D`. This method is faster than the
      QR approach, however it is susceptible to numerical instabilities
      if the Jacobian matrix is rank deficient or near-singular. In
      these cases, an attempt is made to reduce the condition number
      of the matrix using Jacobi preconditioning, but for highly
      ill-conditioned problems the QR approach is better. If it is
      known that the Jacobian matrix is well conditioned, this method
      is accurate and will perform faster than the QR approach.

   .. var:: gsl_multifit_nlinear_solver * gsl_multifit_nlinear_solver_mcholesky
            gsl_multilarge_nlinear_solver * gsl_multilarge_nlinear_solver_mcholesky

      This method solves the alternate normal equations problem

      .. only:: not texinfo

         .. math:: \left( J^T J + \mu D^T D \right) \delta = -J^T f

      .. only:: texinfo

         ::

            ( J^T J + \mu D^T D ) \delta = -J^T f

      by using a modified Cholesky decomposition of the matrix
      :math:`J^T J + \mu D^T D`. This is more suitable for the dogleg
      methods where the parameter :math:`\mu = 0`, and the matrix
      :math:`J^T J` may be ill-conditioned or indefinite causing the standard Cholesky
      decomposition to fail. This method is based on Level 2 BLAS
      and is thus slower than the standard Cholesky decomposition, which
      is based on Level 3 BLAS.

   .. var:: gsl_multifit_nlinear_solver * gsl_multifit_nlinear_solver_svd

      This method solves the system using a singular value
      decomposition of the Jacobian :math:`J`. This method will
      produce the most reliable solutions for ill-conditioned Jacobians
      but is also the slowest solver method.

.. type:: gsl_multifit_nlinear_fdtype

   The parameter :data:`fdtype` specifies whether to use forward or centered
   differences when approximating the Jacobian. This is only
   used when an analytic Jacobian is not provided to the solver.
   This parameter may be set to one of the following choices.

   .. macro:: GSL_MULTIFIT_NLINEAR_FWDIFF

      This specifies a forward finite difference to approximate
      the Jacobian matrix. The Jacobian matrix will be calculated as

      .. only:: not texinfo

         .. math:: J_{ij} = {1 \over \Delta_j} \left( f_i(x + \Delta_j e_j) - f_i(x) \right)

      .. only:: texinfo

         ::

            J_ij = 1 / \Delta_j ( f_i(x + \Delta_j e_j) - f_i(x) )

      where :math:`\Delta_j = h |x_j|` and :math:`e_j` is the standard
      :math:`j`-th Cartesian unit basis vector so that
      :math:`x + \Delta_j e_j` represents a small (forward) perturbation of
      the :math:`j`-th parameter by an amount :math:`\Delta_j`. The perturbation
      :math:`\Delta_j` is proportional to the current value :math:`|x_j|` which
      helps to calculate an accurate Jacobian when the various parameters have
      different scale sizes. The value of :math:`h` is specified by the :code:`h_df`
      parameter. The accuracy of this method is :math:`O(h)`, and evaluating this
      matrix requires an additional :math:`p` function evaluations.

   .. macro:: GSL_MULTIFIT_NLINEAR_CTRDIFF

      This specifies a centered finite difference to approximate
      the Jacobian matrix. The Jacobian matrix will be calculated as

      .. only:: not texinfo

         .. math:: J_{ij} = {1 \over \Delta_j} \left( f_i(x + {1 \over 2} \Delta_j e_j) - f_i(x - {1 \over 2} \Delta_j e_j) \right)

      .. only:: texinfo

         ::

            J_ij = 1 / \Delta_j ( f_i(x + 1/2 \Delta_j e_j) - f_i(x - 1/2 \Delta_j e_j) )

      See above for a description of :math:`\Delta_j`. The accuracy of this
      method is :math:`O(h^2)`, but evaluating this
      matrix requires an additional :math:`2p` function evaluations.

:code:`double factor_up`

When a step is accepted, the trust region radius will be increased
by this factor. The default value is :math:`3`.

:code:`double factor_down`

When a step is rejected, the trust region radius will be decreased
by this factor. The default value is :math:`2`.

:code:`double avmax`

When using geodesic acceleration to solve a nonlinear least squares problem,
an important parameter to monitor is the ratio of the acceleration term
to the velocity term,

.. only:: not texinfo

   .. math:: { ||a|| \over ||v|| }

.. only:: texinfo

   ::

      |a| / |v|

If this ratio is small, it means the acceleration correction
is contributing very little to the step. This could be because
the problem is not "nonlinear" enough to benefit from
the acceleration. If the ratio is large (:math:`> 1`) it
means that the acceleration is larger than the velocity,
which shouldn't happen since the step represents a truncated
series and so the second order term :math:`a` should be smaller than
the first order term :math:`v` to guarantee convergence.
Therefore any steps with a ratio larger than the parameter
:data:`avmax` are rejected. :data:`avmax` is set to 0.75 by default.
For problems which experience difficulty converging, this threshold
could be lowered.

:code:`double h_df`

This parameter specifies the step size for approximating the
Jacobian matrix with finite differences. It is set to
:math:`\sqrt{\epsilon}` by default, where :math:`\epsilon`
is :macro:`GSL_DBL_EPSILON`.

:code:`double h_fvv`

When using geodesic acceleration, the user must either supply
a function to calculate :math:`f_{vv}(x)` or the library
can estimate this second directional derivative using a finite
difference method. When using finite differences, the library
must calculate :math:`f(x + h v)` where :math:`h` represents
a small step in the velocity direction. The parameter
:data:`h_fvv` defines this step size and is set to 0.02 by
default.

Initializing the Solver
=======================

.. type:: gsl_multifit_nlinear_type

   This structure specifies the type of algorithm which will be used
   to solve a nonlinear least squares problem. It may be selected from the
   following choices,

   .. var:: gsl_multifit_nlinear_type * gsl_multifit_nlinear_trust

      This specifies a trust region method. It is currently the only implemented
      nonlinear least squares method.

.. function:: gsl_multifit_nlinear_workspace * gsl_multifit_nlinear_alloc (const gsl_multifit_nlinear_type * T, const gsl_multifit_nlinear_parameters * params, const size_t n, const size_t p)
              gsl_multilarge_nlinear_workspace * gsl_multilarge_nlinear_alloc (const gsl_multilarge_nlinear_type * T, const gsl_multilarge_nlinear_parameters * params, const size_t n, const size_t p)

   These functions return a pointer to a newly allocated instance of a
   derivative solver of type :data:`T` for :data:`n` observations and :data:`p`
   parameters. The :data:`params` input specifies a tunable set of
   parameters which will affect important details in each iteration
   of the trust region subproblem algorithm. It is recommended to start
   with the suggested default parameters (see
   :func:`gsl_multifit_nlinear_default_parameters` and
   :func:`gsl_multilarge_nlinear_default_parameters`) and then tune
   the parameters once the code is working correctly. See
   :ref:`sec_tunable-parameters`.
   for descriptions of the various parameters.
   For example, the following code creates an instance of a
   Levenberg-Marquardt solver for 100 data points and 3 parameters,
   using suggested defaults::

      const gsl_multifit_nlinear_type * T = gsl_multifit_nlinear_trust;
      gsl_multifit_nlinear_parameters params = gsl_multifit_nlinear_default_parameters();
      gsl_multifit_nlinear_workspace * w = gsl_multifit_nlinear_alloc (T, &params, 100, 3);

   The number of observations :data:`n` must be greater than or equal to
   parameters :data:`p`.

   If there is insufficient memory to create the solver then the function
   returns a null pointer and the error handler is invoked with an error
   code of :macro:`GSL_ENOMEM`.

.. function:: gsl_multifit_nlinear_parameters gsl_multifit_nlinear_default_parameters (void)
              gsl_multilarge_nlinear_parameters gsl_multilarge_nlinear_default_parameters (void)

   These functions return a set of recommended default parameters
   for use in solving nonlinear least squares problems. The user
   can tune each parameter to improve the performance on their
   particular problem, see :ref:`sec_tunable-parameters`.

.. function:: int gsl_multifit_nlinear_init (const gsl_vector * x, gsl_multifit_nlinear_fdf * fdf, gsl_multifit_nlinear_workspace * w)
              int gsl_multifit_nlinear_winit (const gsl_vector * x, const gsl_vector * wts, gsl_multifit_nlinear_fdf * fdf, gsl_multifit_nlinear_workspace * w)
              int gsl_multilarge_nlinear_init (const gsl_vector * x, gsl_multilarge_nlinear_fdf * fdf, gsl_multilarge_nlinear_workspace * w)

   These functions initialize, or reinitialize, an existing workspace :data:`w`
   to use the system :data:`fdf` and the initial guess
   :data:`x`. See :ref:`sec_providing-function-minimized`
   for a description of the :data:`fdf` structure.

   Optionally, a weight vector :data:`wts` can be given to perform
   a weighted nonlinear regression. Here, the weighting matrix is
   :math:`W = \diag(w_1,w_2,...,w_n)`.

.. function:: void gsl_multifit_nlinear_free (gsl_multifit_nlinear_workspace * w)
              void gsl_multilarge_nlinear_free (gsl_multilarge_nlinear_workspace * w)

   These functions free all the memory associated with the workspace :data:`w`.

.. function:: const char * gsl_multifit_nlinear_name (const gsl_multifit_nlinear_workspace * w)
              const char * gsl_multilarge_nlinear_name (const gsl_multilarge_nlinear_workspace * w)

   These functions return a pointer to the name of the solver.  For example::

      printf ("w is a '%s' solver\n", gsl_multifit_nlinear_name (w));

   would print something like :code:`w is a 'trust-region' solver`.

.. function:: const char * gsl_multifit_nlinear_trs_name (const gsl_multifit_nlinear_workspace * w)
              const char * gsl_multilarge_nlinear_trs_name (const gsl_multilarge_nlinear_workspace * w)

   These functions return a pointer to the name of the trust region subproblem
   method.  For example::

      printf ("w is a '%s' solver\n", gsl_multifit_nlinear_trs_name (w));

   would print something like :code:`w is a 'levenberg-marquardt' solver`.

.. _sec_providing-function-minimized:

Providing the Function to be Minimized
======================================

The user must provide :math:`n` functions of :math:`p` variables for the
minimization algorithm to operate on.  In order to allow for
arbitrary parameters the functions are defined by the following data
types:

.. type:: gsl_multifit_nlinear_fdf

   This data type defines a general system of functions with arbitrary parameters,
   the corresponding Jacobian matrix of derivatives, and optionally the
   second directional derivative of the functions for geodesic acceleration.

   :code:`int (* f) (const gsl_vector * x, void * params, gsl_vector * f)`

      This function should store the :math:`n` components of the vector
      :math:`f(x)` in :data:`f` for argument :data:`x` and arbitrary parameters :data:`params`,
      returning an appropriate error code if the function cannot be computed.

   :code:`int (* df) (const gsl_vector * x, void * params, gsl_matrix * J)`

      This function should store the :data:`n`-by-:data:`p` matrix result

      .. only:: not texinfo

         .. math:: J_{ij} = \partial f_i(x) / \partial x_j

      .. only:: texinfo

         ::

            J_ij = d f_i(x) / d x_j

      in :data:`J` for argument :data:`x` 
      and arbitrary parameters :data:`params`, returning an appropriate error code if the
      matrix cannot be computed. If an analytic Jacobian is unavailable, or too expensive
      to compute, this function pointer may be set to :code:`NULL`, in which
      case the Jacobian will be internally computed using finite difference approximations
      of the function :data:`f`.

   :code:`int (* fvv) (const gsl_vector * x, const gsl_vector * v, void * params, gsl_vector * fvv)`

      When geodesic acceleration is enabled, this function should store the
      :math:`n` components of the vector
      :math:`f_{vv}(x) = \sum_{\alpha\beta} v_{\alpha} v_{\beta} {\partial \over \partial x_{\alpha}} {\partial \over \partial x_{\beta}} f(x)`,
      representing second directional derivatives of the function to be minimized,
      into the output :data:`fvv`. The parameter vector is provided in :data:`x` and
      the velocity vector is provided in :data:`v`, both of which have :math:`p`
      components. The arbitrary parameters are given in :data:`params`. If
      analytic expressions for :math:`f_{vv}(x)` are unavailable or too difficult
      to compute, this function pointer may be set to :code:`NULL`, in which case
      :math:`f_{vv}(x)` will be computed internally using a finite difference
      approximation.

   :code:`size_t n`

       the number of functions, i.e. the number of components of the
       vector :data:`f`.

   :code:`size_t p`

      the number of independent variables, i.e. the number of components of
      the vector :data:`x`.

   :code:`void * params`

      a pointer to the arbitrary parameters of the function.

   :code:`size_t nevalf`

      This does not need to be set by the user. It counts the number of
      function evaluations and is initialized by the :code:`_init` function.

   :code:`size_t nevaldf`

      This does not need to be set by the user. It counts the number of
      Jacobian evaluations and is initialized by the :code:`_init` function.

   :code:`size_t nevalfvv`

      This does not need to be set by the user. It counts the number of
      :math:`f_{vv}(x)` evaluations and is initialized by the :code:`_init` function.

.. type:: gsl_multilarge_nlinear_fdf

   This data type defines a general system of functions with arbitrary parameters,
   a function to compute :math:`J u` or :math:`J^T u` for a given vector :math:`u`,
   the normal equations matrix :math:`J^T J`,
   and optionally the second directional derivative of the functions for geodesic acceleration.

   :code:`int (* f) (const gsl_vector * x, void * params, gsl_vector * f)`

      This function should store the :math:`n` components of the vector
      :math:`f(x)` in :data:`f` for argument :data:`x` and arbitrary parameters :data:`params`,
      returning an appropriate error code if the function cannot be computed.

   :code:`int (* df) (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x, const gsl_vector * u, void * params, gsl_vector * v, gsl_matrix * JTJ)`

      If :data:`TransJ` is equal to :code:`CblasNoTrans`, then this function should
      compute the matrix-vector product :math:`J u` and store the result in :data:`v`.
      If :data:`TransJ` is equal to :code:`CblasTrans`, then this function should
      compute the matrix-vector product :math:`J^T u` and store the result in :data:`v`.
      Additionally, the normal equations matrix :math:`J^T J` should be stored in the
      lower half of :data:`JTJ`. The input matrix :data:`JTJ` could be set to :code:`NULL`,
      for example by iterative methods which do not require this matrix, so the user
      should check for this prior to constructing the matrix.
      The input :data:`params` contains the arbitrary parameters.

   :code:`int (* fvv) (const gsl_vector * x, const gsl_vector * v, void * params, gsl_vector * fvv)`

      When geodesic acceleration is enabled, this function should store the
      :math:`n` components of the vector
      :math:`f_{vv}(x) = \sum_{\alpha\beta} v_{\alpha} v_{\beta} {\partial \over \partial x_{\alpha}} {\partial \over \partial x_{\beta}} f(x)`,
      representing second directional derivatives of the function to be minimized,
      into the output :data:`fvv`. The parameter vector is provided in :data:`x` and
      the velocity vector is provided in :data:`v`, both of which have :math:`p`
      components. The arbitrary parameters are given in :data:`params`. If
      analytic expressions for :math:`f_{vv}(x)` are unavailable or too difficult
      to compute, this function pointer may be set to :code:`NULL`, in which case
      :math:`f_{vv}(x)` will be computed internally using a finite difference
      approximation.

   :code:`size_t n`

      the number of functions, i.e. the number of components of the
      vector :data:`f`.

   :code:`size_t p`

      the number of independent variables, i.e. the number of components of
      the vector :data:`x`.

   :code:`void * params`

      a pointer to the arbitrary parameters of the function.

   :code:`size_t nevalf`

      This does not need to be set by the user. It counts the number of
      function evaluations and is initialized by the :code:`_init` function.

   :code:`size_t nevaldfu`

      This does not need to be set by the user. It counts the number of
      Jacobian matrix-vector evaluations (:math:`J u` or :math:`J^T u`) and
      is initialized by the :code:`_init` function.

   :code:`size_t nevaldf2`

      This does not need to be set by the user. It counts the number of
      :math:`J^T J` evaluations and is initialized by the :code:`_init` function.

   :code:`size_t nevalfvv`

      This does not need to be set by the user. It counts the number of
      :math:`f_{vv}(x)` evaluations and is initialized by the :code:`_init` function.

Note that when fitting a non-linear model against experimental data,
the data is passed to the functions above using the
:data:`params` argument and the trial best-fit parameters through the
:data:`x` argument.

Iteration
=========

The following functions drive the iteration of each algorithm.  Each
function performs one iteration of the trust region method and updates
the state of the solver.

.. function:: int gsl_multifit_nlinear_iterate (gsl_multifit_nlinear_workspace * w)
              int gsl_multilarge_nlinear_iterate (gsl_multilarge_nlinear_workspace * w)

   These functions perform a single iteration of the solver :data:`w`.  If
   the iteration encounters an unexpected problem then an error code will
   be returned.  The solver workspace maintains a current estimate of the
   best-fit parameters at all times.

The solver workspace :data:`w` contains the following entries, which can
be used to track the progress of the solution:

:code:`gsl_vector * x`

  The current position, length :math:`p`.

:code:`gsl_vector * f`

  The function residual vector at the current position :math:`f(x)`, length
  :math:`n`.

:code:`gsl_matrix * J`

  The Jacobian matrix at the current position :math:`J(x)`, size
  :math:`n`-by-:math:`p` (only for :code:`gsl_multifit_nlinear` interface).

:code:`gsl_vector * dx`

  The difference between the current position and the previous position,
  i.e. the last step :math:`\delta`, taken as a vector, length :math:`p`.

These quantities can be accessed with the following functions,

.. function:: gsl_vector * gsl_multifit_nlinear_position (const gsl_multifit_nlinear_workspace * w)
              gsl_vector * gsl_multilarge_nlinear_position (const gsl_multilarge_nlinear_workspace * w)

   These functions return the current position :math:`x` (i.e. best-fit
   parameters) of the solver :data:`w`.

.. function:: gsl_vector * gsl_multifit_nlinear_residual (const gsl_multifit_nlinear_workspace * w)
              gsl_vector * gsl_multilarge_nlinear_residual (const gsl_multilarge_nlinear_workspace * w)

   These functions return the current residual vector :math:`f(x)` of the
   solver :data:`w`.  For weighted systems, the residual vector includes the
   weighting factor :math:`\sqrt{W}`.

.. function:: gsl_matrix * gsl_multifit_nlinear_jac (const gsl_multifit_nlinear_workspace * w)

   This function returns a pointer to the :math:`n`-by-:math:`p` Jacobian matrix for the
   current iteration of the solver :data:`w`. This function is available only for the
   :code:`gsl_multifit_nlinear` interface.

.. function:: size_t gsl_multifit_nlinear_niter (const gsl_multifit_nlinear_workspace * w)
              size_t gsl_multilarge_nlinear_niter (const gsl_multilarge_nlinear_workspace * w)

   These functions return the number of iterations performed so far.
   The iteration counter is updated on each call to the
   :code:`_iterate` functions above, and reset to 0 in the
   :code:`_init` functions.

.. function:: int gsl_multifit_nlinear_rcond (double * rcond, const gsl_multifit_nlinear_workspace * w)
              int gsl_multilarge_nlinear_rcond (double * rcond, const gsl_multilarge_nlinear_workspace * w)

   This function estimates the reciprocal condition number
   of the Jacobian matrix at the current position :math:`x` and
   stores it in :data:`rcond`. The computed value is only an estimate
   to give the user a guideline as to the conditioning of their particular
   problem. Its calculation is based on which factorization
   method is used (Cholesky, QR, or SVD). 

   * For the Cholesky solver, the matrix :math:`J^T J` is factored at each
     iteration. Therefore this function will estimate the 1-norm condition number
     :math:`rcond^2 = 1/(||J^T J||_1 \cdot ||(J^T J)^{-1}||_1)`

   * For the QR solver, :math:`J` is factored as :math:`J = Q R` at each
     iteration. For simplicity, this function calculates the 1-norm conditioning of
     only the :math:`R` factor, :math:`rcond = 1 / (||R||_1 \cdot ||R^{-1}||_1)`.
     This can be computed efficiently since :math:`R` is upper triangular.

   * For the SVD solver, in order to efficiently solve the trust region
     subproblem, the matrix which is factored is :math:`J D^{-1}`, instead of
     :math:`J` itself. The resulting singular values are used to provide
     the 2-norm reciprocal condition number, as :math:`rcond = \sigma_{min} / \sigma_{max}`.
     Note that when using |More| scaling, :math:`D \ne I` and the resulting
     :data:`rcond` estimate may be significantly different from the true
     :data:`rcond` of :math:`J` itself.

.. function:: double gsl_multifit_nlinear_avratio (const gsl_multifit_nlinear_workspace * w)
              double gsl_multilarge_nlinear_avratio (const gsl_multilarge_nlinear_workspace * w)

   This function returns the current ratio :math:`|a| / |v|` of the acceleration correction term to
   the velocity step term. The acceleration term is computed only by the
   :type:`gsl_multifit_nlinear_trs_lmaccel` and :type:`gsl_multilarge_nlinear_trs_lmaccel` methods, so
   this ratio will be zero for other TRS methods.

.. index::
   single: nonlinear fitting, stopping parameters, convergence

Testing for Convergence
=======================

A minimization procedure should stop when one of the following conditions is
true:

* A minimum has been found to within the user-specified precision.
* A user-specified maximum number of iterations has been reached.
* An error has occurred.

The handling of these conditions is under user control.  The functions
below allow the user to test the current estimate of the best-fit
parameters in several standard ways.

.. function:: int gsl_multifit_nlinear_test (const double xtol, const double gtol, const double ftol, int * info, const gsl_multifit_nlinear_workspace * w)
              int gsl_multilarge_nlinear_test (const double xtol, const double gtol, const double ftol, int * info, const gsl_multilarge_nlinear_workspace * w)

   These functions test for convergence of the minimization method
   using the following criteria:

   * Testing for a small step size relative to the current parameter vector

     .. math:: |\delta_i| \le xtol (|x_i| + xtol)

     for each :math:`0 <= i < p`. Each element of the step vector :math:`\delta`
     is tested individually in case the different parameters have widely
     different scales. Adding :data:`xtol` to :math:`|x_i|` helps the test avoid
     breaking down in situations where the true solution value :math:`x_i = 0`.
     If this test succeeds, :data:`info` is set to 1 and the function
     returns :macro:`GSL_SUCCESS`.

     A general guideline for selecting the step tolerance is to choose
     :math:`xtol = 10^{-d}` where :math:`d` is the number of accurate
     decimal digits desired in the solution :math:`x`. See Dennis and
     Schnabel for more information.

   * Testing for a small gradient (:math:`g = \nabla \Phi(x) = J^T f`)
     indicating a local function minimum:

     .. only:: not texinfo

        .. math:: \max_i |g_i \times \max(x_i, 1)| \le gtol \times \max(\Phi(x), 1)

     .. only:: texinfo

        ::

           ||g||_inf <= gtol

     This expression tests whether the ratio
     :math:`(\nabla \Phi)_i x_i / \Phi` is small. Testing this scaled gradient
     is a better than :math:`\nabla \Phi` alone since it is a dimensionless
     quantity and so independent of the scale of the problem. The
     :code:`max` arguments help ensure the test doesn't break down in
     regions where :math:`x_i` or :math:`\Phi(x)` are close to 0.
     If this test succeeds, :data:`info` is set to 2 and the function
     returns :macro:`GSL_SUCCESS`.

     A general guideline for choosing the gradient tolerance is to set
     :code:`gtol = GSL_DBL_EPSILON^(1/3)`. See Dennis and Schnabel for
     more information.

   If none of the tests succeed, :data:`info` is set to 0 and the
   function returns :macro:`GSL_CONTINUE`, indicating further iterations
   are required.

High Level Driver
=================

These routines provide a high level wrapper that combines the iteration
and convergence testing for easy use.

.. function:: int gsl_multifit_nlinear_driver (const size_t maxiter, const double xtol, const double gtol, const double ftol, void (* callback)(const size_t iter, void * params, const gsl_multifit_linear_workspace * w), void * callback_params, int * info, gsl_multifit_nlinear_workspace * w)
              int gsl_multilarge_nlinear_driver (const size_t maxiter, const double xtol, const double gtol, const double ftol, void (* callback)(const size_t iter, void * params, const gsl_multilarge_linear_workspace * w), void * callback_params, int * info, gsl_multilarge_nlinear_workspace * w)

   These functions iterate the nonlinear least squares solver :data:`w` for a
   maximum of :data:`maxiter` iterations. After each iteration, the system is
   tested for convergence with the error tolerances :data:`xtol`, :data:`gtol` and :data:`ftol`.
   Additionally, the user may supply a callback function :data:`callback`
   which is called after each iteration, so that the user may save or print
   relevant quantities for each iteration. The parameter :data:`callback_params`
   is passed to the :data:`callback` function. The parameters :data:`callback`
   and :data:`callback_params` may be set to :code:`NULL` to disable this feature.
   Upon successful convergence, the function returns :macro:`GSL_SUCCESS`
   and sets :data:`info` to the reason for convergence (see
   :func:`gsl_multifit_nlinear_test`). If the function has not
   converged after :data:`maxiter` iterations, :macro:`GSL_EMAXITER` is
   returned. In rare cases, during an iteration the algorithm may
   be unable to find a new acceptable step :math:`\delta` to take. In
   this case, :macro:`GSL_ENOPROG` is returned indicating no further
   progress can be made. If your problem is having difficulty converging,
   see :ref:`sec_nlinear-troubleshooting` for further guidance.

.. index::
   single: best-fit parameters, covariance
   single: least squares, covariance of best-fit parameters
   single: covariance matrix, nonlinear fits

Covariance matrix of best fit parameters
========================================

.. function:: int gsl_multifit_nlinear_covar (const gsl_matrix * J, const double epsrel, gsl_matrix * covar)
              int gsl_multilarge_nlinear_covar (gsl_matrix * covar, gsl_multilarge_nlinear_workspace * w)

   This function computes the covariance matrix of best-fit parameters
   using the Jacobian matrix :data:`J` and stores it in :data:`covar`.
   The parameter :data:`epsrel` is used to remove linear-dependent columns
   when :data:`J` is rank deficient.

   The covariance matrix is given by,

   .. math:: C = (J^T J)^{-1}

   or in the weighted case,

   .. math:: C = (J^T W J)^{-1}

   and is computed using the factored form of the Jacobian (Cholesky, QR, or SVD).
   Any columns of :math:`R` which satisfy 

   .. math:: |R_{kk}| \leq epsrel |R_{11}|

   are considered linearly-dependent and are excluded from the covariance
   matrix (the corresponding rows and columns of the covariance matrix are
   set to zero).

   If the minimisation uses the weighted least-squares function
   :math:`f_i = (Y(x, t_i) - y_i) / \sigma_i` then the covariance
   matrix above gives the statistical error on the best-fit parameters
   resulting from the Gaussian errors :math:`\sigma_i` on 
   the underlying data :math:`y_i`.  This can be verified from the relation 
   :math:`\delta f = J \delta c` and the fact that the fluctuations in :math:`f`
   from the data :math:`y_i` are normalised by :math:`\sigma_i` and 
   so satisfy
   
   .. only:: not texinfo
   
      .. math:: \langle \delta f \delta f^T \rangle = I

   .. only:: texinfo

      ::

         <\delta f \delta f^T> = I

   For an unweighted least-squares function :math:`f_i = (Y(x, t_i) - y_i)`
   the covariance matrix above should be multiplied by the variance
   of the residuals about the best-fit :math:`\sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p)`
   to give the variance-covariance
   matrix :math:`\sigma^2 C`.  This estimates the statistical error on the
   best-fit parameters from the scatter of the underlying data.

   For more information about covariance matrices see
   :ref:`Linear Least-Squares Overview <sec_lls-overview>`.

.. _sec_nlinear-troubleshooting:

Troubleshooting
===============

When developing a code to solve a nonlinear least squares problem,
here are a few considerations to keep in mind.

#. The most common difficulty is the accurate implementation of the Jacobian
   matrix. If the analytic Jacobian is not properly provided to the
   solver, this can hinder and many times prevent convergence of the method.
   When developing a new nonlinear least squares code, it often helps
   to compare the program output with the internally computed finite
   difference Jacobian and the user supplied analytic Jacobian. If there
   is a large difference in coefficients, it is likely the analytic
   Jacobian is incorrectly implemented.

#. If your code is having difficulty converging, the next thing to
   check is the starting point provided to the solver. The methods
   of this chapter are local methods, meaning if you provide a starting
   point far away from the true minimum, the method may converge to
   a local minimum or not converge at all. Sometimes it is possible
   to solve a linearized approximation to the nonlinear problem,
   and use the linear solution as the starting point to the nonlinear
   problem.

#. If the various parameters of the coefficient vector :math:`x`
   vary widely in magnitude, then the problem is said to be badly scaled.
   The methods of this chapter do attempt to automatically rescale
   the elements of :math:`x` to have roughly the same order of magnitude,
   but in extreme cases this could still cause problems for convergence.
   In these cases it is recommended for the user to scale their
   parameter vector :math:`x` so that each parameter spans roughly the
   same range, say :math:`[-1,1]`. The solution vector can be backscaled
   to recover the original units of the problem.

Examples
========

The following example programs demonstrate the nonlinear least
squares fitting capabilities.

Exponential Fitting Example
---------------------------

The following example program fits a weighted exponential model with
background to experimental data, :math:`Y = A \exp(-\lambda t) + b`. The
first part of the program sets up the functions :func:`expb_f` and
:func:`expb_df` to calculate the model and its Jacobian.  The appropriate
fitting function is given by,

.. math:: f_i = (A \exp(-\lambda t_i) + b) - y_i

where we have chosen :math:`t_i = i T / (N - 1)`, where :math:`N` is the number
of data points fitted, so that :math:`t_i \in [0, T]`. The Jacobian matrix :math:`J` is
the derivative of these functions with respect to the three parameters
(:math:`A`, :math:`\lambda`, :math:`b`).  It is given by,

.. only:: not texinfo

   .. math:: J_{ij} = {\partial f_i \over \partial x_j}

.. only:: texinfo

   ::

      J_{ij} = d f_i / d x_j

where :math:`x_0 = A`, :math:`x_1 = \lambda` and :math:`x_2 = b`.
The :math:`i`-th row of the Jacobian is therefore

.. only:: not texinfo

   .. math::

      J_{i\cdot} =
      \left(
        \begin{array}{ccc}
          \exp(-\lambda t_i) & -t_i A \exp(-\lambda t_i) & 1
        \end{array}
      \right)

.. only:: texinfo

   ::

      J(i,:) = [ \exp(-\lambda t_i) ; -t_i A \exp(-\lambda t_i) ; 1 ]

The main part of the program sets up a Levenberg-Marquardt solver and
some simulated random data. The data uses the known parameters
(5.0,1.5,1.0) combined with Gaussian noise (standard deviation = 0.1)
with a maximum time :math:`T = 3` and :math:`N = 100` timesteps.
The initial guess for the parameters is
chosen as (1.0, 1.0, 0.0). The iteration terminates when the relative
change in x is smaller than :math:`10^{-8}`, or when the magnitude of
the gradient falls below :math:`10^{-8}`. Here are the results of running
the program::

  iter  0: A = 1.0000, lambda = 1.0000, b = 0.0000, cond(J) =      inf, |f(x)| = 88.4448
  iter  1: A = 4.5109, lambda = 2.5258, b = 1.0704, cond(J) =  26.2686, |f(x)| = 24.0646
  iter  2: A = 4.8565, lambda = 1.7442, b = 1.1669, cond(J) =  23.7470, |f(x)| = 11.9797
  iter  3: A = 4.9356, lambda = 1.5713, b = 1.0767, cond(J) =  17.5849, |f(x)| = 10.7355
  iter  4: A = 4.8678, lambda = 1.4838, b = 1.0252, cond(J) =  16.3428, |f(x)| = 10.5000
  iter  5: A = 4.8118, lambda = 1.4481, b = 1.0076, cond(J) =  15.7925, |f(x)| = 10.4786
  iter  6: A = 4.7983, lambda = 1.4404, b = 1.0041, cond(J) =  15.5840, |f(x)| = 10.4778
  iter  7: A = 4.7967, lambda = 1.4395, b = 1.0037, cond(J) =  15.5396, |f(x)| = 10.4778
  iter  8: A = 4.7965, lambda = 1.4394, b = 1.0037, cond(J) =  15.5344, |f(x)| = 10.4778
  iter  9: A = 4.7965, lambda = 1.4394, b = 1.0037, cond(J) =  15.5339, |f(x)| = 10.4778
  iter 10: A = 4.7965, lambda = 1.4394, b = 1.0037, cond(J) =  15.5339, |f(x)| = 10.4778
  iter 11: A = 4.7965, lambda = 1.4394, b = 1.0037, cond(J) =  15.5339, |f(x)| = 10.4778
  summary from method 'trust-region/levenberg-marquardt'
  number of iterations: 11
  function evaluations: 16
  Jacobian evaluations: 12
  reason for stopping: small gradient
  initial |f(x)| = 88.444756
  final   |f(x)| = 10.477801
  chisq/dof = 1.1318
  A      = 4.79653 +/- 0.18704
  lambda = 1.43937 +/- 0.07390
  b      = 1.00368 +/- 0.03473
  status = success

The approximate values of the parameters are found correctly, and the
chi-squared value indicates a good fit (the chi-squared per degree of
freedom is approximately 1).  In this case the errors on the parameters
can be estimated from the square roots of the diagonal elements of the
covariance matrix. If the chi-squared value shows a poor fit (i.e.
:math:`\chi^2/(n-p) \gg 1`
then the error estimates obtained from the
covariance matrix will be too small.  In the example program the error estimates
are multiplied by :math:`\sqrt{\chi^2/(n-p)}`
in this case, a common way of increasing the
errors for a poor fit.  Note that a poor fit will result from the use
of an inappropriate model, and the scaled error estimates may then
be outside the range of validity for Gaussian errors.

Additionally, we see that the condition number of :math:`J(x)` stays
reasonably small throughout the iteration. This indicates we could
safely switch to the Cholesky solver for speed improvement,
although this particular system is too small to really benefit.

:numref:`fig_fit-exp` shows the fitted curve with the original data.

.. _fig_fit-exp:

.. figure:: /images/fit-exp.png
   :scale: 60%

   Exponential fitted curve with data

.. include:: examples/nlfit.c
   :code:

Geodesic Acceleration Example 1
-------------------------------

The following example program minimizes a modified Rosenbrock function,
which is characterized by a narrow canyon with steep walls. The
starting point is selected high on the canyon wall, so the solver
must first find the canyon bottom and then navigate to the minimum.
The problem is solved both with and without using geodesic acceleration
for comparison. The cost function is given by

.. only:: not texinfo

   .. math::

      \Phi(x) &= {1 \over 2} (f_1^2 + f_2^2) \\
      f_1 &= 100 \left( x_2 - x_1^2 \right) \\
      f_2 &= 1 - x_1

.. only:: texinfo

   ::

      Phi(x) = 1/2 (f1^2 + f2^2)
      f1 = 100 ( x2 - x1^2 )
      f2 = 1 - x1

The Jacobian matrix is

.. only:: not texinfo

   .. math::

      J =
      \left(
        \begin{array}{cc}
          {\partial f_1 \over \partial x_1} & {\partial f_1 \over \partial x_2} \\
          {\partial f_2 \over \partial x_1} & {\partial f_2 \over \partial x_2}
        \end{array}
      \right) =
      \left(
        \begin{array}{cc}
          -200 x_1 & 100 \\
          -1 & 0
        \end{array}
      \right)

.. only:: texinfo

   ::

      J = [ -200*x1 100  ]
          [   -1     0   ]

In order to use geodesic acceleration, the user must provide
the second directional derivative of each residual in the
velocity direction,
:math:`D_v^2 f_i = \sum_{\alpha\beta} v_{\alpha} v_{\beta} \partial_{\alpha} \partial_{\beta} f_i`.
The velocity vector :math:`v` is provided by the solver. For this example,
these derivatives are

.. only:: not texinfo

   .. math::

      f_{vv} =
      D_v^2
      \left(
        \begin{array}{c}
          f_1 \\
          f_2
        \end{array}
      \right) =
      \left(
        \begin{array}{c}
          -200 v_1^2 \\
          0
        \end{array}
      \right)

.. only:: texinfo

   ::

      fvv = [ -200 v1^2 ]
            [     0     ]

The solution of this minimization problem is

.. only:: not texinfo

   .. math::

      x^{*} &=
      \left(
        \begin{array}{c}
          1 \\
          1
        \end{array}
      \right) \\
      \Phi(x^{*}) &= 0

.. only:: texinfo

   ::

      x* = [ 1 ; 1 ]
      Phi(x*) = 0

The program output is shown below::

  === Solving system without acceleration ===
  NITER         = 53
  NFEV          = 56
  NJEV          = 54
  NAEV          = 0
  initial cost  = 2.250225000000e+04
  final cost    = 6.674986031430e-18
  final x       = (9.999999974165e-01, 9.999999948328e-01)
  final cond(J) = 6.000096055094e+02
  === Solving system with acceleration ===
  NITER         = 15
  NFEV          = 17
  NJEV          = 16
  NAEV          = 16
  initial cost  = 2.250225000000e+04
  final cost    = 7.518932873279e-19
  final x       = (9.999999991329e-01, 9.999999982657e-01)
  final cond(J) = 6.000097233278e+02

.. _fig_nlfit2:

.. figure:: /images/nlfit2.png

   Paths taken by solver for Rosenbrock function

We can see that enabling geodesic acceleration requires less
than a third of the number of Jacobian evaluations in order to locate
the minimum. The path taken by both methods is shown in :numref:`fig_nlfit2`.
The contours show the cost function
:math:`\Phi(x_1,x_2)`. We see that both methods quickly
find the canyon bottom, but the geodesic acceleration method
navigates along the bottom to the solution with significantly
fewer iterations.

The program is given below.

.. include:: examples/nlfit2.c
   :code:

Geodesic Acceleration Example 2
-------------------------------

The following example fits a set of data to a Gaussian model
using the Levenberg-Marquardt method with geodesic acceleration.
The cost function is

.. only:: not texinfo

   .. math::

      \Phi(x) &= {1 \over 2} \sum_i f_i^2 \\
      f_i &= y_i - Y(a,b,c;t_i)

.. only:: texinfo

   ::

      Phi(x) = 1/2 \sum_i f_i^2
      f_i = y_i - Y(a,b,c;t_i)

where :math:`y_i` is the measured data point at time :math:`t_i`, and
the model is specified by

.. only:: not texinfo

   .. math::

      Y(a,b,c;t) = a \exp{
      \left[
      -{1 \over 2}
      \left(
      { t - b \over c }
      \right)^2
      \right]
      }

.. only:: texinfo

   ::

      Y(a,b,c;t) = a exp(-1/2 ((t-b)/c)^2)

The parameters :math:`a,b,c` represent the amplitude, mean, and width of the Gaussian
respectively. The program below generates the :math:`y_i` data on :math:`[0,1]` using
the values :math:`a = 5`, :math:`b = 0.4`, :math:`c = 0.15` and adding random noise.
The :math:`i`-th row of the Jacobian is

.. only:: not texinfo

   .. math::

      J_{i,:} =
      \left(
        \begin{array}{ccc}
          {\partial f_i \over \partial a} & {\partial f_i \over \partial b} & {\partial f_i \over \partial c}
        \end{array}
      \right) =
      \left(
        \begin{array}{ccc}
          -e_i & -{a \over c} z_i e_i & -{a \over c} z_i^2 e_i
        \end{array}
      \right)

.. only:: texinfo

   ::

      J(i,:) = ( -e_i  -(a/c)*z_i*e_i  -(a/c)*z_i^2*e_i )

where

.. only:: not texinfo

   .. math::

      z_i &= { t_i - b \over c} \\
      e_i &= \exp{\left( -{1 \over 2} z_i^2 \right)}

.. only:: texinfo

   ::

      z_i = (t_i - b) / c
      e_i = \exp(-1/2 z_i^2)

In order to use geodesic acceleration, we need the second directional derivative
of the residuals in the velocity direction,
:math:`D_v^2 f_i = \sum_{\alpha\beta} v_{\alpha} v_{\beta} \partial_{\alpha} \partial_{\beta} f_i`,
where :math:`v` is provided by the solver. To compute this, it is helpful to make a table of
all second derivatives of the residuals :math:`f_i` with respect to each combination of model parameters.
This table is

.. only:: not texinfo

   .. math::

      \begin{array}{cccc}
        & {\partial \over \partial a} & {\partial \over \partial b} & {\partial \over \partial c} \cr
        {\partial \over \partial a} & 0 & -{z_i \over c} e_i & -{z_i^2 \over c} e_i \cr
        {\partial \over \partial b} & & {a \over c^2} \left( 1 - z_i^2 \right) e_i & {a \over c^2} z_i \left( 2 - z_i^2 \right) e_i \cr
        {\partial \over \partial c} & & & {a \over c^2} z_i^2 \left( 3 - z_i^2 \right) e_i
      \end{array}

The lower half of the table is omitted since it is symmetric. Then, the second directional derivative
of :math:`f_i` is

.. only:: not texinfo

   .. math:: D_v^2 f_i = v_a^2 \partial_a^2 f_i + 2 v_a v_b \partial_a \partial_b f_i + 2 v_a v_c \partial_a \partial_c f_i + v_b^2 \partial_b^2 f_i + 2 v_b v_c \partial_b \partial_c f_i + v_c^2 \partial_c^2 f_i

.. only:: texinfo

   ::

      D_v^2 f_i = v_a^2 (d/da)^2 f_i + 2 v_a v_b (d/da) (d/db) f_i + 2 v_a v_c (d/da) (d/dc) f_i + v_b^2 (d/db)^2 f_i + 2 v_b v_c (d/db) (d/dc) f_i + v_c^2 (d/dc)^2 f_i

The factors of 2 come from the symmetry of the mixed second partial derivatives.
The iteration is started using the initial guess :math:`a = 1, b = 0, c = 1`.
The program output is shown below::

  iter  0: a = 1.0000, b = 0.0000, c = 1.0000, |a|/|v| = 0.0000 cond(J) =      inf, |f(x)| = 35.4785
  iter  1: a = 1.5708, b = 0.5321, c = 0.5219, |a|/|v| = 0.3093 cond(J) =  29.0443, |f(x)| = 31.1042
  iter  2: a = 1.7387, b = 0.4040, c = 0.4568, |a|/|v| = 0.1199 cond(J) =   3.5256, |f(x)| = 28.7217
  iter  3: a = 2.2340, b = 0.3829, c = 0.3053, |a|/|v| = 0.3308 cond(J) =   4.5121, |f(x)| = 23.8074
  iter  4: a = 3.2275, b = 0.3952, c = 0.2243, |a|/|v| = 0.2784 cond(J) =   8.6499, |f(x)| = 15.6003
  iter  5: a = 4.3347, b = 0.3974, c = 0.1752, |a|/|v| = 0.2029 cond(J) =  15.1732, |f(x)| = 7.5908
  iter  6: a = 4.9352, b = 0.3992, c = 0.1536, |a|/|v| = 0.1001 cond(J) =  26.6621, |f(x)| = 4.8402
  iter  7: a = 5.0716, b = 0.3994, c = 0.1498, |a|/|v| = 0.0166 cond(J) =  34.6922, |f(x)| = 4.7103
  iter  8: a = 5.0828, b = 0.3994, c = 0.1495, |a|/|v| = 0.0012 cond(J) =  36.5422, |f(x)| = 4.7095
  iter  9: a = 5.0831, b = 0.3994, c = 0.1495, |a|/|v| = 0.0000 cond(J) =  36.6929, |f(x)| = 4.7095
  iter 10: a = 5.0831, b = 0.3994, c = 0.1495, |a|/|v| = 0.0000 cond(J) =  36.6975, |f(x)| = 4.7095
  iter 11: a = 5.0831, b = 0.3994, c = 0.1495, |a|/|v| = 0.0000 cond(J) =  36.6976, |f(x)| = 4.7095
  NITER         = 11
  NFEV          = 18
  NJEV          = 12
  NAEV          = 17
  initial cost  = 1.258724737288e+03
  final cost    = 2.217977560180e+01
  final x       = (5.083101559156e+00, 3.994484109594e-01, 1.494898e-01)
  final cond(J) = 3.669757713403e+01

We see the method converges after 11 iterations. For comparison the standard
Levenberg-Marquardt method requires 26 iterations and so the Gaussian fitting
problem benefits substantially from the geodesic acceleration correction. The
column marked :code:`|a|/|v|` above shows the ratio of the acceleration term
to the velocity term as the iteration progresses. Larger values of this
ratio indicate that the geodesic acceleration correction term is contributing
substantial information to the solver relative to the standard LM velocity step.

The data and fitted model are shown in :numref:`fig_nlfit2b`.

.. _fig_nlfit2b:

.. figure:: /images/nlfit2b.png

   Gaussian model fitted to data

The program is given below.

.. include:: examples/nlfit2b.c
   :code:

Comparing TRS Methods Example
-----------------------------

The following program compares all available nonlinear least squares
trust-region subproblem (TRS) methods on the Branin function, a common
optimization test problem. The cost function is

.. only:: not texinfo

   .. math::

      \Phi(x) &= {1 \over 2} (f_1^2 + f_2^2) \\
      f_1 &= x_2 + a_1 x_1^2 + a_2 x_1 + a_3 \\
      f_2 &= \sqrt{a_4} \sqrt{1 + (1 - a_5) \cos{x_1}}

.. only:: texinfo

   ::

      \Phi(x) &= 1/2 (f_1^2 + f_2^2)
      f_1 &= x_2 + a_1 x_1^2 + a_2 x_1 + a_3
      f_2 &= sqrt(a_4) sqrt(1 + (1 - a_5) cos(x_1))

with :math:`a_1 = -{5.1 \over 4 \pi^2}, a_2 = {5 \over \pi}, a_3 = -6, a_4 = 10, a_5 = {1 \over 8\pi}`.
There are three minima of this function in the range
:math:`(x_1,x_2) \in [-5,15] \times [-5,15]`. The program
below uses the starting point :math:`(x_1,x_2) = (6,14.5)`
and calculates the solution with all available nonlinear
least squares TRS methods. The program output is shown below::

  Method                    NITER  NFEV  NJEV  Initial Cost  Final cost   Final cond(J) Final x        
  levenberg-marquardt       20     27    21    1.9874e+02    3.9789e-01   6.1399e+07    (-3.14e+00, 1.23e+01)
  levenberg-marquardt+accel 27     36    28    1.9874e+02    3.9789e-01   1.4465e+07    (3.14e+00, 2.27e+00)
  dogleg                    23     64    23    1.9874e+02    3.9789e-01   5.0692e+08    (3.14e+00, 2.28e+00)
  double-dogleg             24     69    24    1.9874e+02    3.9789e-01   3.4879e+07    (3.14e+00, 2.27e+00)
  2D-subspace               23     54    24    1.9874e+02    3.9789e-01   2.5142e+07    (3.14e+00, 2.27e+00)

The first row of output above corresponds to standard Levenberg-Marquardt, while
the second row includes geodesic acceleration. We see that the standard LM method
converges to the minimum at :math:`(-\pi,12.275)` and also uses the least number
of iterations and Jacobian evaluations. All other methods converge to the minimum
:math:`(\pi,2.275)` and perform similarly in terms of number of Jacobian evaluations.
We see that :math:`J` is fairly ill-conditioned
at both minima, indicating that the QR (or SVD) solver is the best choice for this problem.
Since there are only two parameters in this optimization problem, we can easily
visualize the paths taken by each method, which are shown in :numref:`fig_nlfit3`.
The figure shows contours of the cost function :math:`\Phi(x_1,x_2)` which exhibits
three global minima in the range :math:`[-5,15] \times [-5,15]`. The paths taken
by each solver are shown as colored lines.

.. _fig_nlfit3:

.. figure:: /images/nlfit3.png
   :scale: 60%

   Paths taken for different TRS methods for the Branin function

The program is given below.

.. include:: examples/nlfit3.c
   :code:

Large Nonlinear Least Squares Example
-------------------------------------

The following program illustrates the large nonlinear least
squares solvers on a system with significant sparse structure
in the Jacobian. The cost function is

.. only:: not texinfo

   .. math::

      \Phi(x) &= {1 \over 2} \sum_{i=1}^{p+1} f_i^2 \\
      f_i &= \sqrt{\alpha} (x_i - 1), \quad 1 \le i \le p \\
      f_{p+1} &= ||x||^2 - {1 \over 4}

.. only:: texinfo

   ::

      \Phi(x) &= 1/2 \sum_{i=1}^{p+1} f_i^2
      f_i &= \sqrt{\alpha} (x_i - 1), 1 \le i \le p
      f_{p+1} &= ||x||^2 - 1/4

with :math:`\alpha = 10^{-5}`. The residual :math:`f_{p+1}` imposes a constraint on the :math:`p`
parameters :math:`x`, to ensure that :math:`||x||^2 \approx {1 \over 4}`.
The :math:`(p+1)`-by-:math:`p` Jacobian for this system is

.. only:: not texinfo

   .. math::

      J(x) =
      \left(
        \begin{array}{c}
          \sqrt{\alpha} I_p \\
          2 x^T
        \end{array}
      \right)

.. only:: texinfo

   ::

     J(x) = [ \sqrt{alpha} I_p; 2 x^T ]

and the normal equations matrix is

.. math:: J^T J = \alpha I_p + 4 x x^T

Finally, the second directional derivative of :math:`f` for the
geodesic acceleration method is

.. only:: not texinfo

   .. math::

      f_{vv} = D_v^2 f =
      \left(
        \begin{array}{c}
          0 \\
          2 ||v||^2
        \end{array}
      \right)

.. only:: texinfo

   ::

      fvv = [     0     ]
            [ 2 ||v||^2 ]

Since the upper :math:`p`-by-:math:`p` block of :math:`J` is diagonal,
this sparse structure should be exploited in the nonlinear solver.
For comparison, the following program solves the system for :math:`p = 2000`
using the dense direct Cholesky solver based on the normal equations matrix
:math:`J^T J`, as well as the iterative Steihaug-Toint solver, based on
sparse matrix-vector products :math:`J u` and :math:`J^T u`. The
program output is shown below::

  Method                    NITER NFEV NJUEV NJTJEV NAEV Init Cost  Final cost cond(J) Final |x|^2 Time (s)  
  levenberg-marquardt       25    31   26    26     0    7.1218e+18 1.9555e-02 447.50  2.5044e-01  46.28
  levenberg-marquardt+accel 22    23   45    23     22   7.1218e+18 1.9555e-02 447.64  2.5044e-01  33.92
  dogleg                    37    87   36    36     0    7.1218e+18 1.9555e-02 447.59  2.5044e-01  56.05
  double-dogleg             35    88   34    34     0    7.1218e+18 1.9555e-02 447.62  2.5044e-01  52.65
  2D-subspace               37    88   36    36     0    7.1218e+18 1.9555e-02 447.71  2.5044e-01  59.75
  steihaug-toint            35    88   345   0      0    7.1218e+18 1.9555e-02 inf     2.5044e-01  0.09

The first five rows use methods based on factoring the dense :math:`J^T J` matrix
while the last row uses the iterative Steihaug-Toint method. While the number
of Jacobian matrix-vector products (NJUEV) is less for the dense methods, the added time
to construct and factor the :math:`J^T J` matrix (NJTJEV) results in a much larger runtime than the
iterative method (see last column).

The program is given below.

.. include:: examples/nlfit4.c
   :code:

References and Further Reading
==============================

The following publications are relevant to the algorithms described
in this section,

* J.J. |More|, *The Levenberg-Marquardt Algorithm: Implementation and
  Theory*, Lecture Notes in Mathematics, v630 (1978), ed G. Watson.

* H. B. Nielsen, "Damping Parameter in Marquardt's Method",
  IMM Department of Mathematical Modeling, DTU, Tech. Report IMM-REP-1999-05
  (1999).

* K. Madsen and H. B. Nielsen, "Introduction to Optimization and Data
  Fitting", IMM Department of Mathematical Modeling, DTU, 2010.

* J. E. Dennis and R. B. Schnabel, Numerical Methods for Unconstrained
  Optimization and Nonlinear Equations, SIAM, 1996.

* M. K. Transtrum, B. B. Machta, and J. P. Sethna,
  Geometry of nonlinear least squares with applications to sloppy models and optimization,
  Phys. Rev. E 83, 036701, 2011.

* M. K. Transtrum and J. P. Sethna, Improvements to the Levenberg-Marquardt
  algorithm for nonlinear least-squares minimization, arXiv:1201.5885, 2012.

* J.J. |More|, B.S. Garbow, K.E. Hillstrom, "Testing Unconstrained
  Optimization Software", ACM Transactions on Mathematical Software, Vol
  7, No 1 (1981), p 17--41.

* H. B. Nielsen, "UCTP Test Problems for Unconstrained Optimization",
  IMM Department of Mathematical Modeling, DTU, Tech. Report IMM-REP-2000-17
  (2000).
