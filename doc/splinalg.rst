.. index::
   single: sparse linear algebra
   single: linear algebra, sparse

*********************
Sparse Linear Algebra
*********************

This chapter describes functions for solving sparse linear systems
of equations. The library provides linear algebra routines which
operate directly on the :type:`gsl_spmatrix` and :type:`gsl_vector`
objects.

The functions described in this chapter are declared in the header file
:file:`gsl_splinalg.h`.

.. index::
   single: sparse linear algebra, overview

Overview
========

This chapter is primarily concerned with the solution of the
linear system

.. math:: A x = b

where :math:`A` is a general square :math:`n`-by-:math:`n` non-singular
sparse matrix, :math:`x` is an unknown :math:`n`-by-:math:`1` vector, and
:math:`b` is a given :math:`n`-by-1 right hand side vector. There exist
many methods for solving such sparse linear systems, which broadly
fall into either direct or iterative categories. Direct methods include
LU and QR decompositions, while iterative methods start with an
initial guess for the vector :math:`x` and update the guess through
iteration until convergence. GSL does not currently provide any
direct sparse solvers.

.. index::
   single: sparse matrices, iterative solvers
   single: sparse linear algebra, iterative solvers
   single: sparse, iterative solvers

Sparse Iterative Solvers
========================

Overview
--------

Many practical iterative methods of solving large :math:`n`-by-:math:`n`
sparse linear systems involve projecting an approximate solution for
:data:`x` onto a subspace of :math:`{\bf R}^n`. If we define a :math:`m`-dimensional
subspace :math:`{\cal K}` as the subspace of approximations to the solution
:data:`x`, then :math:`m` constraints must be imposed to determine
the next approximation. These :math:`m` constraints define another
:math:`m`-dimensional subspace denoted by :math:`{\cal L}`. The
subspace dimension :math:`m` is typically chosen to be much smaller than
:math:`n` in order to reduce the computational
effort needed to generate the next approximate solution vector.
The many iterative algorithms which exist differ mainly
in their choice of :math:`{\cal K}` and :math:`{\cal L}`.

Types of Sparse Iterative Solvers
---------------------------------

The sparse linear algebra library provides the following types
of iterative solvers:

.. type:: gsl_splinalg_itersolve_type

   .. index:: gmres

   .. var:: gsl_splinalg_itersolve_type * gsl_splinalg_itersolve_gmres

      This specifies the Generalized Minimum Residual Method (GMRES).
      This is a projection method using :math:`{\cal K} = {\cal K}_m`
      and :math:`{\cal L} = A {\cal K}_m` where :math:`{\cal K}_m` is
      the :math:`m`-th Krylov subspace

      .. only:: not texinfo

         .. math:: {\cal K}_m = span \left\{ r_0, A r_0, ..., A^{m-1} r_0 \right\}

      .. only:: texinfo

         ::

            K_m = span( r_0, A r_0, ..., A^(m-1) r_0)

      and :math:`r_0 = b - A x_0` is the residual vector of the initial guess
      :math:`x_0`. If :math:`m` is set equal to :math:`n`, then the Krylov
      subspace is :math:`{\bf R}^n` and GMRES will provide the exact solution
      :data:`x`.  However, the goal is for the method to arrive at a very good
      approximation to :data:`x` using a much smaller subspace :math:`{\cal K}_m`. By
      default, the GMRES method selects :math:`m = MIN(n,10)` but the user
      may specify a different value for :math:`m`. The GMRES storage
      requirements grow as :math:`O(n(m+1))` and the number of flops
      grow as :math:`O(4 m^2 n - 4 m^3 / 3)`.

      In the below function :func:`gsl_splinalg_itersolve_iterate`, one
      GMRES iteration is defined as projecting the approximate solution
      vector onto each Krylov subspace :math:`{\cal K}_1, ..., {\cal K}_m`,
      and so :math:`m` can be kept smaller by "restarting" the method
      and calling :func:`gsl_splinalg_itersolve_iterate` multiple times,
      providing the updated approximation :data:`x` to each new call. If
      the method is not adequately converging, the user may try increasing
      the parameter :math:`m`.

      GMRES is considered a robust general purpose iterative solver, however
      there are cases where the method stagnates if the matrix is not
      positive-definite and fails to reduce the residual until the very last
      projection onto the subspace :math:`{\cal K}_n = {\bf R}^n`. In these
      cases, preconditioning the linear system can help, but GSL does not
      currently provide any preconditioners.

Iterating the Sparse Linear System
----------------------------------

The following functions are provided to allocate storage for the
sparse linear solvers and iterate the system to a solution.

.. function:: gsl_splinalg_itersolve * gsl_splinalg_itersolve_alloc (const gsl_splinalg_itersolve_type * T, const size_t n, const size_t m)

   This function allocates a workspace for the iterative solution of
   :data:`n`-by-:data:`n` sparse matrix systems. The iterative solver type
   is specified by :data:`T`. The argument :data:`m` specifies the size
   of the solution candidate subspace :math:`{\cal K}_m`. The dimension
   :data:`m` may be set to 0 in which case a reasonable default value is used.

.. function:: void gsl_splinalg_itersolve_free (gsl_splinalg_itersolve * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: const char * gsl_splinalg_itersolve_name (const gsl_splinalg_itersolve * w)

   This function returns a string pointer to the name of the solver.

.. function:: int gsl_splinalg_itersolve_iterate (const gsl_spmatrix * A, const gsl_vector * b, const double tol, gsl_vector * x, gsl_splinalg_itersolve * w)

   This function performs one iteration of the iterative method for
   the sparse linear system specfied by the matrix :data:`A`, right hand
   side vector :data:`b` and solution vector :data:`x`. On input, :data:`x`
   must be set to an initial guess for the solution. On output,
   :data:`x` is updated to give the current solution estimate. The
   parameter :data:`tol` specifies the relative tolerance between the residual
   norm and norm of :data:`b` in order to check for convergence.
   When the following condition is satisfied:

   .. only:: not texinfo

      .. math:: || A x - b || \le tol \times || b ||

   .. only:: texinfo

      ::

         || A x - b || <= tol * || b ||

   the method has converged, the function returns :macro:`GSL_SUCCESS` and
   the final solution is provided in :data:`x`. Otherwise, the function
   returns :macro:`GSL_CONTINUE` to signal that more iterations are
   required. Here, :math:`|| \cdot ||` represents the Euclidean norm.
   The input matrix :data:`A` may be in triplet or compressed format.

.. function:: double gsl_splinalg_itersolve_normr (const gsl_splinalg_itersolve * w)

   This function returns the current residual norm
   :math:`||r|| = ||A x - b||`, which is updated after each call to
   :func:`gsl_splinalg_itersolve_iterate`.

.. index::
   single: sparse linear algebra, examples

Examples
========

This example program demonstrates the sparse linear algebra routines on
the solution of a simple 1D Poisson equation on :math:`[0,1]`:

.. only:: not texinfo

   .. math:: {d^2 u(x) \over dx^2} = f(x) = -\pi^2 \sin{(\pi x)}

.. only:: texinfo

   ::

      u''(x) = f(x) = -\pi^2 \sin(\pi x)

with boundary conditions :math:`u(0) = u(1) = 0`. The analytic solution of
this simple problem is :math:`u(x) = \sin{\pi x}`. We will solve this
problem by finite differencing the left hand side to give

.. only:: not texinfo

   .. math:: {1 \over h^2} \left( u_{i+1} - 2 u_i + u_{i-1} \right) = f_i

.. only:: texinfo

   ::

      1/h^2 ( u_(i+1) - 2 u_i + u_(i-1) ) = f_i

Defining a grid of :math:`N` points, :math:`h = 1/(N-1)`. In the finite
difference equation above, :math:`u_0 = u_{N-1} = 0` are known from
the boundary conditions, so we will only put the equations for
:math:`i = 1, ..., N-2` into the matrix system. The resulting
:math:`(N-2) \times (N-2)` matrix equation is

.. only:: not texinfo

   .. math::

      {1 \over h^2}
      \left(
        \begin{array}{cccccc}
          -2 & 1 & 0 & 0 & \ldots & 0 \\
          1 & -2 & 1 & 0 & \ldots & 0 \\
          0 & 1 & -2 & 1 & \ldots & 0 \\
          \vdots & \vdots & \ddots & \ddots & \ddots & \vdots \\
          0 & \ldots & \ldots & 1 & -2 & 1 \\
          0 & \ldots & \ldots & \ldots & 1 & -2
        \end{array}
      \right)
      \left(
        \begin{array}{c}
          u_1 \\
          u_2 \\
          u_3 \\
          \vdots \\
          u_{N-3} \\
          u_{N-2}
        \end{array}
      \right) =
      \left(
        \begin{array}{c}
          f_1 \\
          f_2 \\
          f_3 \\
          \vdots \\
          f_{N-3} \\
          f_{N-2}
        \end{array}
      \right)

An example program which constructs and solves this system is given below.
The system is solved using the iterative GMRES solver. Here is
the output of the program::

  iter 0 residual = 4.297275996844e-11
  Converged

showing that the method converged in a single iteration.
The calculated solution is shown in :numref:`fig_splinalg-poisson`.

.. _fig_splinalg-poisson:

.. figure:: /images/sparse_poisson.png

   Solution of PDE

.. include:: examples/poisson.c
   :code:

.. index::
   single: sparse linear algebra, references

References and Further Reading
==============================

The implementation of the GMRES iterative solver closely follows
the publications

* H. F. Walker, Implementation of the GMRES method using
  Householder transformations, SIAM J. Sci. Stat. Comput.
  9(1), 1988.

* Y. Saad, Iterative methods for sparse linear systems, 2nd edition,
  SIAM, 2003.
