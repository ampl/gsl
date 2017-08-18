.. index::
   single: basis splines, B-splines
   single: splines, basis

.. _chap_basis-splines:

*************
Basis Splines
*************

This chapter describes functions for the computation of smoothing
basis splines (B-splines). A smoothing spline differs from an
interpolating spline in that the resulting curve is not required to
pass through each datapoint.  For information about
interpolating splines, see :ref:`sec_interpolation`.

The header file :file:`gsl_bspline.h` contains the prototypes for the
bspline functions and related declarations.

.. index::
   single: basis splines, overview

Overview
========

B-splines are commonly used as basis functions to fit smoothing
curves to large data sets. To do this, the abscissa axis is
broken up into some number of intervals, where the endpoints
of each interval are called *breakpoints*. These breakpoints
are then converted to *knots* by imposing various continuity
and smoothness conditions at each interface. Given a nondecreasing
knot vector

.. math:: t = \{t_0, t_1, \dots, t_{n+k-1}\}

the :math:`n` basis splines of order :math:`k` are defined by

.. only:: not texinfo

   .. math::

      B_{i,1}(x) &=
        \left\{
          \begin{array}{cc}
            1, & t_i \le x < t_{i+1} \\
            0, & else
          \end{array}
        \right. \\
      B_{i,k}(x) &= {(x - t_i) \over (t_{i+k-1} - t_i)} B_{i,k-1}(x) +
                    {(t_{i+k} - x) \over (t_{i+k} - t_{i+1})} B_{i+1,k-1}(x)

.. only:: texinfo

   ::

      B_(i,1)(x) = (1, t_i <= x < t_(i+1)
                   (0, else
      B_(i,k)(x) = [(x - t_i)/(t_(i+k-1) - t_i)] B_(i,k-1)(x)
                    + [(t_(i+k) - x)/(t_(i+k) - t_(i+1))] B_(i+1,k-1)(x)

for :math:`i = 0, \ldots, n-1`. The common case of cubic B-splines
is given by :math:`k = 4`. The above recurrence relation can be
evaluated in a numerically stable way by the de Boor algorithm.

If we define appropriate knots on an interval :math:`[a,b]` then
the B-spline basis functions form a complete set on that interval.
Therefore we can expand a smoothing function as

.. math:: f(x) = \sum_{i=0}^{n-1} c_i B_{i,k}(x)

given enough :math:`(x_j, f(x_j))` data pairs. The coefficients
:math:`c_i` can be readily obtained from a least-squares fit.

.. index::
   single: basis splines, initializing

Initializing the B-splines solver
=================================

.. type:: gsl_bspline_workspace

   The computation of B-spline functions requires a preallocated
   workspace.

.. function:: gsl_bspline_workspace * gsl_bspline_alloc (const size_t k, const size_t nbreak)

   This function allocates a workspace for computing B-splines of order
   :data:`k`. The number of breakpoints is given by :data:`nbreak`. This
   leads to :math:`n = nbreak + k - 2` basis functions. Cubic B-splines
   are specified by :math:`k = 4`. The size of the workspace is
   :math:`O(2k^2 + 5k + nbreak)`.

.. function:: void gsl_bspline_free (gsl_bspline_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. index::
   single: knots, basis splines

Constructing the knots vector
=============================

.. function:: int gsl_bspline_knots (const gsl_vector * breakpts, gsl_bspline_workspace * w)

   This function computes the knots associated with the given breakpoints
   and stores them internally in :code:`w->knots`.

.. function:: int gsl_bspline_knots_uniform (const double a, const double b, gsl_bspline_workspace * w)

   This function assumes uniformly spaced breakpoints on :math:`[a,b]`
   and constructs the corresponding knot vector using the previously
   specified :data:`nbreak` parameter. The knots are stored in
   :code:`w->knots`.

.. index::
   single: basis splines, evaluation

Evaluation of B-splines
=======================

.. function:: int gsl_bspline_eval (const double x, gsl_vector * B, gsl_bspline_workspace * w)

   This function evaluates all B-spline basis functions at the position
   :data:`x` and stores them in the vector :data:`B`, so that the :math:`i`-th element
   is :math:`B_i(x)`. The vector :data:`B` must be of length
   :math:`n = nbreak + k - 2`.  This value may also be obtained by calling
   :func:`gsl_bspline_ncoeffs`.
   Computing all the basis functions at once is more efficient than
   computing them individually, due to the nature of the defining
   recurrence relation.

.. function:: int gsl_bspline_eval_nonzero (const double x, gsl_vector * Bk, size_t * istart, size_t * iend, gsl_bspline_workspace * w)

   This function evaluates all potentially nonzero B-spline basis
   functions at the position :data:`x` and stores them in the vector :data:`Bk`, so
   that the :math:`i`-th element is :math:`B_{(istart+i)}(x)`.
   The last element of :data:`Bk` is :math:`B_{iend}(x)`.
   The vector :data:`Bk` must be
   of length :math:`k`.  By returning only the nonzero basis functions,
   this function
   allows quantities involving linear combinations of the :math:`B_i(x)`
   to be computed without unnecessary terms
   (such linear combinations occur, for example,
   when evaluating an interpolated function).

.. function:: size_t gsl_bspline_ncoeffs (gsl_bspline_workspace * w)

   This function returns the number of B-spline coefficients given by
   :math:`n = nbreak + k - 2`.

.. index::
   single: basis splines, derivatives

Evaluation of B-spline derivatives
==================================

.. function:: int gsl_bspline_deriv_eval (const double x, const size_t nderiv, gsl_matrix * dB, gsl_bspline_workspace * w)

   This function evaluates all B-spline basis function derivatives of orders
   :math:`0` through :data:`nderiv` (inclusive) at the position :data:`x`
   and stores them in the matrix :data:`dB`.  The :math:`(i,j)`-th element of :data:`dB`
   is :math:`d^jB_i(x)/dx^j`.  The matrix :data:`dB` must be
   of size :math:`n = nbreak + k - 2` by :math:`nderiv + 1`.
   The value :math:`n` may also be obtained
   by calling :func:`gsl_bspline_ncoeffs`.  Note that function evaluations
   are included as the zeroth order derivatives in :data:`dB`.
   Computing all the basis function derivatives at once is more efficient
   than computing them individually, due to the nature of the defining
   recurrence relation.

.. function:: int gsl_bspline_deriv_eval_nonzero (const double x, const size_t nderiv, gsl_matrix * dB, size_t * istart, size_t * iend, gsl_bspline_workspace * w)

   This function evaluates all potentially nonzero B-spline basis function
   derivatives of orders :math:`0` through :data:`nderiv` (inclusive) at
   the position :data:`x` and stores them in the matrix :data:`dB`.  The
   :math:`(i,j)`-th element of :data:`dB` is :math:`d^jB_{(istart+i)}(x)/dx^j`.
   The last row of :data:`dB` contains :math:`d^jB_{iend}(x)/dx^j`.
   The matrix :data:`dB` must be
   of size :math:`k` by at least :math:`nderiv + 1`.  Note that function
   evaluations are included as the zeroth order derivatives in :data:`dB`.
   By returning only the nonzero basis functions, this function allows
   quantities involving linear combinations of the :math:`B_i(x)` and
   their derivatives to be computed without unnecessary terms.

.. index::
   single: basis splines, Greville abscissae
   single: basis splines, Marsden-Schoenberg points

Working with the Greville abscissae
===================================

The Greville abscissae are defined to be the mean location of :math:`k-1`
consecutive knots in the knot vector for each basis spline function of order
:math:`k`.  With the first and last knots in the :type:`gsl_bspline_workspace`
knot vector excluded, there are :func:`gsl_bspline_ncoeffs` Greville abscissae
for any given B-spline basis.  These values are often used in B-spline
collocation applications and may also be called Marsden-Schoenberg points.

.. function:: double gsl_bspline_greville_abscissa (size_t i, gsl_bspline_workspace * w)

   Returns the location of the :math:`i`-th Greville abscissa for the given
   B-spline basis.  For the ill-defined case when :math:`k = 1`, the implementation
   chooses to return breakpoint interval midpoints.

.. See https://savannah.gnu.org/bugs/index.php?34361
.. @deftypefun int gsl_bspline_knots_greville (const gsl_vector * abscissae, gsl_bspline_workspace * w, double * abserr);
.. Given target Greville abscissae values in :data:`abscissae` and a workspace
.. :data:`w` where @code{abscissae->size == gsl_bspline_ncoeffs(w)}, this functions
.. computes and stores the knots required for the workspace to best approximate
.. the target abscissae.  The approximation is optimal in that the first and last
.. values in :data:`abscissae` are preserved exactly while the 2-norm of the error
.. in any other abscissae is minimized.  If not-@code{NULL}, the sum of the
.. absolute approximation errors over each abscissa is returned in :data:`abserr`.
..
.. The workspace order must satisfy :math:`k > 1` and :data:`abscissae` should be
.. monotonically increasing.  Beware that when @code{w->nbreak} is small relative
.. to @code{w->k} the best approximation may still be of poor quality for
.. non-uniformly spaced :data:`abscissae`.  This function has memory and runtime
.. overhead that scales like a QR-based linear least squares solution on a
.. @code{(abscissae->size - 2)} by @code{(w->nbreak - 2)} problem.
.. @end deftypefun

.. index::
   single: basis splines, examples

Examples
========

The following program computes a linear least squares fit to data using
cubic B-spline basis functions with uniform breakpoints. The data is
generated from the curve :math:`y(x) = \cos{(x)} \exp{(-x/10)}` on
the interval :math:`[0, 15]` with Gaussian noise added.

.. include:: examples/bspline.c
   :code:

The output is shown below::

  $ ./a.out > bspline.txt
  chisq/dof = 1.118217e+00, Rsq = 0.989771

The data and fitted model are shown in :numref:`fig_bspline`.

.. _fig_bspline:

.. figure:: /images/bspline.png
   :scale: 60%

   Data (black) and fitted model (red)

References and Further Reading
==============================

Further information on the algorithms described in this section can be
found in the following book,

* C. de Boor, *A Practical Guide to Splines* (1978), Springer-Verlag,
  ISBN 0-387-90356-9.

Further information of Greville abscissae and B-spline collocation
can be found in the following paper,

* Richard W. Johnson, Higher order B-spline collocation at the Greville
  abscissae.  *Applied Numerical Mathematics*. vol.: 52, 2005, 63--75.

A large collection of B-spline routines is available in the
PPPACK library available at http://www.netlib.org/pppack,
which is also part of SLATEC.
