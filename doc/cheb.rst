.. index::
   single: Chebyshev series
   single: fitting, using Chebyshev polynomials
   single: interpolation, using Chebyshev polynomials

************************
Chebyshev Approximations
************************

This chapter describes routines for computing Chebyshev approximations
to univariate functions.  A Chebyshev approximation is a truncation of
the series :math:`f(x) = \sum c_n T_n(x)`, where the Chebyshev
polynomials :math:`T_n(x) = \cos(n \arccos x)` provide an orthogonal
basis of polynomials on the interval :math:`[-1,1]` with the weight
function :math:`1 / \sqrt{1-x^2}`.
The first few Chebyshev polynomials are,
:math:`T_0(x) = 1`, :math:`T_1(x) = x`, :math:`T_2(x) = 2 x^2 - 1`.
For further information see Abramowitz & Stegun, Chapter 22. 

The functions described in this chapter are declared in the header file
:file:`gsl_chebyshev.h`.

Definitions
===========

.. type:: gsl_cheb_series

   A Chebyshev series  is stored using the following structure::

      typedef struct
      {
        double * c;   /* coefficients  c[0] .. c[order] */
        int order;    /* order of expansion             */
        double a;     /* lower interval point           */
        double b;     /* upper interval point           */
        ...
      } gsl_cheb_series

The approximation is made over the range :math:`[a,b]` using
:data:`order` + 1 terms, including the coefficient :math:`c[0]`.  The series
is computed using the following convention,

.. only:: not texinfo

   .. math:: f(x) = {c_0 \over 2} + \sum_{n=1} c_n T_n(x)

.. only:: texinfo

   ::

      f(x) = (c_0 / 2) + \sum_{n=1} c_n T_n(x)

which is needed when accessing the coefficients directly.

Creation and Calculation of Chebyshev Series
============================================

.. function:: gsl_cheb_series * gsl_cheb_alloc (const size_t n)

   This function allocates space for a Chebyshev series of order :data:`n`
   and returns a pointer to a new :type:`gsl_cheb_series` struct.

.. function:: void gsl_cheb_free (gsl_cheb_series * cs)

   This function frees a previously allocated Chebyshev series :data:`cs`.

.. function:: int gsl_cheb_init (gsl_cheb_series * cs, const gsl_function * f, const double a, const double b)

   This function computes the Chebyshev approximation :data:`cs` for the
   function :data:`f` over the range :math:`(a,b)` to the previously specified
   order.  The computation of the Chebyshev approximation is an
   :math:`O(n^2)` process, and requires :math:`n` function evaluations.

Auxiliary Functions
===================
The following functions provide information about an existing
Chebyshev series.

.. function:: size_t gsl_cheb_order (const gsl_cheb_series * cs)

   This function returns the order of Chebyshev series :data:`cs`.

.. function:: size_t gsl_cheb_size (const gsl_cheb_series * cs)
              double * gsl_cheb_coeffs (const gsl_cheb_series * cs)

   These functions return the size of the Chebyshev coefficient array
   :code:`c[]` and a pointer to its location in memory for the Chebyshev
   series :data:`cs`.

Chebyshev Series Evaluation
===========================

.. function:: double gsl_cheb_eval (const gsl_cheb_series * cs, double x)

   This function evaluates the Chebyshev series :data:`cs` at a given point :data:`x`.

.. function:: int gsl_cheb_eval_err (const gsl_cheb_series * cs, const double x, double * result, double * abserr)

   This function computes the Chebyshev series :data:`cs` at a given point
   :data:`x`, estimating both the series :data:`result` and its absolute error
   :data:`abserr`.  The error estimate is made from the first neglected term
   in the series.

.. function:: double gsl_cheb_eval_n (const gsl_cheb_series * cs, size_t order, double x)

   This function evaluates the Chebyshev series :data:`cs` at a given point
   :data:`x`, to (at most) the given order :data:`order`.

.. function:: int gsl_cheb_eval_n_err (const gsl_cheb_series * cs, const size_t order, const double x, double * result, double * abserr)

   This function evaluates a Chebyshev series :data:`cs` at a given point
   :data:`x`, estimating both the series :data:`result` and its absolute error
   :data:`abserr`, to (at most) the given order :data:`order`.  The error
   estimate is made from the first neglected term in the series.

.. @deftypefun double gsl_cheb_eval_mode (const gsl_cheb_series * cs, double x, gsl_mode_t mode)
.. @end deftypefun

.. @deftypefun int gsl_cheb_eval_mode_err (const gsl_cheb_series * cs, const double x, gsl_mode_t mode, double * result, double * abserr)
.. Evaluate a Chebyshev series at a given point, using the default
.. order for double precision mode(s) and the single precision
.. order for other modes.
.. @end deftypefun

Derivatives and Integrals
=========================

The following functions allow a Chebyshev series to be differentiated or
integrated, producing a new Chebyshev series.  Note that the error
estimate produced by evaluating the derivative series will be
underestimated due to the contribution of higher order terms being
neglected.

.. function:: int gsl_cheb_calc_deriv (gsl_cheb_series * deriv, const gsl_cheb_series * cs)

   This function computes the derivative of the series :data:`cs`, storing
   the derivative coefficients in the previously allocated :data:`deriv`.
   The two series :data:`cs` and :data:`deriv` must have been allocated with
   the same order.

.. function:: int gsl_cheb_calc_integ (gsl_cheb_series * integ, const gsl_cheb_series * cs)

   This function computes the integral of the series :data:`cs`, storing the
   integral coefficients in the previously allocated :data:`integ`.  The two
   series :data:`cs` and :data:`integ` must have been allocated with the same
   order.  The lower limit of the integration is taken to be the left hand
   end of the range :data:`a`.

Examples
========

The following example program computes Chebyshev approximations to a
step function.  This is an extremely difficult approximation to make,
due to the discontinuity, and was chosen as an example where
approximation error is visible.  For smooth functions the Chebyshev
approximation converges extremely rapidly and errors would not be
visible.

.. include:: examples/cheb.c
   :code:

:numref:`fig_cheb` shows output from the program with the original function, 10-th order
approximation and 40-th order approximation, all sampled at intervals of
0.001 in :math:`x`.

.. _fig_cheb:

.. figure:: /images/cheb.png
   :scale: 60%

   Chebyshev approximations to a step function

References and Further Reading
==============================

The following paper describes the use of Chebyshev series,

* R. Broucke, "Ten Subroutines for the Manipulation of Chebyshev Series
  [C1] (Algorithm 446)". *Communications of the ACM* 16(4), 254--256
  (1973)
