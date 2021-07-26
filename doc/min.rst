.. index::
   single: optimization, see minimization
   single: maximization, see minimization
   single: minimization, one-dimensional
   single: finding minima
   single: nonlinear functions, minimization

****************************
One Dimensional Minimization
****************************

This chapter describes routines for finding minima of arbitrary
one-dimensional functions.  The library provides low level components
for a variety of iterative minimizers and convergence tests.  These can be
combined by the user to achieve the desired solution, with full access
to the intermediate steps of the algorithms.  Each class of methods uses
the same framework, so that you can switch between minimizers at runtime
without needing to recompile your program.  Each instance of a minimizer
keeps track of its own state, allowing the minimizers to be used in
multi-threaded programs.

The header file :file:`gsl_min.h` contains prototypes for the
minimization functions and related declarations.  To use the minimization
algorithms to find the maximum of a function simply invert its sign.

.. index::
   single: minimization, overview

Overview
========

The minimization algorithms begin with a bounded region known to contain
a minimum.  The region is described by a lower bound :math:`a` and an
upper bound :math:`b`, with an estimate of the location of the minimum
:math:`x`, as shown in :numref:`fig_min-interval`.

.. _fig_min-interval:

.. figure:: /images/min-interval.png
   :scale: 60%

   Function with lower and upper bounds with an estimate of the minimum.

The value of the function at :math:`x` must be less than the value of the
function at the ends of the interval,

.. math:: f(a) > f(x) < f(b)

This condition guarantees that a minimum is contained somewhere within
the interval.  On each iteration a new point :math:`x'` is selected using
one of the available algorithms.  If the new point is a better estimate
of the minimum, i.e.: where :math:`f(x') < f(x)`, then the current
estimate of the minimum :math:`x` is updated.  The new point also allows
the size of the bounded interval to be reduced, by choosing the most
compact set of points which satisfies the constraint :math:`f(a) > f(x) < f(b)`.
The interval is reduced until it encloses the true minimum to a
desired tolerance.  This provides a best estimate of the location of the
minimum and a rigorous error estimate.

Several bracketing algorithms are available within a single framework.
The user provides a high-level driver for the algorithm, and the
library provides the individual functions necessary for each of the
steps.  There are three main phases of the iteration.  The steps are,

* initialize minimizer state, :data:`s`, for algorithm :data:`T`
* update :data:`s` using the iteration :data:`T`
* test :data:`s` for convergence, and repeat iteration if necessary

The state for the minimizers is held in a :type:`gsl_min_fminimizer`
struct.  The updating procedure uses only function evaluations (not
derivatives).

.. index::
   single: minimization, caveats

Caveats
=======

Note that minimization functions can only search for one minimum at a
time.  When there are several minima in the search area, the first
minimum to be found will be returned; however it is difficult to predict
which of the minima this will be. *In most cases, no error will be
reported if you try to find a minimum in an area where there is more
than one.*

With all minimization algorithms it can be difficult to determine the
location of the minimum to full numerical precision.  The behavior of the
function in the region of the minimum :math:`x^*` can be approximated by
a Taylor expansion,

.. only:: not texinfo

   .. math:: y = f(x^*) + {1 \over 2} f''(x^*) (x - x^*)^2

.. only:: texinfo

   ::

      y = f(x^*) + (1/2) f''(x^*) (x - x^*)^2

and the second term of this expansion can be lost when added to the
first term at finite precision.  This magnifies the error in locating
:math:`x^*`, making it proportional to :math:`\sqrt \epsilon` (where
:math:`\epsilon` is the relative accuracy of the floating point numbers).
For functions with higher order minima, such as :math:`x^4`, the
magnification of the error is correspondingly worse.  The best that can
be achieved is to converge to the limit of numerical accuracy in the
function values, rather than the location of the minimum itself.

Initializing the Minimizer
==========================

.. type:: gsl_min_fminimizer

   This is a workspace for minimizing functions.

.. function:: gsl_min_fminimizer * gsl_min_fminimizer_alloc (const gsl_min_fminimizer_type * T)

   This function returns a pointer to a newly allocated instance of a
   minimizer of type :data:`T`.  For example, the following code
   creates an instance of a golden section minimizer::

      const gsl_min_fminimizer_type * T = gsl_min_fminimizer_goldensection;
      gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);

   If there is insufficient memory to create the minimizer then the function
   returns a null pointer and the error handler is invoked with an error
   code of :macro:`GSL_ENOMEM`.

.. function:: int gsl_min_fminimizer_set (gsl_min_fminimizer * s, gsl_function * f, double x_minimum, double x_lower, double x_upper)

   This function sets, or resets, an existing minimizer :data:`s` to use the
   function :data:`f` and the initial search interval [:data:`x_lower`,
   :data:`x_upper`], with a guess for the location of the minimum
   :data:`x_minimum`.

   If the interval given does not contain a minimum, then the function
   returns an error code of :macro:`GSL_EINVAL`.

.. function:: int gsl_min_fminimizer_set_with_values (gsl_min_fminimizer * s, gsl_function * f, double x_minimum, double f_minimum, double x_lower, double f_lower, double x_upper, double f_upper)

   This function is equivalent to :func:`gsl_min_fminimizer_set` but uses
   the values :data:`f_minimum`, :data:`f_lower` and :data:`f_upper` instead of
   computing :code:`f(x_minimum)`, :code:`f(x_lower)` and :code:`f(x_upper)`.

.. function:: void gsl_min_fminimizer_free (gsl_min_fminimizer * s)

   This function frees all the memory associated with the minimizer
   :data:`s`.

.. function:: const char * gsl_min_fminimizer_name (const gsl_min_fminimizer * s)

   This function returns a pointer to the name of the minimizer.  For example::

      printf ("s is a '%s' minimizer\n", gsl_min_fminimizer_name (s));

   would print something like :code:`s is a 'brent' minimizer`.

.. index::
   single: minimization, providing a function to minimize

Providing the function to minimize
==================================

You must provide a continuous function of one variable for the
minimizers to operate on.  In order to allow for general parameters the
functions are defined by a :type:`gsl_function` data type
(:ref:`providing-function-to-solve`).

Iteration
=========

The following functions drive the iteration of each algorithm.  Each
function performs one iteration to update the state of any minimizer of the
corresponding type.  The same functions work for all minimizers so that
different methods can be substituted at runtime without modifications to
the code.

.. function:: int gsl_min_fminimizer_iterate (gsl_min_fminimizer * s)

   This function performs a single iteration of the minimizer :data:`s`.  If the
   iteration encounters an unexpected problem then an error code will be
   returned,

   :macro:`GSL_EBADFUNC`

      the iteration encountered a singular point where the function evaluated
      to :code:`Inf` or :code:`NaN`.

   :macro:`GSL_FAILURE`

      the algorithm could not improve the current best approximation or
      bounding interval.

The minimizer maintains a current best estimate of the position of the
minimum at all times, and the current interval bounding the minimum.
This information can be accessed with the following auxiliary functions,

.. function:: double gsl_min_fminimizer_x_minimum (const gsl_min_fminimizer * s)

   This function returns the current estimate of the position of the
   minimum for the minimizer :data:`s`.

.. function:: double gsl_min_fminimizer_x_upper (const gsl_min_fminimizer * s)
              double gsl_min_fminimizer_x_lower (const gsl_min_fminimizer * s)

   These functions return the current upper and lower bound of the interval
   for the minimizer :data:`s`.

.. function:: double gsl_min_fminimizer_f_minimum (const gsl_min_fminimizer * s)
              double gsl_min_fminimizer_f_upper (const gsl_min_fminimizer * s)
              double gsl_min_fminimizer_f_lower (const gsl_min_fminimizer * s)

   These functions return the value of the function at the current estimate
   of the minimum and at the upper and lower bounds of the interval for the
   minimizer :data:`s`.

.. index::
   single: minimization, stopping parameters

Stopping Parameters
===================

A minimization procedure should stop when one of the following
conditions is true:

* A minimum has been found to within the user-specified precision.
* A user-specified maximum number of iterations has been reached.
* An error has occurred.

The handling of these conditions is under user control.  The function
below allows the user to test the precision of the current result.

.. function:: int gsl_min_test_interval (double x_lower, double x_upper,  double epsabs, double epsrel)

   This function tests for the convergence of the interval [:data:`x_lower`,
   :data:`x_upper`] with absolute error :data:`epsabs` and relative error
   :data:`epsrel`.  The test returns :macro:`GSL_SUCCESS` if the following
   condition is achieved,

   .. only:: not texinfo

      .. math:: |a - b| < \hbox{\it epsabs} + \hbox{\it epsrel\/}\, \min(|a|,|b|)

   .. only:: texinfo

      ::

         |a - b| < epsabs + epsrel min(|a|,|b|) 

   when the interval :math:`x = [a,b]` does not include the origin.  If the
   interval includes the origin then :math:`\min(|a|,|b|)` is replaced by
   zero (which is the minimum value of :math:`|x|` over the interval).  This
   ensures that the relative error is accurately estimated for minima close
   to the origin.

   This condition on the interval also implies that any estimate of the
   minimum :math:`x_m` in the interval satisfies the same condition with respect
   to the true minimum :math:`x_m^*`,

   .. only:: not texinfo

      .. math:: |x_m - x_m^*| < \hbox{\it epsabs} + \hbox{\it epsrel\/}\, x_m^*

   .. only:: texinfo

      ::

         |x_m - x_m^*| < epsabs + epsrel x_m^*

   assuming that the true minimum :math:`x_m^*` is contained within the interval.

Minimization Algorithms
=======================

The minimization algorithms described in this section require an initial
interval which is guaranteed to contain a minimum---if :math:`a` and
:math:`b` are the endpoints of the interval and :math:`x` is an estimate
of the minimum then :math:`f(a) > f(x) < f(b)`.  This ensures that the
function has at least one minimum somewhere in the interval.  If a valid
initial interval is used then these algorithm cannot fail, provided the
function is well-behaved.

.. type:: gsl_min_fminimizer_type

   .. index::
      single: golden section algorithm for finding minima
      single: minimum finding, golden section algorithm

   .. var:: gsl_min_fminimizer_type * gsl_min_fminimizer_goldensection

      The *golden section algorithm* is the simplest method of bracketing
      the minimum of a function.  It is the slowest algorithm provided by the
      library, with linear convergence.

      On each iteration, the algorithm first compares the subintervals from
      the endpoints to the current minimum.  The larger subinterval is divided
      in a golden section (using the famous ratio :math:`(3-\sqrt 5)/2 \approx 0.3819660`
      and the value of the function at this new point is
      calculated.  The new value is used with the constraint :math:`f(a') > f(x') < f(b')`
      to a select new interval containing the minimum, by
      discarding the least useful point.  This procedure can be continued
      indefinitely until the interval is sufficiently small.  Choosing the
      golden section as the bisection ratio can be shown to provide the
      fastest convergence for this type of algorithm.

   .. index::
      single: Brent's method for finding minima
      single: minimum finding, Brent's method

   .. var:: gsl_min_fminimizer_type * gsl_min_fminimizer_brent

      The *Brent minimization algorithm* combines a parabolic
      interpolation with the golden section algorithm.  This produces a fast
      algorithm which is still robust.

      The outline of the algorithm can be summarized as follows: on each
      iteration Brent's method approximates the function using an
      interpolating parabola through three existing points.  The minimum of the
      parabola is taken as a guess for the minimum.  If it lies within the
      bounds of the current interval then the interpolating point is accepted,
      and used to generate a smaller interval.  If the interpolating point is
      not accepted then the algorithm falls back to an ordinary golden section
      step.  The full details of Brent's method include some additional checks
      to improve convergence.

   .. index:: safeguarded step-length algorithm

   .. var:: gsl_min_fminimizer_type * gsl_min_fminimizer_quad_golden

      This is a variant of Brent's algorithm which uses the safeguarded
      step-length algorithm of Gill and Murray.

Examples
========

The following program uses the Brent algorithm to find the minimum of
the function :math:`f(x) = \cos(x) + 1`, which occurs at :math:`x = \pi`.
The starting interval is :math:`(0,6)`, with an initial guess for the
minimum of :math:`2`.

.. include:: examples/min.c
   :code:

Here are the results of the minimization procedure.

.. include:: examples/min.txt
   :code:

References and Further Reading
==============================

Further information on Brent's algorithm is available in the following
book,

* Richard Brent, *Algorithms for minimization without derivatives*,
  Prentice-Hall (1973), republished by Dover in paperback (2002), ISBN
  0-486-41998-3.
