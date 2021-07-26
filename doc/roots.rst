.. index::
   single: root finding
   single: zero finding
   single: finding roots
   single: finding zeros
   single: roots
   single: solving a nonlinear equation
   single: nonlinear equation, solutions of

****************************
One Dimensional Root-Finding
****************************

This chapter describes routines for finding roots of arbitrary
one-dimensional functions.  The library provides low level components
for a variety of iterative solvers and convergence tests.  These can be
combined by the user to achieve the desired solution, with full access
to the intermediate steps of the iteration.  Each class of methods uses
the same framework, so that you can switch between solvers at runtime
without needing to recompile your program.  Each instance of a solver
keeps track of its own state, allowing the solvers to be used in
multi-threaded programs.

The header file :file:`gsl_roots.h` contains prototypes for the root
finding functions and related declarations.

.. index::
   single: root finding, overview

Overview
========

One-dimensional root finding algorithms can be divided into two classes,
*root bracketing* and *root polishing*.  Algorithms which proceed
by bracketing a root are guaranteed to converge.  Bracketing algorithms
begin with a bounded region known to contain a root.  The size of this
bounded region is reduced, iteratively, until it encloses the root to a
desired tolerance.  This provides a rigorous error estimate for the
location of the root.

The technique of *root polishing* attempts to improve an initial
guess to the root.  These algorithms converge only if started "close
enough" to a root, and sacrifice a rigorous error bound for speed.  By
approximating the behavior of a function in the vicinity of a root they
attempt to find a higher order improvement of an initial guess.  When the
behavior of the function is compatible with the algorithm and a good
initial guess is available a polishing algorithm can provide rapid
convergence.

In GSL both types of algorithm are available in similar frameworks.  The
user provides a high-level driver for the algorithms, and the library
provides the individual functions necessary for each of the steps.
There are three main phases of the iteration.  The steps are,

* initialize solver state, :data:`s`, for algorithm :data:`T`

* update :data:`s` using the iteration :data:`T`

* test :data:`s` for convergence, and repeat iteration if necessary

The state for bracketing solvers is held in a :type:`gsl_root_fsolver`
struct.  The updating procedure uses only function evaluations (not
derivatives).  The state for root polishing solvers is held in a
:type:`gsl_root_fdfsolver` struct.  The updates require both the function
and its derivative (hence the name :code:`fdf`) to be supplied by the
user.

.. index::
   single: root finding, caveats

Caveats
=======

Note that root finding functions can only search for one root at a time.
When there are several roots in the search area, the first root to be
found will be returned; however it is difficult to predict which of the
roots this will be. *In most cases, no error will be reported if
you try to find a root in an area where there is more than one.*

Care must be taken when a function may have a multiple root (such as 
:math:`f(x) = (x-x_0)^2` or
:math:`f(x) = (x-x_0)^3`.
It is not possible to use root-bracketing algorithms on
even-multiplicity roots.  For these algorithms the initial interval must
contain a zero-crossing, where the function is negative at one end of
the interval and positive at the other end.  Roots with even-multiplicity
do not cross zero, but only touch it instantaneously.  Algorithms based
on root bracketing will still work for odd-multiplicity roots
(e.g. cubic, quintic, ...). 
Root polishing algorithms generally work with higher multiplicity roots,
but at a reduced rate of convergence.  In these cases the *Steffenson
algorithm* can be used to accelerate the convergence of multiple roots.

While it is not absolutely required that :math:`f` have a root within the
search region, numerical root finding functions should not be used
haphazardly to check for the *existence* of roots.  There are better
ways to do this.  Because it is easy to create situations where numerical
root finders can fail, it is a bad idea to throw a root finder at a
function you do not know much about.  In general it is best to examine
the function visually by plotting before searching for a root.

Initializing the Solver
=======================

.. type:: gsl_root_fsolver

   This is a workspace for finding roots using methods which do not require
   derivatives.

.. type:: gsl_root_fdfsolver

   This is a workspace for finding roots using methods which require
   derivatives.

.. function:: gsl_root_fsolver * gsl_root_fsolver_alloc (const gsl_root_fsolver_type * T)

   This function returns a pointer to a newly allocated instance of a
   solver of type :data:`T`.  For example, the following code creates an
   instance of a bisection solver::

      const gsl_root_fsolver_type * T = gsl_root_fsolver_bisection;
      gsl_root_fsolver * s = gsl_root_fsolver_alloc (T);

   If there is insufficient memory to create the solver then the function
   returns a null pointer and the error handler is invoked with an error
   code of :macro:`GSL_ENOMEM`.

.. function:: gsl_root_fdfsolver * gsl_root_fdfsolver_alloc (const gsl_root_fdfsolver_type * T)

   This function returns a pointer to a newly allocated instance of a
   derivative-based solver of type :data:`T`.  For example, the following
   code creates an instance of a Newton-Raphson solver::

      const gsl_root_fdfsolver_type * T = gsl_root_fdfsolver_newton;
      gsl_root_fdfsolver * s = gsl_root_fdfsolver_alloc (T);

   If there is insufficient memory to create the solver then the function
   returns a null pointer and the error handler is invoked with an error
   code of :macro:`GSL_ENOMEM`.

.. function:: int gsl_root_fsolver_set (gsl_root_fsolver * s, gsl_function * f, double x_lower, double x_upper)

   This function initializes, or reinitializes, an existing solver :data:`s`
   to use the function :data:`f` and the initial search interval
   [:data:`x_lower`, :data:`x_upper`].

.. function:: int gsl_root_fdfsolver_set (gsl_root_fdfsolver * s, gsl_function_fdf * fdf, double root)

   This function initializes, or reinitializes, an existing solver :data:`s`
   to use the function and derivative :data:`fdf` and the initial guess
   :data:`root`.

.. function:: void gsl_root_fsolver_free (gsl_root_fsolver * s)
              void gsl_root_fdfsolver_free (gsl_root_fdfsolver * s)

   These functions free all the memory associated with the solver :data:`s`.

.. function:: const char * gsl_root_fsolver_name (const gsl_root_fsolver * s)
              const char * gsl_root_fdfsolver_name (const gsl_root_fdfsolver * s)

   These functions return a pointer to the name of the solver.  For example::

      printf ("s is a '%s' solver\n", gsl_root_fsolver_name (s));

   would print something like :code:`s is a 'bisection' solver`.

.. index::
   single: root finding, providing a function to solve

.. _providing-function-to-solve:

Providing the function to solve
===============================

You must provide a continuous function of one variable for the root
finders to operate on, and, sometimes, its first derivative.  In order
to allow for general parameters the functions are defined by the
following data types:

.. type:: gsl_function 

   This data type defines a general function with parameters. 

   :code:`double (* function) (double x, void * params)`

      this function should return the value
      :math:`f(x,params)` for argument :data:`x` and parameters :data:`params`

   :code:`void * params`

      a pointer to the parameters of the function

Here is an example for the general quadratic function,

.. math:: f(x) = a x^2 + b x + c

with :math:`a = 3`, :math:`b = 2`, :math:`c = 1`.  The following code
defines a :type:`gsl_function` :code:`F` which you could pass to a root
finder as a function pointer::

  struct my_f_params { double a; double b; double c; };

  double
  my_f (double x, void * p)
    {
      struct my_f_params * params = (struct my_f_params *)p;
      double a = (params->a);
      double b = (params->b);
      double c = (params->c);

      return  (a * x + b) * x + c;
    }

  gsl_function F;
  struct my_f_params params = { 3.0, 2.0, 1.0 };

  F.function = &my_f;
  F.params = &params;

The function :math:`f(x)` can be evaluated using the macro
:code:`GSL_FN_EVAL(&F,x)` defined in :file:`gsl_math.h`.

.. type:: gsl_function_fdf

   This data type defines a general function with parameters and its first
   derivative.

   :code:`double (* f) (double x, void * params)`

      this function should return the value of
      :math:`f(x,params)` for argument :data:`x` and parameters :data:`params`

   :code:`double (* df) (double x, void * params)`

      this function should return the value of the derivative of :data:`f` with
      respect to :data:`x`,
      :math:`f'(x,params)`, for argument :data:`x` and parameters :data:`params`

   :code:`void (* fdf) (double x, void * params, double * f, double * df)`

      this function should set the values of the function :data:`f` to 
      :math:`f(x,params)`
      and its derivative :data:`df` to
      :math:`f'(x,params)`
      for argument :data:`x` and parameters :data:`params`.  This function
      provides an optimization of the separate functions for :math:`f(x)` and
      :math:`f'(x)`---it is always faster to compute the function and its
      derivative at the same time.

   :code:`void * params`

      a pointer to the parameters of the function

Here is an example where 
:math:`f(x) = \exp(2x)`::

  double
  my_f (double x, void * params)
  {
     return exp (2 * x);
  }

  double
  my_df (double x, void * params)
  {
     return 2 * exp (2 * x);
  }

  void
  my_fdf (double x, void * params, 
          double * f, double * df)
  {
     double t = exp (2 * x);

     *f = t;
     *df = 2 * t;   /* uses existing value */
  }

  gsl_function_fdf FDF;

  FDF.f = &my_f;
  FDF.df = &my_df;
  FDF.fdf = &my_fdf;
  FDF.params = 0;

The function :math:`f(x)` can be evaluated using the macro
:code:`GSL_FN_FDF_EVAL_F(&FDF,x)` and the derivative :math:`f'(x)` can
be evaluated using the macro :code:`GSL_FN_FDF_EVAL_DF(&FDF,x)`.  Both
the function :math:`y = f(x)` and its derivative :math:`dy = f'(x)` can
be evaluated at the same time using the macro
:code:`GSL_FN_FDF_EVAL_F_DF(&FDF,x,y,dy)`.  The macro stores
:math:`f(x)` in its :data:`y` argument and :math:`f'(x)` in its :data:`dy`
argument---both of these should be pointers to :code:`double`.

.. index::
   single: root finding, search bounds
   single: root finding, initial guess

Search Bounds and Guesses
=========================

You provide either search bounds or an initial guess; this section
explains how search bounds and guesses work and how function arguments
control them.

A guess is simply an :math:`x` value which is iterated until it is within
the desired precision of a root.  It takes the form of a :code:`double`.

Search bounds are the endpoints of an interval which is iterated until
the length of the interval is smaller than the requested precision.  The
interval is defined by two values, the lower limit and the upper limit.
Whether the endpoints are intended to be included in the interval or not
depends on the context in which the interval is used.

Iteration
=========

The following functions drive the iteration of each algorithm.  Each
function performs one iteration to update the state of any solver of the
corresponding type.  The same functions work for all solvers so that
different methods can be substituted at runtime without modifications to
the code.

.. function:: int gsl_root_fsolver_iterate (gsl_root_fsolver * s)
              int gsl_root_fdfsolver_iterate (gsl_root_fdfsolver * s)

   These functions perform a single iteration of the solver :data:`s`.  If the
   iteration encounters an unexpected problem then an error code will be
   returned,

   :code:`GSL_EBADFUNC`

      the iteration encountered a singular point where the function or its
      derivative evaluated to :code:`Inf` or :code:`NaN`.

   :code:`GSL_EZERODIV`

      the derivative of the function vanished at the iteration point,
      preventing the algorithm from continuing without a division by zero.

The solver maintains a current best estimate of the root at all
times.  The bracketing solvers also keep track of the current best
interval bounding the root.  This information can be accessed with the
following auxiliary functions,

.. function:: double gsl_root_fsolver_root (const gsl_root_fsolver * s)
              double gsl_root_fdfsolver_root (const gsl_root_fdfsolver * s)

   These functions return the current estimate of the root for the solver :data:`s`.

.. function:: double gsl_root_fsolver_x_lower (const gsl_root_fsolver * s)
              double gsl_root_fsolver_x_upper (const gsl_root_fsolver * s)

   These functions return the current bracketing interval for the solver :data:`s`.

.. index::
   single: root finding, stopping parameters

Search Stopping Parameters
==========================

A root finding procedure should stop when one of the following conditions is
true:

* A root has been found to within the user-specified precision.
* A user-specified maximum number of iterations has been reached.
* An error has occurred.

The handling of these conditions is under user control.  The functions
below allow the user to test the precision of the current result in
several standard ways.

.. function:: int gsl_root_test_interval (double x_lower, double x_upper, double epsabs, double epsrel)

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
   ensures that the relative error is accurately estimated for roots close
   to the origin.

   This condition on the interval also implies that any estimate of the
   root :math:`r` in the interval satisfies the same condition with respect
   to the true root :math:`r^*`,

   .. only:: not texinfo

      .. math:: |r - r^*| < \hbox{\it epsabs} + \hbox{\it epsrel\/}\, r^*

   .. only:: texinfo

      ::

         |r - r^*| < epsabs + epsrel r^*

   assuming that the true root :math:`r^*` is contained within the interval.

.. function:: int gsl_root_test_delta (double x1, double x0, double epsabs, double epsrel)

   This function tests for the convergence of the sequence :data:`x0`,
   :data:`x1` with absolute error :data:`epsabs` and relative error
   :data:`epsrel`.  The test returns :macro:`GSL_SUCCESS` if the following
   condition is achieved,

   .. only:: not texinfo

      .. math:: |x_1 - x_0| < \hbox{\it epsabs} + \hbox{\it epsrel\/}\, |x_1|

   .. only:: texinfo

      ::

         |x_1 - x_0| < epsabs + epsrel |x_1|

   and returns :macro:`GSL_CONTINUE` otherwise.

.. function:: int gsl_root_test_residual (double f, double epsabs)

   This function tests the residual value :data:`f` against the absolute
   error bound :data:`epsabs`.  The test returns :macro:`GSL_SUCCESS` if the
   following condition is achieved,

   .. only:: not texinfo

      .. math:: |f| < \hbox{\it epsabs}

   .. only:: texinfo

      ::

         |f| < epsabs

   and returns :macro:`GSL_CONTINUE` otherwise.  This criterion is suitable
   for situations where the precise location of the root, :math:`x`, is
   unimportant provided a value can be found where the residual,
   :math:`|f(x)|`, is small enough.

Root Bracketing Algorithms
==========================

The root bracketing algorithms described in this section require an
initial interval which is guaranteed to contain a root---if :math:`a`
and :math:`b` are the endpoints of the interval then :math:`f(a)` must
differ in sign from :math:`f(b)`.  This ensures that the function crosses
zero at least once in the interval.  If a valid initial interval is used
then these algorithm cannot fail, provided the function is well-behaved.

Note that a bracketing algorithm cannot find roots of even degree, since
these do not cross the :math:`x`-axis.

.. type:: gsl_root_fsolver_type

   .. index::
      single: bisection algorithm for finding roots
      single: root finding, bisection algorithm

   .. var:: gsl_root_fsolver_type * gsl_root_fsolver_bisection

      The *bisection algorithm* is the simplest method of bracketing the
      roots of a function.   It is the slowest algorithm provided by
      the library, with linear convergence.

      On each iteration, the interval is bisected and the value of the
      function at the midpoint is calculated.  The sign of this value is used
      to determine which half of the interval does not contain a root.  That
      half is discarded to give a new, smaller interval containing the
      root.  This procedure can be continued indefinitely until the interval is
      sufficiently small.

      At any time the current estimate of the root is taken as the midpoint of
      the interval.

      .. eps file "roots-bisection.eps"
      .. @iftex
      .. @sp 1
      .. @center @image{roots-bisection,3.4in}

      .. @quotation
      .. Four iterations of bisection, where :math:`a_n` is :math:`n`-th position of
      .. the beginning of the interval and :math:`b_n` is the :math:`n`-th position
      .. of the end.  The midpoint of each interval is also indicated.
      .. @end quotation
      .. @end iftex

   .. index::
      single: false position algorithm for finding roots
      single: root finding, false position algorithm

   .. var:: gsl_root_fsolver_type * gsl_root_fsolver_falsepos

      The *false position algorithm* is a method of finding roots based on
      linear interpolation.  Its convergence is linear, but it is usually
      faster than bisection.

      On each iteration a line is drawn between the endpoints :math:`(a,f(a))`
      and :math:`(b,f(b))` and the point where this line crosses the
      :math:`x`-axis taken as a "midpoint".  The value of the function at
      this point is calculated and its sign is used to determine which side of
      the interval does not contain a root.  That side is discarded to give a
      new, smaller interval containing the root.  This procedure can be
      continued indefinitely until the interval is sufficiently small.

      The best estimate of the root is taken from the linear interpolation of
      the interval on the current iteration.

      .. eps file "roots-false-position.eps"
      .. @iftex
      .. @image{roots-false-position,4in}
      .. @quotation
      .. Several iterations of false position, where :math:`a_n` is :math:`n`-th
      .. position of the beginning of the interval and :math:`b_n` is the
      .. :math:`n`-th position of the end.
      .. @end quotation
      .. @end iftex

   .. index::
      single: Brent's method for finding roots
      single: root finding, Brent's method

   .. var:: gsl_root_fsolver_type * gsl_root_fsolver_brent

      The *Brent-Dekker method* (referred to here as *Brent's method*)
      combines an interpolation strategy with the bisection algorithm.  This
      produces a fast algorithm which is still robust.

      On each iteration Brent's method approximates the function using an
      interpolating curve.  On the first iteration this is a linear
      interpolation of the two endpoints.  For subsequent iterations the
      algorithm uses an inverse quadratic fit to the last three points, for
      higher accuracy.  The intercept of the interpolating curve with the
      :math:`x`-axis is taken as a guess for the root.  If it lies within the
      bounds of the current interval then the interpolating point is accepted,
      and used to generate a smaller interval.  If the interpolating point is
      not accepted then the algorithm falls back to an ordinary bisection
      step.

      The best estimate of the root is taken from the most recent
      interpolation or bisection.

Root Finding Algorithms using Derivatives
=========================================

The root polishing algorithms described in this section require an
initial guess for the location of the root.  There is no absolute
guarantee of convergence---the function must be suitable for this
technique and the initial guess must be sufficiently close to the root
for it to work.  When these conditions are satisfied then convergence is
quadratic.

These algorithms make use of both the function and its derivative. 

.. type:: gsl_root_fdfsolver_type

   .. index::
      single: Newton's method for finding roots
      single: root finding, Newton's method

   .. var:: gsl_root_fdfsolver_type * gsl_root_fdfsolver_newton

      Newton's Method is the standard root-polishing algorithm.  The algorithm
      begins with an initial guess for the location of the root.  On each
      iteration, a line tangent to the function :math:`f` is drawn at that
      position.  The point where this line crosses the :math:`x`-axis becomes
      the new guess.  The iteration is defined by the following sequence,

      .. only:: not texinfo

         .. math:: x_{i+1} = x_i - {f(x_i) \over f'(x_i)}

      .. only:: texinfo

         ::

            x_{i+1} = x_i - f(x_i)/f'(x_i)

      Newton's method converges quadratically for single roots, and linearly
      for multiple roots.

      .. eps file "roots-newtons-method.eps"
      .. @iftex
      .. @sp 1
      .. @center @image{roots-newtons-method,3.4in}

      .. @quotation
      .. Several iterations of Newton's Method, where :math:`g_n` is the
      .. :math:`n`-th guess.
      .. @end quotation
      .. @end iftex

   .. index::
      single: secant method for finding roots
      single: root finding, secant method

   .. var:: gsl_root_fdfsolver_type * gsl_root_fdfsolver_secant

      The *secant method* is a simplified version of Newton's method which does
      not require the computation of the derivative on every step.

      On its first iteration the algorithm begins with Newton's method, using
      the derivative to compute a first step,

      .. only:: not texinfo

         .. math:: x_1 = x_0 - {f(x_0) \over f'(x_0)}

      .. only:: texinfo

         ::

            x_1 = x_0 - f(x_0)/f'(x_0)

      Subsequent iterations avoid the evaluation of the derivative by
      replacing it with a numerical estimate, the slope of the line through
      the previous two points,

      .. only:: not texinfo

         .. math::

            x_{i+1} = x_i - {f(x_i) \over f'_{est}}
             ~\hbox{where}~
             f'_{est} =  {f(x_{i}) - f(x_{i-1}) \over x_i - x_{i-1}}

      .. only:: texinfo
      
         ::

            x_{i+1} = x_i f(x_i) / f'_{est} where
             f'_{est} = (f(x_i) - f(x_{i-1})/(x_i - x_{i-1})

      When the derivative does not change significantly in the vicinity of the
      root the secant method gives a useful saving.  Asymptotically the secant
      method is faster than Newton's method whenever the cost of evaluating
      the derivative is more than 0.44 times the cost of evaluating the
      function itself.  As with all methods of computing a numerical
      derivative the estimate can suffer from cancellation errors if the
      separation of the points becomes too small.

      On single roots, the method has a convergence of order :math:`(1 + \sqrt
      5)/2` (approximately :math:`1.62`).  It converges linearly for multiple
      roots.  

      .. eps file "roots-secant-method.eps"
      .. @iftex
      .. @tex
      .. \input epsf
      .. \medskip
      .. \centerline{\epsfxsize=5in\epsfbox{roots-secant-method.eps}}
      .. @end tex
      .. @quotation
      .. Several iterations of Secant Method, where :math:`g_n` is the :math:`n`-th
      .. guess.
      .. @end quotation
      .. @end iftex

   .. index::
      single: Steffenson's method for finding roots
      single: root finding, Steffenson's method

   .. var:: gsl_root_fdfsolver_type * gsl_root_fdfsolver_steffenson

      The *Steffenson Method* [#f1]_
      provides the fastest
      convergence of all the routines.  It combines the basic Newton
      algorithm with an Aitken "delta-squared" acceleration.  If the
      Newton iterates are :math:`x_i` then the acceleration procedure
      generates a new sequence :math:`R_i`,

      .. only:: not texinfo

         .. math:: R_i = x_i - {(x_{i+1} - x_i)^2 \over (x_{i+2} - 2 x_{i+1} + x_i)}

      .. only:: texinfo

         ::

            R_i = x_i - (x_{i+1} - x_i)^2 / (x_{i+2} - 2 x_{i+1} + x_{i})

      which converges faster than the original sequence under reasonable
      conditions.  The new sequence requires three terms before it can produce
      its first value so the method returns accelerated values on the second
      and subsequent iterations.  On the first iteration it returns the
      ordinary Newton estimate.  The Newton iterate is also returned if the
      denominator of the acceleration term ever becomes zero.

      As with all acceleration procedures this method can become unstable if
      the function is not well-behaved. 

Examples
========

For any root finding algorithm we need to prepare the function to be
solved.  For this example we will use the general quadratic equation
described earlier.  We first need a header file (:file:`demo_fn.h`) to
define the function parameters,

.. include:: examples/demo_fn.h
   :code:

We place the function definitions in a separate file (:file:`demo_fn.c`),

.. include:: examples/demo_fn.c
   :code:

The first program uses the function solver :data:`gsl_root_fsolver_brent`
for Brent's method and the general quadratic defined above to solve the
following equation,

.. math:: x^2 - 5 = 0

with solution :math:`x = \sqrt 5 = 2.236068...`

.. include:: examples/roots.c
   :code:

Here are the results of the iterations::

  $ ./a.out 
  using brent method
   iter [    lower,     upper]      root        err  err(est)
      1 [1.0000000, 5.0000000] 1.0000000 -1.2360680 4.0000000
      2 [1.0000000, 3.0000000] 3.0000000 +0.7639320 2.0000000
      3 [2.0000000, 3.0000000] 2.0000000 -0.2360680 1.0000000
      4 [2.2000000, 3.0000000] 2.2000000 -0.0360680 0.8000000
      5 [2.2000000, 2.2366300] 2.2366300 +0.0005621 0.0366300
  Converged:                            
      6 [2.2360634, 2.2366300] 2.2360634 -0.0000046 0.0005666

If the program is modified to use the bisection solver instead of
Brent's method, by changing :data:`gsl_root_fsolver_brent` to
:data:`gsl_root_fsolver_bisection` the slower convergence of the
Bisection method can be observed::

  $ ./a.out 
  using bisection method
   iter [    lower,     upper]      root        err  err(est)
      1 [0.0000000, 2.5000000] 1.2500000 -0.9860680 2.5000000
      2 [1.2500000, 2.5000000] 1.8750000 -0.3610680 1.2500000
      3 [1.8750000, 2.5000000] 2.1875000 -0.0485680 0.6250000
      4 [2.1875000, 2.5000000] 2.3437500 +0.1076820 0.3125000
      5 [2.1875000, 2.3437500] 2.2656250 +0.0295570 0.1562500
      6 [2.1875000, 2.2656250] 2.2265625 -0.0095055 0.0781250
      7 [2.2265625, 2.2656250] 2.2460938 +0.0100258 0.0390625
      8 [2.2265625, 2.2460938] 2.2363281 +0.0002601 0.0195312
      9 [2.2265625, 2.2363281] 2.2314453 -0.0046227 0.0097656
     10 [2.2314453, 2.2363281] 2.2338867 -0.0021813 0.0048828
     11 [2.2338867, 2.2363281] 2.2351074 -0.0009606 0.0024414
  Converged:                            
     12 [2.2351074, 2.2363281] 2.2357178 -0.0003502 0.0012207

The next program solves the same function using a derivative solver
instead.

.. include:: examples/rootnewt.c
   :code:

Here are the results for Newton's method::

  $ ./a.out 
  using newton method
  iter        root        err   err(est)
      1  3.0000000 +0.7639320 -2.0000000
      2  2.3333333 +0.0972654 -0.6666667
      3  2.2380952 +0.0020273 -0.0952381
  Converged:      
      4  2.2360689 +0.0000009 -0.0020263

Note that the error can be estimated more accurately by taking the
difference between the current iterate and next iterate rather than the
previous iterate.  The other derivative solvers can be investigated by
changing :data:`gsl_root_fdfsolver_newton` to
:data:`gsl_root_fdfsolver_secant` or
:data:`gsl_root_fdfsolver_steffenson`.

References and Further Reading
==============================

For information on the Brent-Dekker algorithm see the following two
papers,

* R. P. Brent, "An algorithm with guaranteed convergence for finding a
  zero of a function", *Computer Journal*, 14 (1971) 422--425

* J. C. P. Bus and T. J. Dekker, "Two Efficient Algorithms with Guaranteed
  Convergence for Finding a Zero of a Function", *ACM Transactions of
  Mathematical Software*, Vol.: 1 No.: 4 (1975) 330--345

.. rubric:: Footnotes

.. [#f1] J.F. Steffensen (1873--1961). The spelling used in the name of the
         function is slightly incorrect, but has been preserved to avoid incompatibility.
