.. index::
   single: quadrature
   single: numerical integration (quadrature)
   single: integration, numerical (quadrature)
   single: QUADPACK

*********************
Numerical Integration
*********************

.. include:: include.rst

This chapter describes routines for performing numerical integration
(quadrature) of a function in one dimension.  There are routines for
adaptive and non-adaptive integration of general functions, with
specialised routines for specific cases.  These include integration over
infinite and semi-infinite ranges, singular integrals, including
logarithmic singularities, computation of Cauchy principal values and
oscillatory integrals.  The library reimplements the algorithms used in
|quadpack|, a numerical integration package written by Piessens,
de Doncker-Kapenga, Ueberhuber and Kahaner.  Fortran code for |quadpack| is
available on Netlib.  Also included are non-adaptive, fixed-order
Gauss-Legendre integration routines with high precision coefficients, as
well as fixed-order quadrature rules for a variety of weighting functions
from IQPACK.

The functions described in this chapter are declared in the header file
:file:`gsl_integration.h`.

Introduction
============

Each algorithm computes an approximation to a definite integral of the
form,

.. math:: I = \int_a^b f(x) w(x) dx

where :math:`w(x)` is a weight function (for general integrands :math:`w(x) = 1`).
The user provides absolute and relative error bounds 
:math:`(epsabs, epsrel)` which specify the following accuracy requirement,

.. only:: not texinfo

   .. math:: |RESULT - I| \leq \max{(epsabs, epsrel |I|)}

.. only:: texinfo

   .. math:: |RESULT - I| <= max(epsabs, epsrel |I|)

where :math:`RESULT` is the numerical approximation obtained by the
algorithm.  The algorithms attempt to estimate the absolute error
:math:`ABSERR = |RESULT - I|` in such a way that the following inequality
holds,

.. only:: not texinfo

   .. math:: |RESULT - I| \leq ABSERR \leq \max{(epsabs, epsrel |I|)}

.. only:: texinfo

   .. math:: |RESULT - I| <= ABSERR <= max(epsabs, epsrel |I|)

In short, the routines return the first approximation 
which has an absolute error smaller than
:math:`epsabs` or a relative error smaller than :math:`epsrel`.   

Note that this is an *either-or* constraint, 
not simultaneous.  To compute to a specified absolute error, set
:math:`epsrel` to zero.  To compute to a specified relative error,
set :math:`epsabs` to zero.  
The routines will fail to converge if the error bounds are too
stringent, but always return the best approximation obtained up to
that stage.

The algorithms in |quadpack| use a naming convention based on the
following letters::

  Q - quadrature routine

  N - non-adaptive integrator
  A - adaptive integrator

  G - general integrand (user-defined)
  W - weight function with integrand

  S - singularities can be more readily integrated
  P - points of special difficulty can be supplied
  I - infinite range of integration
  O - oscillatory weight function, cos or sin
  F - Fourier integral
  C - Cauchy principal value

The algorithms are built on pairs of quadrature rules, a higher order
rule and a lower order rule.  The higher order rule is used to compute
the best approximation to an integral over a small range.  The
difference between the results of the higher order rule and the lower
order rule gives an estimate of the error in the approximation.

.. index:: Gauss-Kronrod quadrature

Integrands without weight functions
-----------------------------------

The algorithms for general functions (without a weight function) are
based on Gauss-Kronrod rules. 

A Gauss-Kronrod rule begins with a classical Gaussian quadrature rule of
order :math:`m`.  This is extended with additional points between each of
the abscissae to give a higher order Kronrod rule of order :math:`2m + 1`.
The Kronrod rule is efficient because it reuses existing function
evaluations from the Gaussian rule.  

The higher order Kronrod rule is used as the best approximation to the
integral, and the difference between the two rules is used as an
estimate of the error in the approximation.

Integrands with weight functions
--------------------------------

.. index::
   single: Clenshaw-Curtis quadrature
   single: Modified Clenshaw-Curtis quadrature

For integrands with weight functions the algorithms use Clenshaw-Curtis
quadrature rules.  

A Clenshaw-Curtis rule begins with an :math:`n`-th order Chebyshev
polynomial approximation to the integrand.  This polynomial can be
integrated exactly to give an approximation to the integral of the
original function.  The Chebyshev expansion can be extended to higher
orders to improve the approximation and provide an estimate of the
error.

Integrands with singular weight functions
-----------------------------------------

The presence of singularities (or other behavior) in the integrand can
cause slow convergence in the Chebyshev approximation.  The modified
Clenshaw-Curtis rules used in |quadpack| separate out several common
weight functions which cause slow convergence.  

These weight functions are integrated analytically against the Chebyshev
polynomials to precompute *modified Chebyshev moments*.  Combining
the moments with the Chebyshev approximation to the function gives the
desired integral.  The use of analytic integration for the singular part
of the function allows exact cancellations and substantially improves
the overall convergence behavior of the integration.

QNG non-adaptive Gauss-Kronrod integration
==========================================
.. index:: QNG quadrature algorithm

The QNG algorithm is a non-adaptive procedure which uses fixed
Gauss-Kronrod-Patterson abscissae to sample the integrand at a maximum of 87
points.  It is provided for fast integration of smooth functions.

.. function:: int gsl_integration_qng (const gsl_function * f, double a, double b, double epsabs, double epsrel, double * result, double * abserr, size_t * neval)

   This function applies the Gauss-Kronrod 10-point, 21-point, 43-point and
   87-point integration rules in succession until an estimate of the
   integral of :math:`f` over :math:`(a,b)` is achieved within the desired
   absolute and relative error limits, :data:`epsabs` and :data:`epsrel`.  The
   function returns the final approximation, :data:`result`, an estimate of
   the absolute error, :data:`abserr` and the number of function evaluations
   used, :data:`neval`.  The Gauss-Kronrod rules are designed in such a way
   that each rule uses all the results of its predecessors, in order to
   minimize the total number of function evaluations.


QAG adaptive integration
========================
.. index:: QAG quadrature algorithm

The QAG algorithm is a simple adaptive integration procedure.  The
integration region is divided into subintervals, and on each iteration
the subinterval with the largest estimated error is bisected.  This
reduces the overall error rapidly, as the subintervals become
concentrated around local difficulties in the integrand.  These
subintervals are managed by the following struct,

.. type:: gsl_integration_workspace

   This workspace handles the memory for the subinterval ranges, results and error
   estimates.

.. index:: gsl_integration_workspace

.. function:: gsl_integration_workspace * gsl_integration_workspace_alloc (size_t n) 

   This function allocates a workspace sufficient to hold :data:`n` double
   precision intervals, their integration results and error estimates.
   One workspace may be used multiple times as all necessary reinitialization
   is performed automatically by the integration routines.

.. function:: void gsl_integration_workspace_free (gsl_integration_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_integration_qag (const gsl_function * f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace * workspace,  double * result, double * abserr)

   This function applies an integration rule adaptively until an estimate
   of the integral of :math:`f` over :math:`(a,b)` is achieved within the
   desired absolute and relative error limits, :data:`epsabs` and
   :data:`epsrel`.  The function returns the final approximation,
   :data:`result`, and an estimate of the absolute error, :data:`abserr`.  The
   integration rule is determined by the value of :data:`key`, which should
   be chosen from the following symbolic names,

   ========================== ===
   Symbolic Name              Key
   ========================== ===
   :macro:`GSL_INTEG_GAUSS15` 1
   :macro:`GSL_INTEG_GAUSS21` 2
   :macro:`GSL_INTEG_GAUSS31` 3
   :macro:`GSL_INTEG_GAUSS41` 4
   :macro:`GSL_INTEG_GAUSS51` 5
   :macro:`GSL_INTEG_GAUSS61` 6
   ========================== ===

   corresponding to the 15, 21, 31, 41, 51 and 61 point Gauss-Kronrod
   rules.  The higher-order rules give better accuracy for smooth functions,
   while lower-order rules save time when the function contains local
   difficulties, such as discontinuities.

   On each iteration the adaptive integration strategy bisects the interval
   with the largest error estimate.  The subintervals and their results are
   stored in the memory provided by :data:`workspace`.  The maximum number of
   subintervals is given by :data:`limit`, which may not exceed the allocated
   size of the workspace.

QAGS adaptive integration with singularities
============================================
.. index:: QAGS quadrature algorithm

The presence of an integrable singularity in the integration region
causes an adaptive routine to concentrate new subintervals around the
singularity.  As the subintervals decrease in size the successive
approximations to the integral converge in a limiting fashion.  This
approach to the limit can be accelerated using an extrapolation
procedure.  The QAGS algorithm combines adaptive bisection with the Wynn
epsilon-algorithm to speed up the integration of many types of
integrable singularities.

.. function:: int gsl_integration_qags (const gsl_function * f, double a, double b, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

   This function applies the Gauss-Kronrod 21-point integration rule
   adaptively until an estimate of the integral of :math:`f` over
   :math:`(a,b)` is achieved within the desired absolute and relative error
   limits, :data:`epsabs` and :data:`epsrel`.  The results are extrapolated
   using the epsilon-algorithm, which accelerates the convergence of the
   integral in the presence of discontinuities and integrable
   singularities.  The function returns the final approximation from the
   extrapolation, :data:`result`, and an estimate of the absolute error,
   :data:`abserr`.  The subintervals and their results are stored in the
   memory provided by :data:`workspace`.  The maximum number of subintervals
   is given by :data:`limit`, which may not exceed the allocated size of the
   workspace.

QAGP adaptive integration with known singular points
====================================================
.. index::
   single: QAGP quadrature algorithm
   single: singular points, specifying positions in quadrature

.. function:: int gsl_integration_qagp (const gsl_function * f, double * pts, size_t npts, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

   This function applies the adaptive integration algorithm QAGS taking
   account of the user-supplied locations of singular points.  The array
   :data:`pts` of length :data:`npts` should contain the endpoints of the
   integration ranges defined by the integration region and locations of
   the singularities.  For example, to integrate over the region
   :math:`(a,b)` with break-points at :math:`x_1, x_2, x_3` (where 
   :math:`a < x_1 < x_2 < x_3 < b`) the following :data:`pts` array should be used::

     pts[0] = a
     pts[1] = x_1
     pts[2] = x_2
     pts[3] = x_3
     pts[4] = b

   with :data:`npts` = 5.

   If you know the locations of the singular points in the integration
   region then this routine will be faster than :func:`gsl_integration_qags`.

QAGI adaptive integration on infinite intervals
===============================================
.. index:: QAGI quadrature algorithm

.. function:: int gsl_integration_qagi (gsl_function * f, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

   This function computes the integral of the function :data:`f` over the
   infinite interval :math:`(-\infty,+\infty)`.  The integral is mapped onto the
   semi-open interval :math:`(0,1]` using the transformation :math:`x = (1-t)/t`,

   .. math:: \int_{-\infty}^{+\infty} dx f(x) = \int_0^1 dt (f((1-t)/t) + f(-(1-t)/t))/t^2.

   It is then integrated using the QAGS algorithm.  The normal 21-point
   Gauss-Kronrod rule of QAGS is replaced by a 15-point rule, because the
   transformation can generate an integrable singularity at the origin.  In
   this case a lower-order rule is more efficient.

.. function:: int gsl_integration_qagiu (gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

   This function computes the integral of the function :data:`f` over the
   semi-infinite interval :math:`(a,+\infty)`.  The integral is mapped onto the
   semi-open interval :math:`(0,1]` using the transformation :math:`x = a + (1-t)/t`,

   .. math:: \int_{a}^{+\infty} dx f(x) = \int_0^1 dt f(a + (1-t)/t)/t^2

   and then integrated using the QAGS algorithm.

.. function:: int gsl_integration_qagil (gsl_function * f, double b, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

   This function computes the integral of the function :data:`f` over the
   semi-infinite interval :math:`(-\infty,b)`.  The integral is mapped onto the
   semi-open interval :math:`(0,1]` using the transformation :math:`x = b - (1-t)/t`,

   .. math:: \int_{-\infty}^{b} dx f(x) = \int_0^1 dt f(b - (1-t)/t)/t^2

   and then integrated using the QAGS algorithm.

QAWC adaptive integration for Cauchy principal values
=====================================================
.. index::
   single: QAWC quadrature algorithm
   single: Cauchy principal value, by numerical quadrature

.. function:: int gsl_integration_qawc (gsl_function * f, double a, double b, double c, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

   This function computes the Cauchy principal value of the integral of
   :math:`f` over :math:`(a,b)`, with a singularity at :data:`c`,

   .. only:: not texinfo

      .. math::

         I = \int_a^b dx\, {f(x) \over x - c}
           = \lim_{\epsilon \to 0} 
         \left\{
         \int_a^{c-\epsilon} dx\, {f(x) \over x - c}
         +
         \int_{c+\epsilon}^b dx\, {f(x) \over x - c}
         \right\}

   .. only:: texinfo

      .. math:: I = \int_a^b dx f(x) / (x - c)

   The adaptive bisection algorithm of QAG is used, with modifications to
   ensure that subdivisions do not occur at the singular point :math:`x = c`.
   When a subinterval contains the point :math:`x = c` or is close to
   it then a special 25-point modified Clenshaw-Curtis rule is used to control
   the singularity.  Further away from the
   singularity the algorithm uses an ordinary 15-point Gauss-Kronrod
   integration rule.

QAWS adaptive integration for singular functions
================================================
.. index::
   single: QAWS quadrature algorithm
   single: singular functions, numerical integration of

The QAWS algorithm is designed for integrands with algebraic-logarithmic
singularities at the end-points of an integration region.  In order to
work efficiently the algorithm requires a precomputed table of
Chebyshev moments.

.. type:: gsl_integration_qaws_table

   This structure contains precomputed quantities for the QAWS algorithm.

.. function:: gsl_integration_qaws_table * gsl_integration_qaws_table_alloc (double alpha, double beta, int mu, int nu)

   This function allocates space for a :type:`gsl_integration_qaws_table`
   struct describing a singular weight function
   :math:`w(x)` with the parameters :math:`(\alpha, \beta, \mu, \nu)`,

   .. math:: w(x) = (x - a)^\alpha (b - x)^\beta \log^\mu (x - a) \log^\nu (b - x)

   where :math:`\alpha > -1`, :math:`\beta > -1`, and :math:`\mu = 0, 1`,
   :math:`\nu = 0, 1`.  The weight function can take four different forms
   depending on the values of :math:`\mu` and :math:`\nu`,

   ============================================================      =================
   Weight function :math:`w(x)`                                      :math:`(\mu,\nu)`
   ============================================================      =================
   :math:`(x - a)^\alpha (b - x)^\beta`                              :math:`(0,0)`
   :math:`(x - a)^\alpha (b - x)^\beta \log{(x-a)}`                  :math:`(1,0)`
   :math:`(x - a)^\alpha (b - x)^\beta \log{(b-x)}`                  :math:`(0,1)`
   :math:`(x - a)^\alpha (b - x)^\beta \log{(x-a)} \log{(b-x)}`      :math:`(1,1)`
   ============================================================      =================

   The singular points :math:`(a,b)` do not have to be specified until the
   integral is computed, where they are the endpoints of the integration
   range.

   The function returns a pointer to the newly allocated table
   :type:`gsl_integration_qaws_table` if no errors were detected, and 0 in
   the case of error.

.. function:: int gsl_integration_qaws_table_set (gsl_integration_qaws_table * t, double alpha, double beta, int mu, int nu)

   This function modifies the parameters :math:`(\alpha, \beta, \mu, \nu)` of
   an existing :type:`gsl_integration_qaws_table` struct :data:`t`.

.. function:: void gsl_integration_qaws_table_free (gsl_integration_qaws_table * t)

   This function frees all the memory associated with the
   :type:`gsl_integration_qaws_table` struct :data:`t`.

.. function:: int gsl_integration_qaws (gsl_function * f, const double a, const double b, gsl_integration_qaws_table * t, const double epsabs, const double epsrel, const size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

   This function computes the integral of the function :math:`f(x)` over the
   interval :math:`(a,b)` with the singular weight function
   :math:`(x-a)^\alpha (b-x)^\beta \log^\mu (x-a) \log^\nu (b-x)`.  The parameters 
   of the weight function :math:`(\alpha, \beta, \mu, \nu)` are taken from the
   table :data:`t`.  The integral is,

   .. math:: I = \int_a^b dx f(x) (x - a)^\alpha (b - x)^\beta \log^\mu (x - a) \log^\nu (b - x).

   The adaptive bisection algorithm of QAG is used.  When a subinterval
   contains one of the endpoints then a special 25-point modified
   Clenshaw-Curtis rule is used to control the singularities.  For
   subintervals which do not include the endpoints an ordinary 15-point
   Gauss-Kronrod integration rule is used.

QAWO adaptive integration for oscillatory functions
===================================================
.. index::
   single: QAWO quadrature algorithm
   single: oscillatory functions, numerical integration of

The QAWO algorithm is designed for integrands with an oscillatory
factor, :math:`\sin(\omega x)` or :math:`\cos(\omega x)`.  In order to
work efficiently the algorithm requires a table of Chebyshev moments
which must be pre-computed with calls to the functions below.

.. index:: gsl_integration_qawo_table

.. function:: gsl_integration_qawo_table * gsl_integration_qawo_table_alloc (double omega, double L, enum gsl_integration_qawo_enum sine, size_t n)

   This function allocates space for a :type:`gsl_integration_qawo_table`
   struct and its associated workspace describing a sine or cosine weight
   function :math:`w(x)` with the parameters :math:`(\omega, L)`,

   .. only:: not texinfo

      .. math::

         w(x) =
         \left\{
         \begin{array}{c}
         \sin{(\omega x)} \\
         \cos{(\omega x)} \\
         \end{array}
         \right\}

   .. only:: texinfo

      | w(x) = sin(\omega x)
      | w(x) = cos(\omega x)

   The parameter :data:`L` must be the length of the interval over which the
   function will be integrated :math:`L = b - a`.  The choice of sine or
   cosine is made with the parameter :data:`sine` which should be chosen from
   one of the two following symbolic values:

   .. macro:: GSL_INTEG_COSINE

   .. macro:: GSL_INTEG_SINE

   The :type:`gsl_integration_qawo_table` is a table of the trigonometric
   coefficients required in the integration process.  The parameter :data:`n`
   determines the number of levels of coefficients that are computed.  Each
   level corresponds to one bisection of the interval :math:`L`, so that
   :data:`n` levels are sufficient for subintervals down to the length
   :math:`L/2^n`.  The integration routine :func:`gsl_integration_qawo`
   returns the error :code:`GSL_ETABLE` if the number of levels is
   insufficient for the requested accuracy.

.. function:: int gsl_integration_qawo_table_set (gsl_integration_qawo_table * t, double omega, double L, enum gsl_integration_qawo_enum sine)

   This function changes the parameters :data:`omega`, :data:`L` and :data:`sine`
   of the existing workspace :data:`t`.

.. function:: int gsl_integration_qawo_table_set_length (gsl_integration_qawo_table * t, double L)

   This function allows the length parameter :data:`L` of the workspace
   :data:`t` to be changed.

.. function:: void gsl_integration_qawo_table_free (gsl_integration_qawo_table * t)

   This function frees all the memory associated with the workspace :data:`t`.

.. function:: int gsl_integration_qawo (gsl_function * f, const double a, const double epsabs, const double epsrel, const size_t limit, gsl_integration_workspace * workspace, gsl_integration_qawo_table * wf, double * result, double * abserr)

   This function uses an adaptive algorithm to compute the integral of
   :math:`f` over :math:`(a,b)` with the weight function 
   :math:`\sin(\omega x)` or :math:`\cos(\omega x)` defined 
   by the table :data:`wf`,
 
   .. only:: not texinfo

      .. math::

         I = \int_a^b dx f(x)
         \left\{
         \begin{array}{c}
         \sin{(\omega x)} \\
         \cos{(\omega x)}
         \end{array}
         \right\}

   .. only:: texinfo

      | I = \int_a^b dx f(x) sin(omega x)
      | I = \int_a^b dx f(x) cos(omega x)

   The results are extrapolated using the epsilon-algorithm to accelerate
   the convergence of the integral.  The function returns the final
   approximation from the extrapolation, :data:`result`, and an estimate of
   the absolute error, :data:`abserr`.  The subintervals and their results are
   stored in the memory provided by :data:`workspace`.  The maximum number of
   subintervals is given by :data:`limit`, which may not exceed the allocated
   size of the workspace.

   Those subintervals with "large" widths :math:`d` where :math:`d\omega > 4` are
   computed using a 25-point Clenshaw-Curtis integration rule, which handles the
   oscillatory behavior.  Subintervals with a "small" widths where
   :math:`d\omega < 4` are computed using a 15-point Gauss-Kronrod integration.

QAWF adaptive integration for Fourier integrals
===============================================
.. index::
   single: QAWF quadrature algorithm
   single: Fourier integrals, numerical

.. function:: int gsl_integration_qawf (gsl_function * f, const double a, const double epsabs, const size_t limit, gsl_integration_workspace * workspace, gsl_integration_workspace * cycle_workspace, gsl_integration_qawo_table * wf, double * result, double * abserr)

   This function attempts to compute a Fourier integral of the function
   :data:`f` over the semi-infinite interval :math:`[a,+\infty)`

   .. only:: not texinfo

      .. math::

         I = \int_a^{+\infty} dx f(x)
         \left\{
         \begin{array}{c}
         \sin{(\omega x)} \\
         \cos{(\omega x)}
         \end{array}
         \right\}

   .. only:: texinfo

      ::

        I = \int_a^{+\infty} dx f(x) sin(omega x)
        I = \int_a^{+\infty} dx f(x) cos(omega x)

   The parameter :math:`\omega` and choice of :math:`\sin` or :math:`\cos` is
   taken from the table :data:`wf` (the length :data:`L` can take any value,
   since it is overridden by this function to a value appropriate for the
   Fourier integration).  The integral is computed using the QAWO algorithm
   over each of the subintervals,

   .. only:: not texinfo

      .. math::

         C_1 &= [a,a+c] \\
         C_2 &= [a+c,a+2c] \\
         \dots &= \dots \\
         C_k &= [a+(k-1)c,a+kc]

   .. only:: texinfo

      ::

        C_1 = [a, a + c]
        C_2 = [a + c, a + 2 c]
        ... = ...
        C_k = [a + (k-1) c, a + k c]

   where 
   :math:`c = (2 floor(|\omega|) + 1) \pi/|\omega|`.  The width :math:`c` is
   chosen to cover an odd number of periods so that the contributions from
   the intervals alternate in sign and are monotonically decreasing when
   :data:`f` is positive and monotonically decreasing.  The sum of this
   sequence of contributions is accelerated using the epsilon-algorithm.

   This function works to an overall absolute tolerance of
   :data:`abserr`.  The following strategy is used: on each interval
   :math:`C_k` the algorithm tries to achieve the tolerance

   .. math:: TOL_k = u_k abserr

   where 
   :math:`u_k = (1 - p)p^{k-1}` and :math:`p = 9/10`.
   The sum of the geometric series of contributions from each interval
   gives an overall tolerance of :data:`abserr`.

   If the integration of a subinterval leads to difficulties then the
   accuracy requirement for subsequent intervals is relaxed,

   .. only:: not texinfo

      .. math:: TOL_k = u_k \max{(abserr, \max_{i<k}{(E_i)})}

   .. only:: texinfo

      .. math:: TOL_k = u_k max(abserr, max_{i<k}(E_i))

   where :math:`E_k` is the estimated error on the interval :math:`C_k`.

   The subintervals and their results are stored in the memory provided by
   :data:`workspace`.  The maximum number of subintervals is given by
   :data:`limit`, which may not exceed the allocated size of the workspace.
   The integration over each subinterval uses the memory provided by
   :data:`cycle_workspace` as workspace for the QAWO algorithm.

.. index::
   single: cquad, doubly-adaptive integration

CQUAD doubly-adaptive integration
=================================

|cquad| is a new doubly-adaptive general-purpose quadrature
routine which can handle most types of singularities,
non-numerical function values such as :code:`Inf` or :code:`NaN`,
as well as some divergent integrals. It generally requires more
function evaluations than the integration routines in
|quadpack|, yet fails less often for difficult integrands.

The underlying algorithm uses a doubly-adaptive scheme in which
Clenshaw-Curtis quadrature rules of increasing degree are used
to compute the integral in each interval. The :math:`L_2`-norm of
the difference between the underlying interpolatory polynomials
of two successive rules is used as an error estimate. The
interval is subdivided if the difference between two successive
rules is too large or a rule of maximum degree has been reached.
     
.. index:: gsl_integration_cquad_workspace

.. function:: gsl_integration_cquad_workspace * gsl_integration_cquad_workspace_alloc (size_t n) 

   This function allocates a workspace sufficient to hold the data for
   :data:`n` intervals. The number :data:`n` is not the maximum number of
   intervals that will be evaluated. If the workspace is full, intervals
   with smaller error estimates will be discarded. A minimum of 3
   intervals is required and for most functions, a workspace of size 100
   is sufficient.

.. function:: void gsl_integration_cquad_workspace_free (gsl_integration_cquad_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_integration_cquad (const gsl_function * f, double a, double b, double epsabs, double epsrel, gsl_integration_cquad_workspace * workspace,  double * result, double * abserr, size_t * nevals)

   This function computes the integral of :math:`f` over :math:`(a,b)`
   within the desired absolute and relative error limits, :data:`epsabs`
   and :data:`epsrel` using the |cquad| algorithm.  The function returns
   the final approximation, :data:`result`, an estimate of the absolute
   error, :data:`abserr`, and the number of function evaluations required,
   :data:`nevals`.

   The |cquad| algorithm divides the integration region into
   subintervals, and in each iteration, the subinterval with the largest
   estimated error is processed.  The algorithm uses Clenshaw-Curtis
   quadrature rules of degree 4, 8, 16 and 32 over 5, 9, 17 and 33 nodes
   respectively. Each interval is initialized with the lowest-degree
   rule. When an interval is processed, the next-higher degree rule is
   evaluated and an error estimate is computed based on the
   :math:`L_2`-norm of the difference between the underlying interpolating
   polynomials of both rules. If the highest-degree rule has already been
   used, or the interpolatory polynomials differ significantly, the
   interval is bisected.

   The subintervals and their results are stored in the memory 
   provided by :data:`workspace`. If the error estimate or the number of 
   function evaluations is not needed, the pointers :data:`abserr` and :data:`nevals`
   can be set to :code:`NULL`.

Romberg integration
===================

The Romberg integration method estimates the definite integral

.. math:: I = \int_a^b f(x) dx

by applying Richardson extrapolation on the trapezoidal rule, using
equally spaced points with spacing

.. math:: h_k = (b - a) 2^{-k}

for :math:`k = 1, \dots, n`. For each :math:`k`, Richardson extrapolation
is used :math:`k-1` times on previous approximations to improve the order
of accuracy as much as possible. Romberg integration typically works
well (and converges quickly) for smooth integrands with no singularities in
the interval or at the end points.

.. function:: gsl_integration_romberg_workspace * gsl_integration_romberg_alloc(const size_t n)

   This function allocates a workspace for Romberg integration, specifying
   a maximum of :math:`n` iterations, or divisions of the interval. Since
   the number of divisions is :math:`2^n + 1`, :math:`n` can be kept relatively
   small (i.e. :math:`10` or :math:`20`). It is capped at a maximum value of
   :math:`30` to prevent overflow. The size of the workspace is :math:`O(2n)`.

.. function:: void gsl_integration_romberg_free(gsl_integration_romberg_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_integration_romberg(const gsl_function * f, const double a, const double b, const double epsabs, const double epsrel, double * result, size_t * neval, gsl_integration_romberg_workspace * w)

   This function integrates :math:`f(x)`, specified by :data:`f`, from :data:`a` to
   :data:`b`, storing the answer in :data:`result`. At each step in the iteration,
   convergence is tested by checking:

   .. math:: | I_k - I_{k-1} | \le \textrm{max} \left( epsabs, epsrel \times |I_k| \right)

   where :math:`I_k` is the current approximation and :math:`I_{k-1}` is the approximation
   of the previous iteration. If the method does not converge within the previously
   specified :math:`n` iterations, the function stores the best current estimate in
   :data:`result` and returns :macro:`GSL_EMAXITER`. If the method converges, the function
   returns :macro:`GSL_SUCCESS`. The total number of function evaluations is returned
   in :data:`neval`.

Gauss-Legendre integration
==========================

The fixed-order Gauss-Legendre integration routines are provided for fast
integration of smooth functions with known polynomial order.  The
:math:`n`-point Gauss-Legendre rule is exact for polynomials of order
:math:`2n - 1` or less.  For example, these rules are useful when integrating
basis functions to form mass matrices for the Galerkin method.  Unlike other
numerical integration routines within the library, these routines do not accept
absolute or relative error bounds.

.. index:: gsl_integration_glfixed_table

.. function:: gsl_integration_glfixed_table * gsl_integration_glfixed_table_alloc (size_t n)

   This function determines the Gauss-Legendre abscissae and weights necessary for
   an :math:`n`-point fixed order integration scheme.  If possible, high precision
   precomputed coefficients are used.  If precomputed weights are not available,
   lower precision coefficients are computed on the fly.

.. function:: double gsl_integration_glfixed (const gsl_function * f, double a, double b, const gsl_integration_glfixed_table * t)

   This function applies the Gauss-Legendre integration rule
   contained in table :data:`t` and returns the result.

.. function:: int gsl_integration_glfixed_point (double a, double b, size_t i, double * xi, double * wi, const gsl_integration_glfixed_table * t)

   For :data:`i` in :math:`[0, \dots, n - 1]`, this function obtains the
   :data:`i`-th Gauss-Legendre point :data:`xi` and weight :data:`wi` on the interval
   [:data:`a`, :data:`b`].  The points and weights are ordered by increasing point
   value.  A function :math:`f` may be integrated on [:data:`a`, :data:`b`] by summing
   :math:`wi * f(xi)` over :data:`i`.

.. index:: gsl_integration_glfixed_table

.. function:: void gsl_integration_glfixed_table_free (gsl_integration_glfixed_table * t)

   This function frees the memory associated with the table :data:`t`.

Fixed point quadratures
=======================

.. index::
   single: interpolating quadrature
   single: quadrature, interpolating
   single: quadrature, fixed point

The routines in this section approximate an integral by the sum

.. math:: \int_a^b w(x) f(x) dx = \sum_{i=1}^n w_i f(x_i)

where :math:`f(x)` is the function to be integrated and :math:`w(x)` is
a weighting function. The :math:`n` weights :math:`w_i` and nodes :math:`x_i` are carefully chosen
so that the result is exact when :math:`f(x)` is a polynomial of degree :math:`2n - 1`
or less. Once the user chooses the order :math:`n` and weighting function :math:`w(x)`, the
weights :math:`w_i` and nodes :math:`x_i` can be precomputed and used to efficiently evaluate
integrals for any number of functions :math:`f(x)`.

This method works best when :math:`f(x)` is well approximated by a polynomial on the interval
:math:`(a,b)`, and so is not suitable for functions with singularities.
Since the user specifies ahead of time how many quadrature nodes will be used, these
routines do not accept absolute or relative error bounds.  The table below lists
the weighting functions currently supported.

.. _tab_fixed-quadratures:

================ ======================== =========================================== =======================================================
Name             Interval                 Weighting function :math:`w(x)`             Constraints
================ ======================== =========================================== =======================================================
Legendre         :math:`(a,b)`            :math:`1`                                   :math:`b > a`
Chebyshev Type 1 :math:`(a,b)`            :math:`1 / \sqrt{(b - x) (x - a)}`          :math:`b > a`
Gegenbauer       :math:`(a,b)`            :math:`((b - x) (x - a))^{\alpha}`          :math:`\alpha > -1, b > a`
Jacobi           :math:`(a,b)`            :math:`(b - x)^{\alpha} (x - a)^{\beta}`    :math:`\alpha,\beta > -1, b > a`
Laguerre         :math:`(a,\infty)`       :math:`(x-a)^{\alpha} \exp{( -b (x - a) )}` :math:`\alpha > -1, b > 0`
Hermite          :math:`(-\infty,\infty)` :math:`|x-a|^{\alpha} \exp{( -b (x-a)^2 )}` :math:`\alpha > -1, b > 0`
Exponential      :math:`(a,b)`            :math:`|x - (a + b)/2|^{\alpha}`            :math:`\alpha > -1, b > a`
Rational         :math:`(a,\infty)`       :math:`(x - a)^{\alpha} (x + b)^{\beta}`    :math:`\alpha > -1, \alpha + \beta + 2n < 0, a + b > 0`
Chebyshev Type 2 :math:`(a,b)`            :math:`\sqrt{(b - x) (x - a)}`              :math:`b > a`
================ ======================== =========================================== =======================================================

The fixed point quadrature routines use the following workspace to store the nodes and weights,
as well as additional variables for intermediate calculations:

.. type:: gsl_integration_fixed_workspace

   This workspace is used for fixed point quadrature rules and looks like this::

     typedef struct
     {
       size_t n;        /* number of nodes/weights */
       double *weights; /* quadrature weights */
       double *x;       /* quadrature nodes */
       double *diag;    /* diagonal of Jacobi matrix */
       double *subdiag; /* subdiagonal of Jacobi matrix */
       const gsl_integration_fixed_type * type;
     } gsl_integration_fixed_workspace;

.. function:: gsl_integration_fixed_workspace * gsl_integration_fixed_alloc(const gsl_integration_fixed_type * T, const size_t n, const double a, const double b, const double alpha, const double beta)

   This function allocates a workspace for computing integrals with interpolating quadratures using :data:`n`
   quadrature nodes. The parameters :data:`a`, :data:`b`, :data:`alpha`, and :data:`beta` specify the integration
   interval and/or weighting function for the various quadrature types. See the :ref:`table <tab_fixed-quadratures>` above
   for constraints on these parameters. The size of the workspace is :math:`O(4n)`.

   .. type:: gsl_integration_fixed_type

      The type of quadrature used is specified by :data:`T` which can be set to the following choices:

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_legendre

         This specifies Legendre quadrature integration. The parameters :data:`alpha` and
         :data:`beta` are ignored for this type.

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_chebyshev

         This specifies Chebyshev type 1 quadrature integration. The parameters :data:`alpha` and
         :data:`beta` are ignored for this type.

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_gegenbauer

         This specifies Gegenbauer quadrature integration. The parameter :data:`beta` is ignored for this type.

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_jacobi

         This specifies Jacobi quadrature integration.

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_laguerre

         This specifies Laguerre quadrature integration. The parameter :data:`beta` is ignored for this type.

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_hermite

         This specifies Hermite quadrature integration. The parameter :data:`beta` is ignored for this type.

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_exponential

         This specifies exponential quadrature integration. The parameter :data:`beta` is ignored for this type.

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_rational

         This specifies rational quadrature integration.

      .. var:: gsl_integration_fixed_type * gsl_integration_fixed_chebyshev2

         This specifies Chebyshev type 2 quadrature integration. The parameters :data:`alpha` and
         :data:`beta` are ignored for this type.

.. function:: void gsl_integration_fixed_free(gsl_integration_fixed_workspace * w)

   This function frees the memory assocated with the workspace :data:`w`

.. function:: size_t gsl_integration_fixed_n(const gsl_integration_fixed_workspace * w)

   This function returns the number of quadrature nodes and weights.

.. function:: double * gsl_integration_fixed_nodes(const gsl_integration_fixed_workspace * w)

   This function returns a pointer to an array of size :data:`n` containing the quadrature nodes :math:`x_i`.

.. function:: double * gsl_integration_fixed_weights(const gsl_integration_fixed_workspace * w)

   This function returns a pointer to an array of size :data:`n` containing the quadrature weights :math:`w_i`.

.. function:: int gsl_integration_fixed(const gsl_function * func, double * result, const gsl_integration_fixed_workspace * w)

   This function integrates the function :math:`f(x)` provided in :data:`func` using previously
   computed fixed quadrature rules. The integral is approximated as
   
   .. math:: \sum_{i=1}^n w_i f(x_i)

   where :math:`w_i` are the quadrature weights and :math:`x_i` are the quadrature nodes computed
   previously by :func:`gsl_integration_fixed_alloc`. The sum is stored in :data:`result` on output.

Integrating on the unit sphere
==============================

.. index::
   single: quadrature, on unit sphere
   single: quadrature, Lebedev

This section contains routines to calculate the surface integral of a
function over the unit sphere,

.. math:: I[f] = \int d\Omega f(\Omega) = \int_0^{\pi} \sin{\theta} d\theta \int_0^{2\pi} d\phi f(\theta,\phi)

Lebedev developed a quadrature scheme to approximate this integral using a
single sum,

.. math:: I[f] \approx 4 \pi \sum_{i=1}^n w_i f(\theta_i, \phi_i)

for appropriately chosen weights :math:`w_i` and nodes :math:`(\theta_i,\phi_i)`.
The Lebedev nodes are chosen to lie on the unit sphere and be invariant
under the octahedral rotation group with inversion.

The number of quadrature nodes :math:`n` is often chosen in order to exactly
integrate a certain degree spherical harmonic function :math:`Y_l^m`.
A general rule of thumb for integrating spherical harmonics up to degree
and order :math:`L` is to choose the number of nodes as,

.. math:: n \approx \frac{(L+1)^2}{3}

Calculating the Lebedev weights and nodes requires solving a set of nonlinear
equations. These equations have been solved, and the nodes and weights have
been tabulated for integrating spherical harmonics up to degree and order 131.
GSL offers a smaller subset of 32 quadrature rules, which are listed in
the table below.

=================================== ======================================
Spherical Harmonic degree :math:`L` Quadrature weights and nodes :math:`n`
=================================== ======================================
3                                   6
5                                   14
7                                   26
9                                   38
11                                  50
13                                  74
15                                  86
17                                  110
19                                  146
21                                  170
23                                  194
25                                  230
27                                  266
29                                  302
31                                  350
35                                  434
41                                  590
47                                  770
53                                  974
59                                  1202
65                                  1454
71                                  1730
77                                  2030
83                                  2354
89                                  2702
95                                  3074
101                                 3470
107                                 3890
113                                 4334
119                                 4802
125                                 5294
131                                 5810
=================================== ======================================

.. type:: gsl_integration_lebedev_workspace

   This workspace is used for Lebedev quadrature rules and looks like this::

     typedef struct
     {
       size_t n;        /* number of nodes/weights */
       double *weights; /* quadrature weights */
       double *x;       /* x quadrature nodes */
       double *y;       /* y quadrature nodes */
       double *z;       /* z quadrature nodes */
       double *theta;   /* theta quadrature nodes */
       double *phi;     /* phi quadrature nodes */
     } gsl_integration_lebedev_workspace;

   The arrays :data:`x`, :data:`y`, :data:`z` of length :math:`n`
   contain the Cartesian coordinates :math:`(x_i,y_i,z_i)` of the
   Lebedev nodes which lie on the unit sphere. The arrays
   :data:`theta`, :data:`phi` contain the spherical coordinates
   :math:`(\theta_i,\phi_i)` of the same nodes on the unit sphere.

.. function:: gsl_integration_lebedev_workspace * gsl_integration_lebedev_alloc(const size_t n)

   This function allocates a workspace for a Lebedev quadrature
   rule of size :data:`n` and computes the nodes and weights.
   The size of the workspace is :math:`O(6n)`.

   If the input :data:`n` does not match one of the rules in the
   table above, the error code :macro:`GSL_EDOM` is returned.

   Here is some example code for integrating a function :math:`f(\theta,\phi)`
   with Lebedev quadrature::

     const size_t n = 230; /* integrate exactly up to spherical harmonic degree 25 */
     gsl_integration_lebedev_workspace * w = gsl_integration_lebedev_alloc(n);
     double result = 0.0;
     size_t i;

     for (i = 0; i < n; ++i)
       result += w->weights[i] * f(w->theta[i], w->phi[i]);

     result *= 4.0 * M_PI;
     gsl_integration_lebedev_free(w);

.. function:: void gsl_integration_lebedev_free(gsl_integration_lebedev_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: size_t gsl_integration_lebedev_n(const gsl_integration_lebedev_workspace * w)

   This function returns the number of quadrature nodes associated
   with the workspace :data:`w`.

Error codes
===========

In addition to the standard error codes for invalid arguments the
functions can return the following values,

===================== ============================================================================================================
:macro:`GSL_EMAXITER` the maximum number of subdivisions was exceeded.
:macro:`GSL_EROUND`   cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table.
:macro:`GSL_ESING`    a non-integrable singularity or other bad integrand behavior was found in the integration interval.
:macro:`GSL_EDIVERGE` the integral is divergent, or too slowly convergent to be integrated numerically.
:macro:`GSL_EDOM`     error in the values of the input arguments
===================== ============================================================================================================

Examples
========

Adaptive integration example
----------------------------

The integrator :code:`QAGS` will handle a large class of definite
integrals.  For example, consider the following integral, which has an
algebraic-logarithmic singularity at the origin,

.. math:: \int_0^1 x^{-1/2} \log(x) dx = -4

The program below computes this integral to a relative accuracy bound of
:code:`1e-7`.

.. include:: examples/integration.c
   :code:

The results below show that the desired accuracy is achieved after 8
subdivisions. 

.. include:: examples/integration.txt
   :code:

In fact, the extrapolation procedure used by :code:`QAGS` produces an
accuracy of almost twice as many digits.  The error estimate returned by
the extrapolation procedure is larger than the actual error, giving a
margin of safety of one order of magnitude.

Fixed-point quadrature example
------------------------------

In this example, we use a fixed-point quadrature rule to integrate the
integral

.. math::
   
   \int_{-\infty}^{\infty} e^{-x^2} \left( x^m + 1 \right) dx =
     \left\{
       \begin{array}{cc}
         \sqrt{\pi} + \Gamma{\left( \frac{m+1}{2} \right)}, & m \textrm{ even} \\
         \sqrt{\pi}, & m \textrm{ odd}
       \end{array}
     \right.

for integer :math:`m`. Consulting our :ref:`table <tab_fixed-quadratures>` of fixed point quadratures,
we see that this integral can be evaluated with a Hermite quadrature rule,
setting :math:`\alpha = 0, a = 0, b = 1`. Since we are integrating a polynomial
of degree :math:`m`, we need to choose the number of nodes :math:`n \ge (m+1)/2`
to achieve the best results.

First we will try integrating for :math:`m = 10, n = 5`, which does not satisfy
our criteria above::

  $ ./integration2 10 5

The output is,

.. include:: examples/integration2a.txt
   :code:

So, we find a large error. Now we try integrating for :math:`m = 10, n = 6` which
does satisfy the criteria above::

  $ ./integration2 10 6

The output is,

.. include:: examples/integration2b.txt
   :code:

The program is given below.

.. include:: examples/integration2.c
   :code:

References and Further Reading
==============================

The following book is the definitive reference for |quadpack|, and was
written by the original authors.  It provides descriptions of the
algorithms, program listings, test programs and examples.  It also
includes useful advice on numerical integration and many references to
the numerical integration literature used in developing |quadpack|.

* R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D.K. Kahaner.
  |quadpack| A subroutine package for automatic integration
  Springer Verlag, 1983.

The |cquad| integration algorithm is described in the following paper:

* P. Gonnet, "Increasing the Reliability of Adaptive Quadrature Using
  Explicit Interpolants", ACM Transactions on Mathematical
  Software, Volume 37 (2010), Issue 3, Article 26.

The fixed-point quadrature routines are based on IQPACK, described in the
following papers:

* S. Elhay, J. Kautsky, Algorithm 655: IQPACK, FORTRAN Subroutines for the
  Weights of Interpolatory Quadrature, ACM Transactions on Mathematical Software,
  Volume 13, Number 4, December 1987, pages 399-415.

* J. Kautsky, S. Elhay, Calculation of the Weights of Interpolatory Quadratures,
  Numerische Mathematik, Volume 40, Number 3, October 1982, pages 407-422.

The Lebedev quadrature routines are based on the paper:

* Lebedev, V. I. and Laikov, D. N. (1999). A quadrature formula for the sphere of
  the 131st algebraic order of accuracy. In Doklady Mathematics (Vol. 59, No. 3, pp. 477-481).
