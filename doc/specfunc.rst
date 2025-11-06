.. index:: special functions

*****************
Special Functions
*****************

This chapter describes the GSL special function library.  The library
includes routines for calculating the values of Airy functions, Bessel
functions, Clausen functions, Coulomb wave functions, Coupling
coefficients, the Dawson function, Debye functions, Dilogarithms,
Elliptic integrals, Jacobi elliptic functions, Error functions,
Exponential integrals, Fermi-Dirac functions, Gamma functions,
Gegenbauer functions, Hermite polynomials and functions, Hypergeometric functions, Laguerre functions,
Legendre functions and Spherical Harmonics, the Psi (Digamma) Function,
Synchrotron functions, Transport functions, Trigonometric functions and
Zeta functions.  Each routine also computes an estimate of the numerical
error in the calculated value of the function.

The functions in this chapter are declared in individual header files,
such as :file:`gsl_sf_airy.h`, :file:`gsl_sf_bessel.h`, etc.  The complete
set of header files can be included using the file :file:`gsl_sf.h`.

Usage
=====

The special functions are available in two calling conventions, a
*natural form* which returns the numerical value of the function and
an *error-handling form* which returns an error code.  The two types
of function provide alternative ways of accessing the same underlying
code.

The *natural form* returns only the value of the function and can be
used directly in mathematical expressions.  For example, the following
function call will compute the value of the Bessel function
:math:`J_0(x)`::

    double y = gsl_sf_bessel_J0 (x);

There is no way to access an error code or to estimate the error using
this method.  To allow access to this information the alternative
error-handling form stores the value and error in a modifiable argument::

    gsl_sf_result result;
    int status = gsl_sf_bessel_J0_e (x, &result);

The error-handling functions have the suffix :code:`_e`. The returned
status value indicates error conditions such as overflow, underflow or
loss of precision.  If there are no errors the error-handling functions
return :code:`GSL_SUCCESS`.

The gsl_sf_result struct
========================

The error handling form of the special functions always calculate an
error estimate along with the value of the result.  Therefore,
structures are provided for amalgamating a value and error estimate.
These structures are declared in the header file :file:`gsl_sf_result.h`.

The following struct contains value and error fields.

.. type:: gsl_sf_result

   ::

     typedef struct
     {
       double val;
       double err;
     } gsl_sf_result;

   The field :data:`val` contains the value and the field :data:`err` contains
   an estimate of the absolute error in the value.

In some cases, an overflow or underflow can be detected and handled by a
function.  In this case, it may be possible to return a scaling exponent
as well as an error/value pair in order to save the result from
exceeding the dynamic range of the built-in types.  The
following struct contains value and error fields as well
as an exponent field such that the actual result is obtained as
:code:`result * 10^(e10)`.

.. type:: gsl_sf_result_e10

   ::

     typedef struct
     {
       double val;
       double err;
       int    e10;
     } gsl_sf_result_e10;

Modes
=====

The goal of the library is to achieve double precision accuracy wherever
possible.  However the cost of evaluating some special functions to
double precision can be significant, particularly where very high order
terms are required.  In these cases a :code:`mode` argument, of type
:type:`gsl_mode_t` allows the
accuracy of the function to be reduced in order to improve performance.
The following precision levels are available for the mode argument,

.. type:: gsl_mode_t

   .. macro:: GSL_PREC_DOUBLE

      Double-precision, a relative accuracy of approximately :math:`2 * 10^{-16}`.

   .. macro:: GSL_PREC_SINGLE

      Single-precision, a relative accuracy of approximately :math:`10^{-7}`.

   .. macro:: GSL_PREC_APPROX

      Approximate values, a relative accuracy of approximately :math:`5 * 10^{-4}`.

The approximate mode provides the fastest evaluation at the lowest
accuracy.

Airy Functions and Derivatives
==============================
.. include:: specfunc-airy.rst

Bessel Functions
================
.. include:: specfunc-bessel.rst

Clausen Functions
=================
.. include:: specfunc-clausen.rst

Coulomb Functions
=================
.. include:: specfunc-coulomb.rst

Coupling Coefficients
=====================
.. include:: specfunc-coupling.rst

Dawson Function
===============
.. include:: specfunc-dawson.rst

Debye Functions
===============
.. include:: specfunc-debye.rst

.. _dilog-function:

Dilogarithm
===========
.. include:: specfunc-dilog.rst

Elementary Operations
=====================
.. include:: specfunc-elementary.rst

Elliptic Integrals
==================
.. include:: specfunc-ellint.rst

Elliptic Functions (Jacobi)
===========================
.. include:: specfunc-elljac.rst

Error Functions
===============
.. include:: specfunc-erf.rst

Exponential Functions
=====================
.. include:: specfunc-exp.rst

Exponential Integrals
=====================
.. include:: specfunc-expint.rst

Fermi-Dirac Function
====================
.. include:: specfunc-fermi-dirac.rst

Gamma and Beta Functions
========================
.. include:: specfunc-gamma.rst

Gegenbauer Functions
====================
.. include:: specfunc-gegenbauer.rst

Hermite Polynomials and Functions
=================================
.. include:: specfunc-hermite.rst

Hypergeometric Functions
========================
.. include:: specfunc-hyperg.rst

.. _laguerre-functions:

Laguerre Functions
==================
.. include:: specfunc-laguerre.rst

Lambert W Functions
===================
.. include:: specfunc-lambert.rst

Legendre Functions and Spherical Harmonics
==========================================
.. include:: specfunc-legendre.rst

.. Associated Legendre Functions and Spherical Harmonics
.. =====================================================
.. .. include:: specfunc-alf.rst

Logarithm and Related Functions
===============================
.. include:: specfunc-log.rst

Mathieu Functions
=================
.. include:: specfunc-mathieu.rst

Power Function
==============
.. include:: specfunc-pow-int.rst

Psi (Digamma) Function
======================
.. include:: specfunc-psi.rst

Synchrotron Functions
=====================
.. include:: specfunc-synchrotron.rst

Transport Functions
===================
.. include:: specfunc-transport.rst

Trigonometric Functions
=======================
.. include:: specfunc-trig.rst

Zeta Functions
==============
.. include:: specfunc-zeta.rst

Examples
========

Example 1: Bessel function :math:`J_0`
--------------------------------------

The following example demonstrates the use of the error handling form of
the special functions, in this case to compute the Bessel function
:math:`J_0(5.0)`,

.. include:: examples/specfun_e.c
   :code:

Here are the results of running the program,

.. include:: examples/specfun_e.txt
   :code:

The next program computes the same quantity using the natural form of
the function. In this case the error term :data:`result.err` and return
status are not accessible.

.. include:: examples/specfun.c
   :code:

The results of the function are the same,

.. include:: examples/specfun.txt
   :code:

Example 2: Associated Legendre Functions
----------------------------------------

The following example program outputs the spherical harmonic
normalized ALF :math:`Y_2^1(\cos{\theta})` and its first and second
derivatives with respect to :math:`\theta`. The analytic
expressions are,

.. math::
   
   Y_2^1(\cos{\theta}) &= -\frac{1}{2} \sqrt{\frac{15}{2 \pi}} \sin{\theta} \cos{\theta} \\
   \frac{d}{d\theta} Y_2^1(\cos{\theta}) &= -\frac{1}{2} \sqrt{\frac{15}{2 \pi}} \left( \cos^2{\theta} - \sin^2{\theta} \right) \\
   \frac{d^2}{d\theta^2} Y_2^1(\cos{\theta}) &= 2 \sqrt{\frac{15}{2 \pi}} \sin{\theta} \cos{\theta}

.. _fig_specfun2:

.. figure:: /images/specfun2.png

   Spherical harmonic normalized ALF :math:`Y_2^1` and its first and
   second derivatives with respect to :math:`\theta`.

The source code is given below.

.. include:: examples/specfun2.c
   :code:

References and Further Reading
==============================

The library follows the conventions of the following book where possible,

* Handbook of Mathematical Functions, edited by Abramowitz & Stegun,
  Dover,  ISBN 0486612724.

The following papers contain information on the algorithms used 
to compute the special functions,

.. index:: MISCFUN

* Allan J. MacLeod, MISCFUN: A software package to compute uncommon
  special functions.  ACM Trans. Math. Soft., vol.: 22,
  1996, 288--301

* Bosch, W., On the computation of derivatives of Legendre functions,
  Phys. Chem. Earth, 25 (9-11), pg. 655--659, 2000.

* Bunck, B. F., A fast algorithm for evaluation of normalized Hermite
  functions, BIT Numer. Math, 49: 281-295, 2009.

* B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, Springer, 2006.

* G.N. Watson, A Treatise on the Theory of Bessel Functions,
  2nd Edition (Cambridge University Press, 1944).

* G. Nemeth, Mathematical Approximations of Special Functions,
  Nova Science Publishers, ISBN 1-56072-052-2

* B.C. Carlson, Special Functions of Applied Mathematics (1977)

* N. M. Temme, Special Functions: An Introduction to the Classical
  Functions of Mathematical Physics (1996), ISBN 978-0471113133.

* W.J. Thompson, Atlas for Computing Mathematical Functions, John Wiley & Sons,
  New York (1997).

* Y.Y. Luke, Algorithms for the Computation of Mathematical Functions, Academic
  Press, New York (1977).

* S. A. Holmes and W. E. Featherstone, A unified approach to the Clenshaw
  summation and the recursive computation of very high degree and order
  normalised associated Legendre functions, Journal of Geodesy, 76,
  pg. 279-299, 2002.

* D. E. Winch, D. J. Ivers, J. P. R. Turner and R. J. Stening, Geomagnetism
  and Schmidt quasi-normalization, Geophys. J. Int., 160, 487-504, 2005.
