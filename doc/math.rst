.. index::
   single: elementary functions
   single: mathematical functions, elementary

**********************
Mathematical Functions
**********************

This chapter describes basic mathematical functions.  Some of these
functions are present in system libraries, but the alternative versions
given here can be used as a substitute when the system functions are not
available.

The functions and macros described in this chapter are defined in the
header file :file:`gsl_math.h`.

.. index::
   single: mathematical constants, defined as macros
   single: numerical constants, defined as macros
   single: constants, mathematical (defined as macros)
   single: macros for mathematical constants

Mathematical Constants
======================

The library ensures that the standard BSD mathematical constants
are defined. For reference, here is a list of the constants:

.. index::
   single: e, defined as a macro
   single: pi, defined as a macro
   single: Euler's constant, defined as a macro

===================== ===================================
:macro:`M_E`          The base of exponentials, :math:`e`
:macro:`M_LOG2E`      The base-2 logarithm of :math:`e`, :math:`\log_2 (e)`
:macro:`M_LOG10E`     The base-10 logarithm of :math:`e`, :math:`\log_{10} (e)`
:macro:`M_SQRT2`      The square root of two, :math:`\sqrt 2`
:macro:`M_SQRT1_2`    The square root of one-half, :math:`\sqrt{1/2}`
:macro:`M_SQRT3`      The square root of three, :math:`\sqrt 3`
:macro:`M_PI`         The constant pi, :math:`\pi`
:macro:`M_PI_2`       Pi divided by two, :math:`\pi/2`
:macro:`M_PI_4`       Pi divided by four, :math:`\pi/4`
:macro:`M_SQRTPI`     The square root of pi, :math:`\sqrt\pi`
:macro:`M_2_SQRTPI`   Two divided by the square root of pi, :math:`2/\sqrt\pi`
:macro:`M_1_PI`       The reciprocal of pi, :math:`1/\pi`
:macro:`M_2_PI`       Twice the reciprocal of pi, :math:`2/\pi`
:macro:`M_LN10`       The natural logarithm of ten, :math:`\ln(10)`
:macro:`M_LN2`        The natural logarithm of two, :math:`\ln(2)`
:macro:`M_LNPI`       The natural logarithm of pi, :math:`\ln(\pi)`
:macro:`M_EULER`      Euler's constant, :math:`\gamma`
===================== ===================================

.. index::
   single: infinity, defined as a macro
   single: IEEE infinity, defined as a macro

Infinities and Not-a-number
===========================

.. macro:: GSL_POSINF

   This macro contains the IEEE representation of positive infinity,
   :math:`+\infty`. It is computed from the expression :code:`+1.0/0.0`.

.. macro:: GSL_NEGINF

   This macro contains the IEEE representation of negative infinity,
   :math:`-\infty`. It is computed from the expression :code:`-1.0/0.0`.

.. index::
   single: NaN, defined as a macro
   single: Not-a-number, defined as a macro
   single: IEEE NaN, defined as a macro

.. macro:: GSL_NAN

   This macro contains the IEEE representation of the Not-a-Number symbol,
   :code:`NaN`. It is computed from the ratio :code:`0.0/0.0`.

.. function:: int gsl_isnan (const double x)

   This function returns 1 if :data:`x` is not-a-number.

.. function:: int gsl_isinf (const double x)

   This function returns :math:`+1` if :data:`x` is positive infinity,
   :math:`-1` if :data:`x` is negative infinity and 0
   otherwise. [#f1]_

.. function:: int gsl_finite (const double x)

   This function returns 1 if :data:`x` is a real number, and 0 if it is
   infinite or not-a-number.

Elementary Functions
====================

The following routines provide portable implementations of functions
found in the BSD math library.  When native versions are not available
the functions described here can be used instead.  The substitution can
be made automatically if you use :code:`autoconf` to compile your
application (see :ref:`portability-functions`).

.. index::
   single: log1p
   single: logarithm, computed accurately near 1

.. function:: double gsl_log1p (const double x)

   This function computes the value of :math:`\log(1+x)` in a way that is
   accurate for small :data:`x`. It provides an alternative to the BSD math
   function :code:`log1p(x)`.

.. index::
   single: expm1
   single: exponential, difference from 1 computed accurately

.. function:: double gsl_expm1 (const double x)

   This function computes the value of :math:`\exp(x)-1` in a way that is
   accurate for small :data:`x`. It provides an alternative to the BSD math
   function :code:`expm1(x)`.

.. index::
   single: hypot
   single: euclidean distance function, hypot
   single: length, computed accurately using hypot

.. function:: double gsl_hypot (const double x, const double y)

   This function computes the value of
   :math:`\sqrt{x^2 + y^2}` in a way that avoids overflow. It provides an
   alternative to the BSD math function :code:`hypot(x,y)`.

.. index::
   single: euclidean distance function, hypot3
   single: length, computed accurately using hypot3

.. function:: double gsl_hypot3 (const double x, const double y, const double z)

   This function computes the value of
   :math:`\sqrt{x^2 + y^2 + z^2}` in a way that avoids overflow.

.. index::
   single: acosh
   single: hyperbolic cosine, inverse
   single: inverse hyperbolic cosine

.. function:: double gsl_acosh (const double x)

   This function computes the value of :math:`\arccosh{(x)}`. It provides an
   alternative to the standard math function :code:`acosh(x)`.

.. index::
   single: asinh
   single: hyperbolic sine, inverse
   single: inverse hyperbolic sine

.. function:: double gsl_asinh (const double x)

   This function computes the value of :math:`\arcsinh{(x)}`. It provides an
   alternative to the standard math function :code:`asinh(x)`.

.. index::
   single: atanh
   single: hyperbolic tangent, inverse
   single: inverse hyperbolic tangent

.. function:: double gsl_atanh (const double x)

   This function computes the value of :math:`\arctanh{(x)}`. It provides an
   alternative to the standard math function :code:`atanh(x)`.

.. index:: ldexp

.. function:: double gsl_ldexp (double x, int e)

   This function computes the value of :math:`x * 2^e`. It provides an
   alternative to the standard math function :code:`ldexp(x,e)`.

.. index:: frexp

.. function:: double gsl_frexp (double x, int * e)

   This function splits the number :data:`x` into its normalized fraction
   :math:`f` and exponent :math:`e`, such that :math:`x = f * 2^e` and
   :math:`0.5 <= f < 1`. The function returns :math:`f` and stores the
   exponent in :math:`e`. If :math:`x` is zero, both :math:`f` and :math:`e`
   are set to zero. This function provides an alternative to the standard
   math function :code:`frexp(x, e)`.

Small integer powers
====================

A common complaint about the standard C library is its lack of a
function for calculating (small) integer powers.  GSL provides some simple
functions to fill this gap.  For reasons of efficiency, these functions
do not check for overflow or underflow conditions. 

.. function::
   double gsl_pow_int (double x, int n)
   double gsl_pow_uint (double x, unsigned int n)

   These routines computes the power :math:`x^n` for integer :data:`n`.  The
   power is computed efficiently---for example, :math:`x^8` is computed as
   :math:`((x^2)^2)^2`, requiring only 3 multiplications.  A version of this
   function which also computes the numerical error in the result is
   available as :func:`gsl_sf_pow_int_e`.

.. function::
   double gsl_pow_2 (const double x)
   double gsl_pow_3 (const double x)
   double gsl_pow_4 (const double x)
   double gsl_pow_5 (const double x)
   double gsl_pow_6 (const double x)
   double gsl_pow_7 (const double x)
   double gsl_pow_8 (const double x)
   double gsl_pow_9 (const double x)

   These functions can be used to compute small integer powers :math:`x^2`,
   :math:`x^3`, etc. efficiently. The functions will be inlined when 
   :macro:`HAVE_INLINE` is defined, so that use of these functions 
   should be as efficient as explicitly writing the corresponding 
   product expression::

     #include <gsl/gsl_math.h>
     double y = gsl_pow_4 (3.141)  /* compute 3.141**4 */

Testing the Sign of Numbers
===========================

.. macro:: GSL_SIGN (x)

   This macro returns the sign of :data:`x`. It is defined as :code:`((x) >= 0
   ? 1 : -1)`. Note that with this definition the sign of zero is positive
   (regardless of its IEEE sign bit).

Testing for Odd and Even Numbers
================================

.. macro:: GSL_IS_ODD (n)

   This macro evaluates to 1 if :data:`n` is odd and 0 if :data:`n` is
   even. The argument :data:`n` must be of integer type.

.. macro:: GSL_IS_EVEN (n)

   This macro is the opposite of :macro:`GSL_IS_ODD`. It evaluates to 1 if
   :data:`n` is even and 0 if :data:`n` is odd. The argument :data:`n` must be of
   integer type.

Maximum and Minimum functions
=============================

Note that the following macros perform multiple evaluations of their
arguments, so they should not be used with arguments that have side
effects (such as a call to a random number generator).

.. index:: maximum of two numbers

.. macro:: GSL_MAX (a, b)

   This macro returns the maximum of :data:`a` and :data:`b`. It is defined
   as :code:`((a) > (b) ? (a):(b))`.

.. index:: minimum of two numbers

.. macro:: GSL_MIN (a, b)

   This macro returns the minimum of :data:`a` and :data:`b`. It is defined as 
   :code:`((a) < (b) ? (a):(b))`.

.. function:: extern inline double GSL_MAX_DBL (double a, double b)

   This function returns the maximum of the double precision numbers
   :data:`a` and :data:`b` using an inline function. The use of a function
   allows for type checking of the arguments as an extra safety feature. On
   platforms where inline functions are not available the macro
   :macro:`GSL_MAX` will be automatically substituted.

.. function:: extern inline double GSL_MIN_DBL (double a, double b)

   This function returns the minimum of the double precision numbers
   :data:`a` and :data:`b` using an inline function. The use of a function
   allows for type checking of the arguments as an extra safety feature. On
   platforms where inline functions are not available the macro
   :macro:`GSL_MIN` will be automatically substituted.

.. function::
   extern inline int GSL_MAX_INT (int a, int b)
   extern inline int GSL_MIN_INT (int a, int b)

   These functions return the maximum or minimum of the integers :data:`a`
   and :data:`b` using an inline function.  On platforms where inline
   functions are not available the macros :macro:`GSL_MAX` or :macro:`GSL_MIN`
   will be automatically substituted.

.. function::
   extern inline long double GSL_MAX_LDBL (long double a, long double b)
   extern inline long double GSL_MIN_LDBL (long double a, long double b)

   These functions return the maximum or minimum of the long doubles :data:`a`
   and :data:`b` using an inline function.  On platforms where inline
   functions are not available the macros :macro:`GSL_MAX` or :macro:`GSL_MIN`
   will be automatically substituted.

Approximate Comparison of Floating Point Numbers
================================================

It is sometimes useful to be able to compare two floating point numbers
approximately, to allow for rounding and truncation errors.  The following
function implements the approximate floating-point comparison algorithm
proposed by D.E. Knuth in Section 4.2.2 of "Seminumerical
Algorithms" (3rd edition).

.. index::
   single: approximate comparison of floating point numbers
   single: safe comparison of floating point numbers
   single: floating point numbers, approximate comparison

.. function:: int gsl_fcmp (double x, double y, double epsilon)

   This function determines whether :data:`x` and :data:`y` are approximately
   equal to a relative accuracy :data:`epsilon`.

   The relative accuracy is measured using an interval of size :math:`2
   \delta`, where :math:`\delta = 2^k \epsilon` and :math:`k` is the
   maximum base-2 exponent of :math:`x` and :math:`y` as computed by the
   function :func:`frexp`.

   If :math:`x` and :math:`y` lie within this interval, they are considered
   approximately equal and the function returns 0. Otherwise if :math:`x <
   y`, the function returns :math:`-1`, or if :math:`x > y`, the function returns
   :math:`+1`.

   Note that :math:`x` and :math:`y` are compared to relative accuracy, so
   this function is not suitable for testing whether a value is
   approximately zero. 

   The implementation is based on the package :code:`fcmp` by T.C. Belding.

.. rubric:: Footnotes

.. [#f1] Note that the C99 standard only requires the
   system :func:`isinf` function to return a non-zero value, without the
   sign of the infinity.  The implementation in some earlier versions of
   GSL used the system :func:`isinf` function and may have this behavior
   on some platforms.  Therefore, it is advisable to test the sign of
   :data:`x` separately, if needed, rather than relying the sign of the
   return value from :func:`gsl_isinf()`.
