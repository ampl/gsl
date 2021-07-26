.. index:: complex numbers

***************
Complex Numbers
***************

The functions described in this chapter provide support for complex
numbers.  The algorithms take care to avoid unnecessary intermediate
underflows and overflows, allowing the functions to be evaluated over 
as much of the complex plane as possible. 

.. FIXME: this still needs to be
.. done for the csc,sec,cot,csch,sech,coth functions

For multiple-valued functions the branch cuts have been chosen to follow
the conventions of Abramowitz and Stegun.
The functions return principal values which are
the same as those in GNU Calc, which in turn are the same as those in
"Common Lisp, The Language (Second Edition)" [#f1]_
and the HP-28/48 series of calculators.

The complex types are defined in the header file :file:`gsl_complex.h`,
while the corresponding complex functions and arithmetic operations are
defined in :file:`gsl_complex_math.h`.

.. index::
   single: representations of complex numbers
   single: polar form of complex numbers
   single: gsl_complex

Representation of complex numbers
=================================

Complex numbers are represented using the type :code:`gsl_complex`. The
default interface defines :code:`gsl_complex` as::

    typedef struct
    {
      double dat[2];
    } gsl_complex;

The real and imaginary part are stored in contiguous elements of a two
element array. This eliminates any padding between the real and
imaginary parts, :code:`dat[0]` and :code:`dat[1]`, allowing the struct to
be mapped correctly onto packed complex arrays.

If a C compiler is available which supports the C11 standard,
and the :file:`<complex.h>` header file is included *prior* to :file:`gsl_complex.h`,
then :code:`gsl_complex` will be defined to be the native C complex type::

    typedef double complex gsl_complex

This allows users to use :code:`gsl_complex` in ordinary operations such as::

    gsl_complex x = 2 + 5 * I;
    gsl_complex y = x + (3 - 4*I);

.. important::

   Native C support for complex numbers was introduced in the C99 standard,
   and additional functionality was added in C11. When :file:`<complex.h>`
   is included in a user's program prior to :file:`gsl_complex.h`,
   GSL uses the new C11 functionality to define the :macro:`GSL_REAL` and
   :macro:`GSL_IMAG` macros. It does not appear possible to properly define
   these macros using the C99 standard, and so using a C99 compiler will not define
   :code:`gsl_complex` to the native complex type.

   Some compilers, such as the gcc 4.8 series implement only a portion of the C11
   standard and so they may fail to correctly compile GSL code when a user tries
   to turn on native complex functionality. A workaround for this issue is to
   either remove :file:`<complex.h>` from the include list, or add
   :code:`-DGSL_COMPLEX_LEGACY` to the compiler flags, which will use the older
   struct-based definition of :code:`gsl_complex`.

Complex number macros
=====================

The following C macros offer convenient ways to manipulate complex numbers.

.. macro::
   GSL_REAL (z)
   GSL_IMAG (z)

   These macros return a memory location (lvalue) corresponding to the real and imaginary
   parts respectively of the complex number :data:`z`. This allows users to perform
   operations like::

       gsl_complex x, y;

       GSL_REAL(x) = 4;
       GSL_IMAG(x) = 2;

       GSL_REAL(y) = GSL_REAL(x);
       GSL_IMAG(y) = GSL_REAL(x);

   In other words, these macros can both read and write to the real and imaginary
   parts of a complex variable.

.. macro:: GSL_SET_COMPLEX (zp, x, y)

   This macro uses the Cartesian components (:data:`x`, :data:`y`) to set the
   real and imaginary parts of the complex number pointed to by :data:`zp`.
   For example::

     GSL_SET_COMPLEX(&z, 3, 4)

   sets :math:`z` to be :math:`3 + 4i`.

.. .. macro::
..    GSL_SET_REAL (zp,x)
..    GSL_SET_IMAG (zp,y)

..    These macros allow the real and imaginary parts of the complex number
..    pointed to by :data:`zp` to be set independently.

Assigning complex numbers
=========================

.. function:: gsl_complex gsl_complex_rect (double x, double y)

   This function uses the rectangular Cartesian components
   :math:`(x,y)` to return the complex number :math:`z = x + i y`.
   An inline version of this function is used when :macro:`HAVE_INLINE`
   is defined.

.. function:: gsl_complex gsl_complex_polar (double r, double theta)

   This function returns the complex number :math:`z = r \exp(i \theta) = r
   (\cos(\theta) + i \sin(\theta))` from the polar representation
   (:data:`r`, :data:`theta`).

Properties of complex numbers
=============================

.. index:: argument of complex number

.. function:: double gsl_complex_arg (gsl_complex z)

   This function returns the argument of the complex number :data:`z`,
   :math:`\arg(z)`, where :math:`-\pi < \arg(z) <= \pi`.

.. index:: magnitude of complex number

.. function:: double gsl_complex_abs (gsl_complex z)

   This function returns the magnitude of the complex number :data:`z`, :math:`|z|`.

.. function:: double gsl_complex_abs2 (gsl_complex z)

   This function returns the squared magnitude of the complex number
   :data:`z`, :math:`|z|^2`.

.. function:: double gsl_complex_logabs (gsl_complex z)

   This function returns the natural logarithm of the magnitude of the
   complex number :data:`z`, :math:`\log|z|`.  It allows an accurate
   evaluation of :math:`\log|z|` when :math:`|z|` is close to one. The direct
   evaluation of :code:`log(gsl_complex_abs(z))` would lead to a loss of
   precision in this case.

.. index:: complex arithmetic

Complex arithmetic operators
============================

.. function:: gsl_complex gsl_complex_add (gsl_complex a, gsl_complex b)

   This function returns the sum of the complex numbers :data:`a` and
   :data:`b`, :math:`z=a+b`.

.. function:: gsl_complex gsl_complex_sub (gsl_complex a, gsl_complex b)

   This function returns the difference of the complex numbers :data:`a` and
   :data:`b`, :math:`z=a-b`.

.. function:: gsl_complex gsl_complex_mul (gsl_complex a, gsl_complex b)

   This function returns the product of the complex numbers :data:`a` and
   :data:`b`, :math:`z=ab`.

.. function:: gsl_complex gsl_complex_div (gsl_complex a, gsl_complex b)

   This function returns the quotient of the complex numbers :data:`a` and
   :data:`b`, :math:`z=a/b`.

.. function:: gsl_complex gsl_complex_add_real (gsl_complex a, double x)

   This function returns the sum of the complex number :data:`a` and the
   real number :data:`x`, :math:`z=a+x`.

.. function:: gsl_complex gsl_complex_sub_real (gsl_complex a, double x)

   This function returns the difference of the complex number :data:`a` and the
   real number :data:`x`, :math:`z=a-x`.

.. function:: gsl_complex gsl_complex_mul_real (gsl_complex a, double x)

   This function returns the product of the complex number :data:`a` and the
   real number :data:`x`, :math:`z=ax`.

.. function:: gsl_complex gsl_complex_div_real (gsl_complex a, double x)

   This function returns the quotient of the complex number :data:`a` and the
   real number :data:`x`, :math:`z=a/x`.

.. function:: gsl_complex gsl_complex_add_imag (gsl_complex a, double y)

   This function returns the sum of the complex number :data:`a` and the
   imaginary number :math:`iy`, :math:`z=a+iy`.

.. function:: gsl_complex gsl_complex_sub_imag (gsl_complex a, double y)

   This function returns the difference of the complex number :data:`a` and the
   imaginary number :math:`iy`, :math:`z=a-iy`.

.. function:: gsl_complex gsl_complex_mul_imag (gsl_complex a, double y)

   This function returns the product of the complex number :data:`a` and the
   imaginary number :math:`iy`, :math:`z=a*(iy)`.

.. function:: gsl_complex gsl_complex_div_imag (gsl_complex a, double y)

   This function returns the quotient of the complex number :data:`a` and the
   imaginary number :math:`iy`, :math:`z=a/(iy)`.

.. index:: conjugate of complex number

.. function:: gsl_complex gsl_complex_conjugate (gsl_complex z)

   This function returns the complex conjugate of the complex number
   :data:`z`, :math:`z^* = x - i y`.

.. function:: gsl_complex gsl_complex_inverse (gsl_complex z)

   This function returns the inverse, or reciprocal, of the complex number
   :data:`z`, :math:`1/z = (x - i y)/(x^2 + y^2)`.

.. function:: gsl_complex gsl_complex_negative (gsl_complex z)

   This function returns the negative of the complex number
   :data:`z`, :math:`-z = (-x) + i(-y)`.


Elementary Complex Functions
============================

.. index:: square root of complex number

.. function:: gsl_complex gsl_complex_sqrt (gsl_complex z)

   This function returns the square root of the complex number :data:`z`,
   :math:`\sqrt z`. The branch cut is the negative real axis. The result
   always lies in the right half of the complex plane.

.. function:: gsl_complex gsl_complex_sqrt_real (double x)

   This function returns the complex square root of the real number
   :data:`x`, where :data:`x` may be negative.

.. index::
   single: power of complex number
   single: exponentiation of complex number

.. function:: gsl_complex gsl_complex_pow (gsl_complex z, gsl_complex a)

   The function returns the complex number :data:`z` raised to the complex
   power :data:`a`, :math:`z^a`. This is computed as :math:`\exp(\log(z)*a)`
   using complex logarithms and complex exponentials.

.. function:: gsl_complex gsl_complex_pow_real (gsl_complex z, double x)

   This function returns the complex number :data:`z` raised to the real
   power :data:`x`, :math:`z^x`.

.. function:: gsl_complex gsl_complex_exp (gsl_complex z)

   This function returns the complex exponential of the complex number
   :data:`z`, :math:`\exp(z)`.

.. index:: logarithm of complex number

.. function:: gsl_complex gsl_complex_log (gsl_complex z)

   This function returns the complex natural logarithm (base :math:`e`) of
   the complex number :data:`z`, :math:`\log(z)`.  The branch cut is the
   negative real axis. 

.. function:: gsl_complex gsl_complex_log10 (gsl_complex z)

   This function returns the complex base-10 logarithm of
   the complex number :data:`z`, :math:`\log_{10} (z)`.

.. function:: gsl_complex gsl_complex_log_b (gsl_complex z, gsl_complex b)

   This function returns the complex base-:data:`b` logarithm of the complex
   number :data:`z`, :math:`\log_b(z)`. This quantity is computed as the ratio
   :math:`\log(z)/\log(b)`.

.. index:: trigonometric functions of complex numbers

Complex Trigonometric Functions
===============================

.. index::
   single: sin, of complex number

.. function:: gsl_complex gsl_complex_sin (gsl_complex z)

   This function returns the complex sine of the complex number :data:`z`,
   :math:`\sin(z) = (\exp(iz) - \exp(-iz))/(2i)`.

.. index:: cosine of complex number

.. function:: gsl_complex gsl_complex_cos (gsl_complex z)

   This function returns the complex cosine of the complex number :data:`z`,
   :math:`\cos(z) = (\exp(iz) + \exp(-iz))/2`.

.. index:: tangent of complex number

.. function:: gsl_complex gsl_complex_tan (gsl_complex z)

   This function returns the complex tangent of the complex number :data:`z`,
   :math:`\tan(z) = \sin(z)/\cos(z)`.

.. function:: gsl_complex gsl_complex_sec (gsl_complex z)

   This function returns the complex secant of the complex number :data:`z`,
   :math:`\sec(z) = 1/\cos(z)`.

.. function:: gsl_complex gsl_complex_csc (gsl_complex z)

   This function returns the complex cosecant of the complex number :data:`z`,
   :math:`\csc(z) = 1/\sin(z)`.

.. function:: gsl_complex gsl_complex_cot (gsl_complex z)

   This function returns the complex cotangent of the complex number :data:`z`,
   :math:`\cot(z) = 1/\tan(z)`.

.. index:: inverse complex trigonometric functions

Inverse Complex Trigonometric Functions
=======================================

.. function:: gsl_complex gsl_complex_arcsin (gsl_complex z)

   This function returns the complex arcsine of the complex number :data:`z`,
   :math:`\arcsin(z)`. The branch cuts are on the real axis, less than :math:`-1`
   and greater than :math:`1`.

.. function:: gsl_complex gsl_complex_arcsin_real (double z)

   This function returns the complex arcsine of the real number :data:`z`,
   :math:`\arcsin(z)`. For :math:`z` between :math:`-1` and :math:`1`, the
   function returns a real value in the range :math:`[-\pi/2,\pi/2]`. For
   :math:`z` less than :math:`-1` the result has a real part of :math:`-\pi/2`
   and a positive imaginary part.  For :math:`z` greater than :math:`1` the
   result has a real part of :math:`\pi/2` and a negative imaginary part.

.. function:: gsl_complex gsl_complex_arccos (gsl_complex z)

   This function returns the complex arccosine of the complex number :data:`z`,
   :math:`\arccos(z)`. The branch cuts are on the real axis, less than :math:`-1`
   and greater than :math:`1`.

.. function:: gsl_complex gsl_complex_arccos_real (double z)

   This function returns the complex arccosine of the real number :data:`z`,
   :math:`\arccos(z)`. For :math:`z` between :math:`-1` and :math:`1`, the
   function returns a real value in the range :math:`[0,\pi]`. For :math:`z`
   less than :math:`-1` the result has a real part of :math:`\pi` and a
   negative imaginary part.  For :math:`z` greater than :math:`1` the result
   is purely imaginary and positive.

.. function:: gsl_complex gsl_complex_arctan (gsl_complex z)

   This function returns the complex arctangent of the complex number
   :data:`z`, :math:`\arctan(z)`. The branch cuts are on the imaginary axis,
   below :math:`-i` and above :math:`i`.

.. function:: gsl_complex gsl_complex_arcsec (gsl_complex z)

   This function returns the complex arcsecant of the complex number :data:`z`,
   :math:`\arcsec(z) = \arccos(1/z)`.

.. function:: gsl_complex gsl_complex_arcsec_real (double z)

   This function returns the complex arcsecant of the real number :data:`z`,
   :math:`\arcsec(z) = \arccos(1/z)`.

.. function:: gsl_complex gsl_complex_arccsc (gsl_complex z)

   This function returns the complex arccosecant of the complex number :data:`z`,
   :math:`\arccsc(z) = \arcsin(1/z)`.

.. function:: gsl_complex gsl_complex_arccsc_real (double z)

   This function returns the complex arccosecant of the real number :data:`z`,
   :math:`\arccsc(z) = \arcsin(1/z)`.

.. function:: gsl_complex gsl_complex_arccot (gsl_complex z)

   This function returns the complex arccotangent of the complex number :data:`z`,
   :math:`\arccot(z) = \arctan(1/z)`.

.. index::
   single: hyperbolic functions, complex numbers

Complex Hyperbolic Functions
============================

.. function:: gsl_complex gsl_complex_sinh (gsl_complex z)

   This function returns the complex hyperbolic sine of the complex number
   :data:`z`, :math:`\sinh(z) = (\exp(z) - \exp(-z))/2`.

.. function:: gsl_complex gsl_complex_cosh (gsl_complex z)

   This function returns the complex hyperbolic cosine of the complex number
   :data:`z`, :math:`\cosh(z) = (\exp(z) + \exp(-z))/2`.

.. function:: gsl_complex gsl_complex_tanh (gsl_complex z)

   This function returns the complex hyperbolic tangent of the complex number
   :data:`z`, :math:`\tanh(z) = \sinh(z)/\cosh(z)`.

.. function:: gsl_complex gsl_complex_sech (gsl_complex z)

   This function returns the complex hyperbolic secant of the complex
   number :data:`z`, :math:`\sech(z) = 1/\cosh(z)`.

.. function:: gsl_complex gsl_complex_csch (gsl_complex z)

   This function returns the complex hyperbolic cosecant of the complex
   number :data:`z`, :math:`\csch(z) = 1/\sinh(z)`.

.. function:: gsl_complex gsl_complex_coth (gsl_complex z)

   This function returns the complex hyperbolic cotangent of the complex
   number :data:`z`, :math:`\coth(z) = 1/\tanh(z)`.

.. index::
   single: inverse hyperbolic functions, complex numbers

Inverse Complex Hyperbolic Functions
====================================

.. function:: gsl_complex gsl_complex_arcsinh (gsl_complex z)

   This function returns the complex hyperbolic arcsine of the
   complex number :data:`z`, :math:`\arcsinh(z)`.  The branch cuts are on the
   imaginary axis, below :math:`-i` and above :math:`i`.

.. function:: gsl_complex gsl_complex_arccosh (gsl_complex z)

   This function returns the complex hyperbolic arccosine of the complex
   number :data:`z`, :math:`\arccosh(z)`.  The branch cut is on the real
   axis, less than :math:`1`.  Note that in this case we use the negative
   square root in formula 4.6.21 of Abramowitz & Stegun giving
   :math:`\arccosh(z)=\log(z-\sqrt{z^2-1})`.

.. function:: gsl_complex gsl_complex_arccosh_real (double z)

   This function returns the complex hyperbolic arccosine of
   the real number :data:`z`, :math:`\arccosh(z)`.

.. function:: gsl_complex gsl_complex_arctanh (gsl_complex z)

   This function returns the complex hyperbolic arctangent of the complex
   number :data:`z`, :math:`\arctanh(z)`.  The branch cuts are on the real
   axis, less than :math:`-1` and greater than :math:`1`.

.. function:: gsl_complex gsl_complex_arctanh_real (double z)

   This function returns the complex hyperbolic arctangent of the real
   number :data:`z`, :math:`\arctanh(z)`.

.. function:: gsl_complex gsl_complex_arcsech (gsl_complex z)

   This function returns the complex hyperbolic arcsecant of the complex
   number :data:`z`, :math:`\arcsech(z) = \arccosh(1/z)`.

.. function:: gsl_complex gsl_complex_arccsch (gsl_complex z)

   This function returns the complex hyperbolic arccosecant of the complex
   number :data:`z`, :math:`\arccsch(z) = \arcsinh(1/z)`.

.. function:: gsl_complex gsl_complex_arccoth (gsl_complex z)

   This function returns the complex hyperbolic arccotangent of the complex
   number :data:`z`, :math:`\arccoth(z) = \arctanh(1/z)`.

References and Further Reading
==============================

The implementations of the elementary and trigonometric functions are
based on the following papers,

* T. E. Hull, Thomas F. Fairgrieve, Ping Tak Peter Tang,
  "Implementing Complex Elementary Functions Using Exception
  Handling", ACM Transactions on Mathematical Software, Volume 20
  (1994), pp 215--244, Corrigenda, p553

* T. E. Hull, Thomas F. Fairgrieve, Ping Tak Peter Tang,
  "Implementing the complex arcsin and arccosine functions using exception
  handling", ACM Transactions on Mathematical Software, Volume 23
  (1997) pp 299--335

The general formulas and details of branch cuts can be found in the
following books,

* Abramowitz and Stegun, Handbook of Mathematical Functions,
  "Circular Functions in Terms of Real and Imaginary Parts", Formulas
  4.3.55--58,
  "Inverse Circular Functions in Terms of Real and Imaginary Parts",
  Formulas 4.4.37--39,
  "Hyperbolic Functions in Terms of Real and Imaginary Parts",
  Formulas 4.5.49--52,
  "Inverse Hyperbolic Functions---relation to Inverse Circular Functions",
  Formulas 4.6.14--19.

* Dave Gillespie, Calc Manual, Free Software Foundation, ISBN
  1-882114-18-3

.. rubric:: Footnotes

.. [#f1] Note that the first edition uses different definitions.
