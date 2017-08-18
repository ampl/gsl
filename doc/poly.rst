.. index::
   single: polynomials, roots of

***********
Polynomials
***********

.. include:: include.rst

This chapter describes functions for evaluating and solving polynomials.
There are routines for finding real and complex roots of quadratic and
cubic equations using analytic methods.  An iterative polynomial solver
is also available for finding the roots of general polynomials with real
coefficients (of any order).  The functions are declared in the header
file :file:`gsl_poly.h`.

.. index::
   single: polynomial evaluation
   single: evaluation of polynomials

Polynomial Evaluation
=====================

The functions described here evaluate the polynomial 

.. only:: not texinfo

   .. math::

      P(x) = c[0] + c[1] x + c[2] x^2 + \dots + c[len-1] x^{len-1}

.. only:: texinfo

   P(x) = c[0] + c[1] x + c[2] x^2 + ... + c[len-1] x^{len-1}
   
using Horner's method for stability. |inlinefns|

.. function:: double gsl_poly_eval (const double c[], const int len, const double x)

   This function evaluates a polynomial with real coefficients for the real variable :data:`x`.

.. function:: gsl_complex gsl_poly_complex_eval (const double c[], const int len, const gsl_complex z)

   This function evaluates a polynomial with real coefficients for the complex variable :data:`z`.

.. function:: gsl_complex gsl_complex_poly_complex_eval (const gsl_complex c[], const int len, const gsl_complex z)

   This function evaluates a polynomial with complex coefficients for the complex variable :data:`z`.

.. function:: int gsl_poly_eval_derivs (const double c[], const size_t lenc, const double x, double res[], const size_t lenres)

   This function evaluates a polynomial and its derivatives storing the
   results in the array :data:`res` of size :data:`lenres`.  The output array
   contains the values of :math:`d^k P(x)/d x^k` for the specified value of
   :data:`x` starting with :math:`k = 0`.

.. index::
   single: divided differences, polynomials
   single: evaluation of polynomials, in divided difference form

Divided Difference Representation of Polynomials
================================================

The functions described here manipulate polynomials stored in Newton's
divided-difference representation.  The use of divided-differences is
described in Abramowitz & Stegun sections 25.1.4 and 25.2.26, and
Burden and Faires, chapter 3, and discussed briefly below.

Given a function :math:`f(x)`, an :math:`n`\ th degree interpolating polynomial :math:`P_{n}(x)`
can be constructed which agrees with :math:`f` at :math:`n+1` distinct points
:math:`x_0,x_1,...,x_{n}`. This polynomial can be written in a
form known as Newton's divided-difference representation

.. only:: not texinfo

   .. math::

      P_{n}(x) = f(x_0) + \sum_{k=1}^n [x_0,x_1,...,x_k] (x-x_0)(x-x_1) \cdots (x-x_{k-1})

.. only:: texinfo

   P_{n}(x) = f(x_0) + \sum_{k=1}^n [x_0,x_1,...,x_k] (x-x_0)(x-x_1) ... (x-x_{k-1})

where the divided differences :math:`[x_0,x_1,...,x_k]` are defined in section 25.1.4 of
Abramowitz and Stegun. Additionally, it is possible to construct an interpolating
polynomial of degree :math:`2n+1` which also matches the first derivatives of :math:`f`
at the points :math:`x_0,x_1,...,x_n`. This is called the Hermite interpolating
polynomial and is defined as

.. only:: not texinfo

   .. math::

      H_{2n+1}(x) = f(z_0) + \sum_{k=1}^{2n+1} [z_0,z_1,...,z_k] (x-z_0)(x-z_1) \cdots (x-z_{k-1})

.. only:: texinfo

   H_{2n+1}(x) = f(z_0) + \sum_{k=1}^{2n+1} [z_0,z_1,...,z_k] (x-z_0)(x-z_1) ... (x-z_{k-1})

where the elements of :math:`z = \{x_0,x_0,x_1,x_1,...,x_n,x_n\}` are defined by
:math:`z_{2k} = z_{2k+1} = x_k`. The divided-differences :math:`[z_0,z_1,...,z_k]`
are discussed in Burden and Faires, section 3.4.

.. function:: int gsl_poly_dd_init (double dd[], const double xa[], const double ya[], size_t size)

   This function computes a divided-difference representation of the
   interpolating polynomial for the points :math:`(x, y)` stored in
   the arrays :data:`xa` and :data:`ya` of length :data:`size`.  On output the
   divided-differences of (:data:`xa`, :data:`ya`) are stored in the array
   :data:`dd`, also of length :data:`size`. Using the notation above,
   :math:`dd[k] = [x_0,x_1,...,x_k]`.

.. function:: double gsl_poly_dd_eval (const double dd[], const double xa[], const size_t size, const double x)

   This function evaluates the polynomial stored in divided-difference form
   in the arrays :data:`dd` and :data:`xa` of length :data:`size` at the point
   :data:`x`. |inlinefn|

.. function:: int gsl_poly_dd_taylor (double c[], double xp, const double dd[], const double xa[], size_t size, double w[])

   This function converts the divided-difference representation of a
   polynomial to a Taylor expansion.  The divided-difference representation
   is supplied in the arrays :data:`dd` and :data:`xa` of length :data:`size`.
   On output the Taylor coefficients of the polynomial expanded about the
   point :data:`xp` are stored in the array :data:`c` also of length
   :data:`size`.  A workspace of length :data:`size` must be provided in the
   array :data:`w`.

.. function:: int gsl_poly_dd_hermite_init (double dd[], double za[], const double xa[], const double ya[], const double dya[], const size_t size)

   This function computes a divided-difference representation of the
   interpolating Hermite polynomial for the points :math:`(x,y)` stored in
   the arrays :data:`xa` and :data:`ya` of length :data:`size`. Hermite interpolation
   constructs polynomials which also match first derivatives :math:`dy/dx` which are
   provided in the array :data:`dya` also of length :data:`size`. The first derivatives can be
   incorported into the usual divided-difference algorithm by forming a new
   dataset :math:`z = \{x_0,x_0,x_1,x_1,...\}`, which is stored in the array
   :data:`za` of length 2*\ :data:`size` on output. On output the
   divided-differences of the Hermite representation are stored in the array
   :data:`dd`, also of length 2*\ :data:`size`. Using the notation above,
   :math:`dd[k] = [z_0,z_1,...,z_k]`. The resulting Hermite polynomial
   can be evaluated by calling :func:`gsl_poly_dd_eval` and using
   :data:`za` for the input argument :data:`xa`.

.. index::
   single: quadratic equation, solving

Quadratic Equations
===================

.. function:: int gsl_poly_solve_quadratic (double a, double b, double c, double * x0, double * x1)

   This function finds the real roots of the quadratic equation,

   .. math::

      a x^2 + b x + c = 0

   The number of real roots (either zero, one or two) is returned, and
   their locations are stored in :data:`x0` and :data:`x1`.  If no real roots
   are found then :data:`x0` and :data:`x1` are not modified.  If one real root
   is found (i.e. if :math:`a=0`) then it is stored in :data:`x0`.  When two
   real roots are found they are stored in :data:`x0` and :data:`x1` in
   ascending order.  The case of coincident roots is not considered
   special.  For example :math:`(x-1)^2=0` will have two roots, which happen
   to have exactly equal values.

   The number of roots found depends on the sign of the discriminant
   :math:`b^2 - 4 a c`.  This will be subject to rounding and cancellation
   errors when computed in double precision, and will also be subject to
   errors if the coefficients of the polynomial are inexact.  These errors
   may cause a discrete change in the number of roots.  However, for
   polynomials with small integer coefficients the discriminant can always
   be computed exactly.

.. function:: int gsl_poly_complex_solve_quadratic (double a, double b, double c, gsl_complex * z0, gsl_complex * z1)

   This function finds the complex roots of the quadratic equation,

   .. math::

      a z^2 + b z + c = 0

   The number of complex roots is returned (either one or two) and the
   locations of the roots are stored in :data:`z0` and :data:`z1`.  The roots
   are returned in ascending order, sorted first by their real components
   and then by their imaginary components.  If only one real root is found
   (i.e. if :math:`a=0`) then it is stored in :data:`z0`.

.. index::
   single: cubic equation, solving

Cubic Equations
===============

.. function:: int gsl_poly_solve_cubic (double a, double b, double c, double * x0, double * x1, double * x2)

   This function finds the real roots of the cubic equation,

   .. math::

      x^3 + a x^2 + b x + c = 0

   with a leading coefficient of unity.  The number of real roots (either
   one or three) is returned, and their locations are stored in :data:`x0`,
   :data:`x1` and :data:`x2`.  If one real root is found then only :data:`x0`
   is modified.  When three real roots are found they are stored in
   :data:`x0`, :data:`x1` and :data:`x2` in ascending order.  The case of
   coincident roots is not considered special.  For example, the equation
   :math:`(x-1)^3=0` will have three roots with exactly equal values.  As
   in the quadratic case, finite precision may cause equal or
   closely-spaced real roots to move off the real axis into the complex
   plane, leading to a discrete change in the number of real roots.

.. function:: int gsl_poly_complex_solve_cubic (double a, double b, double c, gsl_complex * z0, gsl_complex * z1, gsl_complex * z2)

   This function finds the complex roots of the cubic equation,

   .. math::

      z^3 + a z^2 + b z + c = 0

   The number of complex roots is returned (always three) and the locations
   of the roots are stored in :data:`z0`, :data:`z1` and :data:`z2`.  The roots
   are returned in ascending order, sorted first by their real components
   and then by their imaginary components.

.. index::
   single: general polynomial equations, solving

General Polynomial Equations
============================

The roots of polynomial equations cannot be found analytically beyond
the special cases of the quadratic, cubic and quartic equation.  The
algorithm described in this section uses an iterative method to find the
approximate locations of roots of higher order polynomials.

.. type:: gsl_poly_complex_workspace

   This workspace contains parameters used for finding roots of general polynomials

.. function:: gsl_poly_complex_workspace * gsl_poly_complex_workspace_alloc (size_t n)

   This function allocates space for a :type:`gsl_poly_complex_workspace`
   struct and a workspace suitable for solving a polynomial with :data:`n`
   coefficients using the routine :func:`gsl_poly_complex_solve`.

   The function returns a pointer to the newly allocated
   :type:`gsl_poly_complex_workspace` if no errors were detected, and a null
   pointer in the case of error.

.. function:: void gsl_poly_complex_workspace_free (gsl_poly_complex_workspace * w)

   This function frees all the memory associated with the workspace
   :data:`w`.

.. function:: int gsl_poly_complex_solve (const double * a, size_t n, gsl_poly_complex_workspace * w, gsl_complex_packed_ptr z)

   This function computes the roots of the general polynomial 

   .. only:: not texinfo

      .. math::

         P(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_{n-1} x^{n-1}

   .. only:: texinfo

      P(x) = a_0 + a_1 x + a_2 x^2 + ... + a_{n-1} x^{n-1}
      
   using balanced-QR reduction of the companion matrix.  The parameter :data:`n`
   specifies the length of the coefficient array.  The coefficient of the
   highest order term must be non-zero.  The function requires a workspace
   :data:`w` of the appropriate size.  The :math:`n-1` roots are returned in
   the packed complex array :data:`z` of length :math:`2(n-1)`, alternating
   real and imaginary parts.

   The function returns :data:`GSL_SUCCESS` if all the roots are found. If
   the QR reduction does not converge, the error handler is invoked with
   an error code of :data:`GSL_EFAILED`.  Note that due to finite precision,
   roots of higher multiplicity are returned as a cluster of simple roots
   with reduced accuracy.  The solution of polynomials with higher-order
   roots requires specialized algorithms that take the multiplicity
   structure into account (see e.g. Z. Zeng, Algorithm 835, ACM
   Transactions on Mathematical Software, Volume 30, Issue 2 (2004), pp
   218--236).

Examples
========

To demonstrate the use of the general polynomial solver we will take the
polynomial :math:`P(x) = x^5 - 1` which has these roots:

.. only:: not texinfo

   .. math::

      1, e^{2\pi i / 5}, e^{4\pi i / 5}, e^{6\pi i / 5}, e^{8\pi i / 5}

.. only:: texinfo

   1, e^{2*pi i / 5}, e^{4*pi i / 5}, e^{6*pi i / 5}, e^{8*pi i / 5}

The following program will find these roots.

.. include:: examples/polyroots.c
   :code:

The output of the program is

.. include:: examples/polyroots.txt
   :code:

which agrees with the analytic result, :math:`z_n = \exp(2 \pi n i/5)`.

References and Further Reading
==============================

The balanced-QR method and its error analysis are described in the
following papers,

* R.S. Martin, G. Peters and J.H. Wilkinson, "The QR Algorithm for Real
  Hessenberg Matrices", Numerische Mathematik, 14 (1970), 219--231.

* B.N. Parlett and C. Reinsch, "Balancing a Matrix for Calculation of
  Eigenvalues and Eigenvectors", Numerische Mathematik, 13 (1969),
  293--304.

* A. Edelman and H. Murakami, "Polynomial roots from companion matrix
  eigenvalues", Mathematics of Computation, Vol.: 64, No.: 210
  (1995), 763--776.

The formulas for divided differences are given in the following texts,

* Abramowitz and Stegun, Handbook of Mathematical Functions,
  Sections 25.1.4 and 25.2.26.

* R. L. Burden and J. D. Faires, Numerical Analysis, 9th edition,
  ISBN 0-538-73351-9, 2011.
