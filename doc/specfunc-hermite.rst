.. Version 1: Konrad Griessinger (konradg(at)gmx.net), 12/2013

Hermite polynomials and functions are discussed in Abramowitz & Stegun, Chapter 22 and
Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society.
The Hermite polynomials and functions are defined in the header file :file:`gsl_sf_hermite.h`.

.. index::
   single: Hermite polynomials

Hermite Polynomials
-------------------

The Hermite polynomials exist in two variants: the physicist version
:math:`H_n(x)` and the probabilist version :math:`He_n(x)`.
They are defined by the derivatives

.. only:: not texinfo

   .. math::

      H_n(x) & = (-1)^n e^{x^2} \left({d \over dx}\right)^n e^{-x^2} \\
      He_n(x) & = (-1)^n e^{x^2/2} \left({d \over dx}\right)^n e^{-x^2/2}

.. only:: texinfo

   ::

      H_n(x) = (-1)^n e^{x^2} (d / dx)^n e^{-x^2} 
      He_n(x) = (-1)^n e^{x^2/2} (d / dx)^n e^{-x^2/2} 

They are connected via 

.. only:: not texinfo

   .. math::

      H_n(x) & = 2^{n/2} He_n \left( \sqrt{2} x \right) \\
      He_n(x) & = 2^{-n/2} H_n \left( {x \over \sqrt{2}} \right)

.. only:: texinfo

   ::

      H_n(x) = 2^{n/2} He_n(\sqrt{2} x)
      He_n(x) = 2^{-n/2} H_n(x / \sqrt{2})

and satisfy the ordinary differential equations

.. only:: not texinfo

   .. math::

      H_n^{\prime\prime}(x) - 2x H_n^{\prime}(x) + 2n H_n(x) & = 0 \\
      He_n^{\prime\prime}(x) - x He_n^{\prime}(x) + n He_n(x) & = 0

.. only:: texinfo

   ::

      H_n^{''}(x) - 2x H_n^{'}(x) + 2n H_n(x) = 0
      He_n^{''}(x) - x He_n^{'}(x) + n He_n(x) = 0

.. function:: double gsl_sf_hermite (const int n, const double x)
              int gsl_sf_hermite_e (const int n, const double x, gsl_sf_result * result)

   These routines evaluate the physicist Hermite polynomial :math:`H_n(x)` of order :data:`n` at position :data:`x`.
   If an overflow is detected, :macro:`GSL_EOVRFLW` is returned without calling the error handler.

.. function:: int gsl_sf_hermite_array (const int nmax, const double x, double * result_array)

   This routine evaluates all physicist Hermite polynomials :math:`H_n` up to order :data:`nmax` at position :data:`x`.
   The results are stored in :data:`result_array`.

.. function:: double gsl_sf_hermite_series (const int n, const double x, const double * a)
              int gsl_sf_hermite_series_e (const int n, const double x, const double * a, gsl_sf_result * result)

   These routines evaluate the series :math:`\sum_{j=0}^n a_j H_j(x)` with :math:`H_j` being
   the :math:`j`-th physicist Hermite polynomial using the Clenshaw algorithm.

.. function:: double gsl_sf_hermite_prob (const int n, const double x)
              int gsl_sf_hermite_prob_e (const int n, const double x, gsl_sf_result * result)

   These routines evaluate the probabilist Hermite polynomial :math:`He_n(x)` of order :data:`n` at position :data:`x`.
   If an overflow is detected, :macro:`GSL_EOVRFLW` is returned without calling the error handler.

.. function:: int gsl_sf_hermite_prob_array (const int nmax, const double x, double * result_array)

   This routine evaluates all probabilist Hermite polynomials :math:`He_n(x)` up to order :data:`nmax` at position :data:`x`.
   The results are stored in :data:`result_array`.

.. function:: double gsl_sf_hermite_prob_series (const int n, const double x, const double * a)
              int gsl_sf_hermite_prob_series_e (const int n, const double x, const double * a, gsl_sf_result * result)

   These routines evaluate the series :math:`\sum_{j=0}^n a_j He_j(x)` with :math:`He_j` being the
   :math:`j`-th probabilist Hermite polynomial using the Clenshaw algorithm.

Derivatives of Hermite Polynomials
----------------------------------
.. index::
   single: Hermite polynomials, derivatives

.. function:: double gsl_sf_hermite_deriv (const int m, const int n, const double x)
              int gsl_sf_hermite_deriv_e (const int m, const int n, const double x, gsl_sf_result * result)

   These routines evaluate the :data:`m`-th derivative of the physicist Hermite polynomial :math:`H_n(x)` of order :data:`n`
   at position :data:`x`.

.. function::  int gsl_sf_hermite_array_deriv (const int m, const int nmax, const double x, double * result_array)

   This routine evaluates the :data:`m`-th derivative of all physicist Hermite polynomials :math:`H_n(x)` from
   orders :math:`0, \dots, \text{nmax}` at position :data:`x`.
   The result :math:`d^m/dx^m H_n(x)` is stored in :code:`result_array[n]`. The output
   :data:`result_array` must have length at least :code:`nmax + 1`.

.. function:: int gsl_sf_hermite_deriv_array (const int mmax, const int n, const double x, double * result_array)

   This routine evaluates all derivative orders from :math:`0, \dots, \text{mmax}` of the
   physicist Hermite polynomial of order :data:`n`, :math:`H_n`, at position :data:`x`.
   The result :math:`d^m/dx^m H_n(x)` is stored in :code:`result_array[m]`. The output
   :data:`result_array` must have length at least :code:`mmax + 1`.

.. function:: double gsl_sf_hermite_prob_deriv (const int m, const  int n, const double x)
              int gsl_sf_hermite_prob_deriv_e (const int m, const  int n, const double x, gsl_sf_result * result)

   These routines evaluate the :data:`m`-th derivative of the probabilist Hermite polynomial :math:`He_n(x)`
   of order :data:`n` at position :data:`x`.

.. function:: int gsl_sf_hermite_prob_array_deriv (const int m, const int nmax, const double x, double * result_array)

   This routine evaluates the :data:`m`-th derivative of all probabilist Hermite polynomials :math:`He_n(x)` from
   orders :math:`0, \dots, \text{nmax}` at position :data:`x`.
   The result :math:`d^m/dx^m He_n(x)` is stored in :code:`result_array[n]`. The output
   :data:`result_array` must have length at least :code:`nmax + 1`.

.. function:: int gsl_sf_hermite_prob_deriv_array (const int mmax, const int n, const double x, double * result_array)

   This routine evaluates all derivative orders from :math:`0, \dots, \text{mmax}` of the
   probabilist Hermite polynomial of order :data:`n`, :math:`He_n`, at position :data:`x`.
   The result :math:`d^m/dx^m He_n(x)` is stored in :code:`result_array[m]`. The output
   :data:`result_array` must have length at least :code:`mmax + 1`.

.. index::
   single: Hermite functions

Hermite Functions
-----------------

The Hermite functions are defined by

.. only:: not texinfo

   .. math:: \psi_n(x) = \left( 2^n n! \sqrt{\pi} \right)^{-1/2} e^{-x^2/2} H_n \left( x \right)

.. only:: texinfo

   ::

      \psi_n(x) = ( 2^n n! \sqrt{\pi} )^{-1/2} e^{-x^2/2} H_n(x)

and satisfy the Schrödinger equation for a quantum mechanical harmonic oscillator

.. only:: not texinfo

   .. math:: \psi_n^{\prime\prime}(x) + (2n + 1 - x^2) \psi_n(x) = 0

.. only:: texinfo

   ::

      psi''_n(x) + (2n + 1 - x^2) psi_n(x) = 0

They are orthonormal,

.. math:: \int_{-\infty}^{\infty} \psi_m(x) \psi_n(x) dx = \delta_{mn}

and form an orthonormal basis of :math:`L^2(\mathbb{R})`. The Hermite functions
are also eigenfunctions of the continuous Fourier transform. GSL offers two
methods for evaluating the Hermite functions. The first uses the standard three-term
recurrence relation which has :math:`O(n)` complexity and is the most accurate. The
second uses a Cauchy integral approach due to Bunck (2009) which has :math:`O(\sqrt{n})`
complexity which represents a significant speed improvement for large :math:`n`, although
it is slightly less accurate.

.. function:: double gsl_sf_hermite_func (const int n, const double x)
              int gsl_sf_hermite_func_e (const int n, const double x, gsl_sf_result * result)

   These routines evaluate the Hermite function :math:`\psi_n(x)` of order :data:`n` at position :data:`x`
   using a three term recurrence relation. The algorithm complexity is :math:`O(n)`.

.. function:: double gsl_sf_hermite_func_fast (const int n, const double x)
              int gsl_sf_hermite_func_fast_e (const int n, const double x, gsl_sf_result * result)

   These routines evaluate the Hermite function :math:`\psi_n(x)` of order :data:`n` at position :data:`x`
   using a the Cauchy integral algorithm due to Bunck, 2009. The algorithm complexity is :math:`O(\sqrt{n})`.

.. function:: int gsl_sf_hermite_func_array (const int nmax, const double x, double * result_array)

   This routine evaluates all Hermite functions :math:`\psi_n(x)` for orders :math:`n = 0, \dots, \textrm{nmax}`
   at position :data:`x`, using the recurrence relation algorithm. The results are stored in
   :data:`result_array` which has length at least :code:`nmax + 1`.

.. function:: double gsl_sf_hermite_func_series (const int n, const double x, const double * a)
              int gsl_sf_hermite_func_series_e (const int n, const double x, const double * a, gsl_sf_result * result)

   These routines evaluate the series :math:`\sum_{j=0}^n a_j \psi_j(x)` with :math:`\psi_j` being
   the :math:`j`-th Hermite function using the Clenshaw algorithm.

Derivatives of Hermite Functions
--------------------------------
.. index::
   single: Hermite functions, derivatives

.. function:: double gsl_sf_hermite_func_der (const int m, const int n, const double x)
              int gsl_sf_hermite_func_der_e (const int m, const int n, const double x, gsl_sf_result * result)

   These routines evaluate the :data:`m`-th derivative of the Hermite function :math:`\psi_n(x)` of order :data:`n` at position :data:`x`.

Zeros of Hermite Polynomials and Hermite Functions
--------------------------------------------------
.. index::
   single: Hermite polynomials, zeros
   single: Hermite functions, zeros

These routines calculate the :math:`s`-th zero of the Hermite polynomial/function of order
:math:`n`. Since the zeros are symmetrical around zero, only positive zeros are calculated,
ordered from smallest to largest, starting from index 1. Only for odd polynomial orders a
zeroth zero exists, its value always being zero.

.. function:: double gsl_sf_hermite_zero (const int n, const int s)
              int gsl_sf_hermite_zero_e (const int n, const int s, gsl_sf_result * result)

   These routines evaluate the :data:`s`-th zero of the physicist Hermite polynomial :math:`H_n(x)` of order :data:`n`.

.. function:: double gsl_sf_hermite_prob_zero (const int n, const int s)
              int gsl_sf_hermite_prob_zero_e (const int n, const int s, gsl_sf_result * result)

   These routines evaluate the :data:`s`-th zero of the probabilist Hermite polynomial :math:`He_n(x)` of order :data:`n`.

.. function:: double gsl_sf_hermite_func_zero (const int n, const int s)
              int gsl_sf_hermite_func_zero_e (const int n, const int s, gsl_sf_result * result)

   These routines evaluate the :data:`s`-th zero of the Hermite function :math:`\psi_n(x)` of order :data:`n`.
