.. Version 1: Konrad Griessinger (konradg(at)gmx.net), 12/2013

.. index::
   single: Hermite polynomials
   single: Hermite functions

.. :math:`He_n(x)`
.. @math{H_n(x)}
.. how can you get greek characters in the index in Texinfo?!?
.. @cindex @math{psi_n(x)}

The Hermite polynomials exist in two variants: the probabilists' version :math:`He_n(x)`
and the physicists'version :math:`H_n(x)`. The are defined by the derivatives

.. only:: not texinfo

   .. math::

      He_n(x) & = (-1)^n e^{x^2/2} \left({d \over dx}\right)^n e^{-x^2/2} \\
      H_n(x) & = (-1)^n e^{x^2} \left({d \over dx}\right)^n e^{-x^2}

.. only:: texinfo

   ::

      He_n(x) = (-1)^n e^{x^2/2} (d / dx)^n e^{-x^2/2} 
      H_n(x) = (-1)^n e^{x^2} (d / dx)^n e^{-x^2} 

They are connected via 

.. only:: not texinfo

   .. math::

      He_n(x) & = 2^{-n/2} H_n \left( {x \over \sqrt{2}} \right) \\
      H_n(x) & = 2^{n/2} He_n \left( \sqrt{2} x \right)

.. only:: texinfo

   ::

      He_n(x) = 2^{-n/2} H_n(x / \sqrt{2})
      H_n(x) = 2^{n/2} He_n(\sqrt{2} x)

and satisfy the ordinary differential equations

.. only:: not texinfo

   .. math::

      He_n^{\prime\prime}(x) - x He_n^{\prime}(x) + n He_n(x) & = 0 \\
      H_n^{\prime\prime}(x) - 2x H_n^{\prime}(x) + 2n H_n(x) & = 0

.. only:: texinfo

   ::

      He_n^{''}(x) - x He_n^{'}(x) + n He_n(x) = 0
      H_n^{''}(x) - 2x H_n^{'}(x) + 2n H_n(x) = 0

The closely related Hermite functions are defined by 

.. only:: not texinfo

   .. math:: \psi_n(x) = \left( n! \sqrt{\pi} \right)^{-1/2} e^{-x^2/2} He_n \left( {\sqrt{2} x} \right)

.. only:: texinfo

   ::

      psi_n = (n! sqrt(\pi))^{-1/2} e^{-x^2/2} He_n({sqrt(2) x})

and satisfy the Schrödinger equation for a quantum mechanical harmonic oscillator

.. only:: not texinfo

   .. math:: \psi_n^{\prime\prime}(x) + (2n + 1 - x^2) \psi_n(x) = 0

.. only:: texinfo

   ::

      psi_n^{''}(x) + (2n + 1 - x^2) psi_n(x) = 0

Maybe most importantly, the Hermite functions :math:`\psi_n` are eigenfunctions of the (continuous) Fourier transform.

For further information see Abramowitz & Stegun, Chapter 22 and Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials,
American Mathematical Society. The Hermite polynomials and functions are defined in the header file :file:`gsl_sf_hermite.h`.

Hermite Polynomials
-------------------

.. function:: double gsl_sf_hermite_prob (const int n, const double x)
              int gsl_sf_hermite_prob_e (const int n, const double x, gsl_sf_result * result)

   These routines evaluate the probabilists' Hermite polynomial :math:`He_n(x)` of order :data:`n` at position :data:`x`.

.. function:: int gsl_sf_hermite_prob_array (const int nmax, const double x, double * result_array)

   This routine evaluates all probabilists' Hermite polynomials :math:`He_n(x)` up to order :data:`nmax` at position :data:`x`.
   The results are stored in :data:`result_array`.

.. function:: double gsl_sf_hermite_prob_series (const int n, const double x, const double * a)
              int gsl_sf_hermite_prob_series_e (const int n, const double x, const double * a, gsl_sf_result * result)

   These routines evaluate the series :math:`\sum_{j=0}^n a_j He_j(x)` with :math:`He_j` being the
   :math:`j`-th probabilists' Hermite polynomial using the Clenshaw algorithm.

.. function:: double gsl_sf_hermite_phys (const int n, const double x)
              int gsl_sf_hermite_phys_e (const int n, const double x, gsl_sf_result * result)

   These routines evaluate the physicists' Hermite polynomial :math:`H_n(x)` of order :data:`n` at position :data:`x`.

.. function:: int gsl_sf_hermite_phys_array (const int nmax, const double x, double * result_array)

   This routine evaluates all physicists' Hermite polynomials :math:`H_n` up to order :data:`nmax` at position :data:`x`.
   The results are stored in :data:`result_array`.

.. function:: double gsl_sf_hermite_phys_series (const int n, const double x, const double * a)
              int gsl_sf_hermite_phys_series_e (const int n, const double x, const double * a, gsl_sf_result * result)

   These routines evaluate the series :math:`\sum_{j=0}^n a_j H_j(x)` with :math:`H_j` being
   the :math:`j`-th physicists' Hermite polynomial using the Clenshaw algorithm.

Hermite Functions
-----------------

.. function:: double gsl_sf_hermite_func (const int n, const double x)
              int gsl_sf_hermite_func_e (const int n, const double x, gsl_sf_result * result)

   These routines evaluate the Hermite function :math:`\psi_n(x)` of order :data:`n` at position :data:`x`.

.. function:: int gsl_sf_hermite_func_array (const int nmax, const double x, double * result_array)

   This routine evaluates all Hermite functions :math:`\psi_n(x)` up to order :data:`nmax` at position :data:`x`.
   The results are stored in :data:`result_array`.

.. function:: double gsl_sf_hermite_func_series (const int n, const double x, const double * a)
              int gsl_sf_hermite_func_series_e (const int n, const double x, const double * a, gsl_sf_result * result)

   These routines evaluate the series :math:`\sum_{j=0}^n a_j \psi_j(x)` with :math:`\psi_j` being
   the :math:`j`-th Hermite function using the Clenshaw algorithm.

Derivatives of Hermite Polynomials
----------------------------------
.. index::
   single: Hermite polynomials, derivatives

.. function:: double gsl_sf_hermite_prob_der (const int m, const  int n, const double x)
              int gsl_sf_hermite_prob_der_e (const int m, const  int n, const double x, gsl_sf_result * result)

   These routines evaluate the :data:`m`-th derivative of the probabilists' Hermite polynomial :math:`He_n(x)`
   of order :data:`n` at position :data:`x`.

.. function:: int gsl_sf_hermite_prob_array_der (const int m, const int nmax, const double x, double * result_array)

   This routine evaluates the :data:`m`-th derivative of all probabilists' Hermite polynomials :math:`He_n(x)` up to
   order :data:`nmax` at position :data:`x`. The results are stored in :data:`result_array`.

.. function:: int gsl_sf_hermite_prob_der_array (const int mmax, const int n, const double x, double * result_array)

   This routine evaluates all derivatives (starting from 0) up to the :data:`mmax`-th derivative of the probabilists' Hermite
   polynomial of order :data:`n` :math:`He_n(x)` at position :data:`x`. The results are stored in :data:`result_array`.

.. function:: double gsl_sf_hermite_phys_der (const int m, const int n, const double x)
              int gsl_sf_hermite_phys_der_e (const int m, const int n, const double x, gsl_sf_result * result)

   These routines evaluate the :data:`m`-th derivative of the physicists' Hermite polynomial :math:`H_n(x)` of order :data:`n` at position :data:`x`.

.. function::  int gsl_sf_hermite_phys_array_der (const int m, const int nmax, const double x, double * result_array)

   This routine evaluates the :data:`m`-th derivative of all physicists' Hermite polynomials :math:`H_n` up to order :data:`nmax` at position :data:`x`.
   The results are stored in :data:`result_array`.

.. function:: int gsl_sf_hermite_phys_der_array (const int mmax, const int n, const double x, double * result_array)

   This routine evaluates all derivatives (starting from 0) up to the :data:`mmax`-th derivative of the
   physicists' Hermite polynomial of order :data:`n` :math:`H_n` at position :data:`x`. The results are stored in :data:`result_array`.

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

These routines calculate the :math:`s`-th zero of the Hermite Polynomial/Function of order
:math:`n`. Since the zeros are symmetrical around zero, only positive zeros are calculated,
ordered from smallest to largest, starting from index 1. Only for odd polynomial orders a
zeroth zero exists, its value always being zero.

.. function:: double gsl_sf_hermite_prob_zero (const int n, const int s)
              int gsl_sf_hermite_prob_zero_e (const int n, const int s, gsl_sf_result * result)

   These routines evaluate the :data:`s`-th zero of the probabilists' Hermite polynomial :math:`He_n(x)` of order :data:`n`.

.. function:: double gsl_sf_hermite_phys_zero (const int n, const int s)
              int gsl_sf_hermite_phys_zero_e (const int n, const int s, gsl_sf_result * result)

   These routines evaluate the :data:`s`-th zero of the physicists' Hermite polynomial :math:`H_n(x)` of order :data:`n`.

.. function:: double gsl_sf_hermite_func_zero (const int n, const int s)
              int gsl_sf_hermite_func_zero_e (const int n, const int s, gsl_sf_result * result)

   These routines evaluate the :data:`s`-th zero of the Hermite function :math:`\psi_n(x)` of order :data:`n`.
