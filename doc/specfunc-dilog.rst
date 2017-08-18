.. index:: dilogarithm

The dilogarithm is defined as

.. only:: not texinfo

   .. math:: Li_2(z) = - \int_0^z ds {\log{(1-s)} \over s}

.. only:: texinfo

   .. math:: Li_2(z) = - \int_0^z ds log(1-s) / s

The functions described in this section are declared in the header file
:file:`gsl_sf_dilog.h`.

Real Argument
-------------

.. function:: double gsl_sf_dilog (double x)
              int gsl_sf_dilog_e (double x, gsl_sf_result * result)

   These routines compute the dilogarithm for a real argument. In Lewin's
   notation this is :math:`Li_2(x)`, the real part of the dilogarithm of a
   real :math:`x`.  It is defined by the integral representation

   .. math:: Li_2(x) = - \Re \int_0^x ds \log(1-s) / s

   Note that :math:`\Im(Li_2(x)) = 0` for
   :math:`x \le 1`, and :math:`-\pi\log(x)` for :math:`x > 1`.

   Note that Abramowitz & Stegun refer to the Spence integral
   :math:`S(x) = Li_2(1 - x)` as the dilogarithm rather than :math:`Li_2(x)`.

Complex Argument
----------------

.. function:: int gsl_sf_complex_dilog_e (double r, double theta, gsl_sf_result * result_re, gsl_sf_result * result_im)

   This function computes the full complex-valued dilogarithm for the
   complex argument :math:`z = r \exp(i \theta)`. The real and imaginary
   parts of the result are returned in :data:`result_re`, :data:`result_im`.
