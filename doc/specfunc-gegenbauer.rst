.. index:: Gegenbauer functions

The Gegenbauer polynomials are defined in Abramowitz & Stegun, Chapter
22, where they are known as Ultraspherical polynomials.  The functions
described in this section are declared in the header file
:file:`gsl_sf_gegenbauer.h`.

.. function:: double gsl_sf_gegenpoly_1 (double lambda, double x)
              double gsl_sf_gegenpoly_2 (double lambda, double x)
              double gsl_sf_gegenpoly_3 (double lambda, double x)
              int gsl_sf_gegenpoly_1_e (double lambda, double x, gsl_sf_result * result)
              int gsl_sf_gegenpoly_2_e (double lambda, double x, gsl_sf_result * result)
              int gsl_sf_gegenpoly_3_e (double lambda, double x, gsl_sf_result * result)

   These functions evaluate the Gegenbauer polynomials
   :math:`C^{(\lambda)}_n(x)` using explicit
   representations for :math:`n = 1, 2, 3`.
.. Exceptional Return Values: none

.. function:: double gsl_sf_gegenpoly_n (int n, double lambda, double x)
              int gsl_sf_gegenpoly_n_e (int n, double lambda, double x, gsl_sf_result * result)

   These functions evaluate the Gegenbauer polynomial :math:`C^{(\lambda)}_n(x)`
   for a specific value of :data:`n`,
   :data:`lambda`, :data:`x` subject to :math:`\lambda > -1/2`, :math:`n \ge 0`.
.. Domain: lambda > -1/2, n >= 0
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_gegenpoly_array (int nmax, double lambda, double x, double result_array[])

   This function computes an array of Gegenbauer polynomials
   :math:`C^{(\lambda)}_n(x)`
   for :math:`n = 0, 1, 2, \dots, nmax`, subject
   to :math:`\lambda > -1/2`, :math:`nmax \ge 0`.
.. Conditions: n = 0, 1, 2, ... nmax
.. Domain: lambda > -1/2, nmax >= 0
.. Exceptional Return Values: GSL_EDOM
