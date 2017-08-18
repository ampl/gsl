.. index::
   single: Laguerre functions
   single: confluent hypergeometric function

The generalized Laguerre polynomials, sometimes referred to as
associated Laguerre polynomials, are defined in terms of confluent
hypergeometric functions as

.. only:: not texinfo

   .. math:: L^a_n(x) = {(a+1)_n \over n!} {}_1F_1(-n,a+1,x)

.. only:: texinfo

   L^a_n(x) = ((a+1)_n / n!) 1F1(-n,a+1,x)
   
where :math:`(a)_n` is the :ref:`Pochhammer symbol <pochhammer-symbol>` (rising factorial).
They are related to the plain
Laguerre polynomials :math:`L_n(x)` by :math:`L^0_n(x) = L_n(x)` and 
:math:`L^k_n(x) = (-1)^k (d^k/dx^k) L_{(n+k)}(x)`
For more information see Abramowitz & Stegun, Chapter 22.

The functions described in this section are
declared in the header file :file:`gsl_sf_laguerre.h`.

.. function:: double gsl_sf_laguerre_1 (double a, double x)
              double gsl_sf_laguerre_2 (double a, double x)
              double gsl_sf_laguerre_3 (double a, double x)
              int gsl_sf_laguerre_1_e (double a, double x, gsl_sf_result * result)
              int gsl_sf_laguerre_2_e (double a, double x, gsl_sf_result * result)
              int gsl_sf_laguerre_3_e (double a, double x, gsl_sf_result * result)

   These routines evaluate the generalized Laguerre polynomials
   :math:`L^a_1(x)`, :math:`L^a_2(x)`, :math:`L^a_3(x)` using explicit
   representations.
.. Exceptional Return Values: none

.. function:: double gsl_sf_laguerre_n (const int n, const double a, const double x)
              int gsl_sf_laguerre_n_e (int n, double a, double x, gsl_sf_result * result)

   These routines evaluate the generalized Laguerre polynomials
   :math:`L^a_n(x)` for :math:`a > -1`,
   :math:`n \ge 0`.
.. Domain: a > -1.0, n >= 0
.. Evaluate generalized Laguerre polynomials.
.. Exceptional Return Values: GSL_EDOM
