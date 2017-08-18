The following routines compute the gamma and beta functions in their
full and incomplete forms, as well as various kinds of factorials.
The functions described in this section are declared in the header
file :file:`gsl_sf_gamma.h`.

Gamma Functions
---------------
.. index:: gamma functions

The Gamma function is defined by the following integral,

.. math:: \Gamma(x) = \int_0^{\infty} dt t^{x-1} \exp(-t)

It is related to the factorial function by :math:`\Gamma(n) = (n-1)!`
for positive integer :math:`n`.  Further information on the Gamma function
can be found in Abramowitz & Stegun, Chapter 6.  

.. function:: double gsl_sf_gamma (double x)
              int gsl_sf_gamma_e (double x, gsl_sf_result * result)

   These routines compute the Gamma function :math:`\Gamma(x)`, subject to :math:`x`
   not being a negative integer or zero.  The function is computed using the real
   Lanczos method. The maximum value of :math:`x` such that :math:`\Gamma(x)` is not
   considered an overflow is given by the macro :macro:`GSL_SF_GAMMA_XMAX`
   and is 171.0.
.. exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND

.. index:: logarithm of Gamma function

.. function:: double gsl_sf_lngamma (double x)
              int gsl_sf_lngamma_e (double x, gsl_sf_result * result)

   These routines compute the logarithm of the Gamma function,
   :math:`\log(\Gamma(x))`, subject to :math:`x` not being a negative
   integer or zero.  For :math:`x < 0` the real part of :math:`\log(\Gamma(x))` is
   returned, which is equivalent to :math:`\log(|\Gamma(x)|)`.  The function
   is computed using the real Lanczos method.
.. exceptions: GSL_EDOM, GSL_EROUND

.. function:: int gsl_sf_lngamma_sgn_e (double x, gsl_sf_result * result_lg, double * sgn)

   This routine computes the sign of the gamma function and the logarithm of
   its magnitude, subject to :math:`x` not being a negative integer or zero.  The
   function is computed using the real Lanczos method.  The value of the
   gamma function and its error can be reconstructed using the relation 
   :math:`\Gamma(x) = sgn * \exp(result\_lg)`, taking into account the two 
   components of :data:`result_lg`.
.. exceptions: GSL_EDOM, GSL_EROUND

.. index:: Regulated Gamma function

.. function:: double gsl_sf_gammastar (double x)
              int gsl_sf_gammastar_e (double x, gsl_sf_result * result)

   These routines compute the regulated Gamma Function :math:`\Gamma^*(x)`
   for :math:`x > 0`. The regulated gamma function is given by,

   .. only:: not texinfo

      .. math::

         \Gamma^*(x) &= \Gamma(x)/(\sqrt{2\pi} x^{(x-1/2)} \exp(-x))\cr
                     &= \left(1 + {1 \over 12x} + ...\right) \quad\hbox{for~} x\to \infty\cr

   .. only:: texinfo

      ::

         \Gamma^*(x) = \Gamma(x)/(\sqrt{2\pi} x^{(x-1/2)} \exp(-x))
                     = (1 + (1/12x) + ...)  for x \to \infty

   and is a useful suggestion of Temme.
.. exceptions: GSL_EDOM

.. index:: Reciprocal Gamma function

.. function:: double gsl_sf_gammainv (double x)
              int gsl_sf_gammainv_e (double x, gsl_sf_result * result)

   These routines compute the reciprocal of the gamma function,
   :math:`1/\Gamma(x)` using the real Lanczos method.
.. exceptions: GSL_EUNDRFLW, GSL_EROUND

.. index:: Complex Gamma function

.. function:: int gsl_sf_lngamma_complex_e (double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg)

   This routine computes :math:`\log(\Gamma(z))` for complex :math:`z = z_r + i z_i`
   and :math:`z` not a negative integer or zero, using the complex Lanczos
   method.  The returned parameters are :math:`lnr = \log|\Gamma(z)|` and
   :math:`arg = \arg(\Gamma(z))` in :math:`(-\pi,\pi]`.  Note that the phase
   part (:data:`arg`) is not well-determined when :math:`|z|` is very large,
   due to inevitable roundoff in restricting to :math:`(-\pi,\pi]`.  This
   will result in a :macro:`GSL_ELOSS` error when it occurs.  The absolute
   value part (:data:`lnr`), however, never suffers from loss of precision.
.. exceptions: GSL_EDOM, GSL_ELOSS

Factorials
----------
.. index:: factorial

Although factorials can be computed from the Gamma function, using
the relation :math:`n! = \Gamma(n+1)` for non-negative integer :math:`n`,
it is usually more efficient to call the functions in this section,
particularly for small values of :math:`n`, whose factorial values are
maintained in hardcoded tables.

.. index:: factorial

.. function:: double gsl_sf_fact (unsigned int n)
              int gsl_sf_fact_e (unsigned int n, gsl_sf_result * result)

   These routines compute the factorial :math:`n!`.  The factorial is
   related to the Gamma function by :math:`n! = \Gamma(n+1)`.
   The maximum value of :math:`n` such that :math:`n!` is not
   considered an overflow is given by the macro :macro:`GSL_SF_FACT_NMAX`
   and is 170.
.. exceptions: GSL_EDOM, GSL_EOVRFLW

.. index:: double factorial

.. function:: double gsl_sf_doublefact (unsigned int n)
              int gsl_sf_doublefact_e (unsigned int n, gsl_sf_result * result)

   These routines compute the double factorial :math:`n!! = n(n-2)(n-4) \dots`. 
   The maximum value of :math:`n` such that :math:`n!!` is not
   considered an overflow is given by the macro :macro:`GSL_SF_DOUBLEFACT_NMAX`
   and is 297.
.. exceptions: GSL_EDOM, GSL_EOVRFLW

.. index:: logarithm of factorial

.. function:: double gsl_sf_lnfact (unsigned int n)
              int gsl_sf_lnfact_e (unsigned int n, gsl_sf_result * result)

   These routines compute the logarithm of the factorial of :data:`n`,
   :math:`\log(n!)`.  The algorithm is faster than computing
   :math:`\ln(\Gamma(n+1))` via :func:`gsl_sf_lngamma` for :math:`n < 170`,
   but defers for larger :data:`n`.
.. exceptions: none

.. index:: logarithm of double factorial

.. function:: double gsl_sf_lndoublefact (unsigned int n)
              int gsl_sf_lndoublefact_e (unsigned int n, gsl_sf_result * result)

   These routines compute the logarithm of the double factorial of :data:`n`,
   :math:`\log(n!!)`.
.. exceptions: none

.. index:: combinatorial factor C(m,n)

.. function:: double gsl_sf_choose (unsigned int n, unsigned int m)
              int gsl_sf_choose_e (unsigned int n, unsigned int m, gsl_sf_result * result)

   These routines compute the combinatorial factor :code:`n choose m`
   :math:`= n!/(m!(n-m)!)`
.. exceptions: GSL_EDOM, GSL_EOVRFLW

.. index:: logarithm of combinatorial factor C(m,n)

.. function:: double gsl_sf_lnchoose (unsigned int n, unsigned int m)
              int gsl_sf_lnchoose_e (unsigned int n, unsigned int m, gsl_sf_result * result)

   These routines compute the logarithm of :code:`n choose m`.  This is
   equivalent to the sum :math:`\log(n!) - \log(m!) - \log((n-m)!)`.
.. exceptions: GSL_EDOM 

.. index::
   single: Taylor coefficients, computation of

.. function:: double gsl_sf_taylorcoeff (int n, double x)
              int gsl_sf_taylorcoeff_e (int n, double x, gsl_sf_result * result)

   These routines compute the Taylor coefficient :math:`x^n / n!` for 
   :math:`x \ge 0`, :math:`n \ge 0`
.. exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. _pochhammer-symbol:

Pochhammer Symbol
-----------------

.. index::
   single: Pochhammer symbol
   single:  Apell symbol, see Pochhammer symbol

.. function:: double gsl_sf_poch (double a, double x)
              int gsl_sf_poch_e (double a, double x, gsl_sf_result * result)

   These routines compute the Pochhammer symbol :math:`(a)_x = \Gamma(a + x)/\Gamma(a)`.
   The Pochhammer symbol is also known as the Apell symbol and
   sometimes written as :math:`(a,x)`.  When :math:`a` and :math:`a + x` 
   are negative integers or zero, the limiting value of the ratio is returned. 
.. exceptions:  GSL_EDOM, GSL_EOVRFLW

.. index:: logarithm of Pochhammer symbol

.. function:: double gsl_sf_lnpoch (double a, double x)
              int gsl_sf_lnpoch_e (double a, double x, gsl_sf_result * result)

   These routines compute the logarithm of the Pochhammer symbol,
   :math:`\log((a)_x) = \log(\Gamma(a + x)/\Gamma(a))`.
.. exceptions:  GSL_EDOM

.. function:: int gsl_sf_lnpoch_sgn_e (double a, double x, gsl_sf_result * result, double * sgn)

   These routines compute the sign of the Pochhammer symbol and the
   logarithm of its magnitude.  The computed parameters are :math:`result = \log(|(a)_x|)`
   with a corresponding error term,  and :math:`sgn = \sgn((a)_x)` where :math:`(a)_x = \Gamma(a + x)/\Gamma(a)`.
.. exceptions:  GSL_EDOM

.. index:: relative Pochhammer symbol

.. function:: double gsl_sf_pochrel (double a, double x)
              int gsl_sf_pochrel_e (double a, double x, gsl_sf_result * result)

   These routines compute the relative Pochhammer symbol :math:`((a)_x - 1)/x`
   where :math:`(a)_x = \Gamma(a + x)/\Gamma(a)`.
.. exceptions:  GSL_EDOM

Incomplete Gamma Functions
--------------------------

.. index::
   single: non-normalized incomplete Gamma function
   single: unnormalized incomplete Gamma function

.. function:: double gsl_sf_gamma_inc (double a, double x)
              int gsl_sf_gamma_inc_e (double a, double x, gsl_sf_result * result)

   These functions compute the unnormalized incomplete Gamma Function
   :math:`\Gamma(a,x) = \int_x^\infty dt t^{(a-1)} \exp(-t)`
   for :math:`a` real and :math:`x \ge 0`.
.. exceptions: GSL_EDOM

.. index:: incomplete Gamma function

.. function:: double gsl_sf_gamma_inc_Q (double a, double x)
              int gsl_sf_gamma_inc_Q_e (double a, double x, gsl_sf_result * result)

   These routines compute the normalized incomplete Gamma Function
   :math:`Q(a,x) = 1/\Gamma(a) \int_x^\infty dt t^{(a-1)} \exp(-t)`
   for :math:`a > 0`, :math:`x \ge 0`.
.. exceptions: GSL_EDOM

.. index:: complementary incomplete Gamma function

.. function:: double gsl_sf_gamma_inc_P (double a, double x)
              int gsl_sf_gamma_inc_P_e (double a, double x, gsl_sf_result * result)

   These routines compute the complementary normalized incomplete Gamma Function
   :math:`P(a,x) = 1 - Q(a,x) = 1/\Gamma(a) \int_0^x dt t^{(a-1)} \exp(-t)`
   for :math:`a > 0`, :math:`x \ge 0`.

   Note that Abramowitz & Stegun call :math:`P(a,x)` the incomplete gamma
   function (section 6.5).
.. exceptions: GSL_EDOM

Beta Functions
--------------

.. index:: Beta function

.. function:: double gsl_sf_beta (double a, double b)
              int gsl_sf_beta_e (double a, double b, gsl_sf_result * result)

   These routines compute the Beta Function, :math:`B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b)`
   subject to :math:`a` and :math:`b` not being negative integers.
.. exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. index:: logarithm of Beta function

.. function:: double gsl_sf_lnbeta (double a, double b)
              int gsl_sf_lnbeta_e (double a, double b, gsl_sf_result * result)

   These routines compute the logarithm of the Beta Function, :math:`\log(B(a,b))`
   subject to :math:`a` and :math:`b` not being negative integers.
.. exceptions: GSL_EDOM

Incomplete Beta Function
------------------------

.. index::
   single: incomplete Beta function, normalized
   single: normalized incomplete Beta function
   single: Beta function, incomplete normalized 

.. function:: double gsl_sf_beta_inc (double a, double b, double x)
              int gsl_sf_beta_inc_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the normalized incomplete Beta function
   :math:`I_x(a,b) = B_x(a,b) / B(a,b)` where
   
   .. math:: B_x(a,b) = \int_0^x t^{a-1} (1-t)^{b-1} dt

   for :math:`0 \le x \le 1`.
   For :math:`a > 0`, :math:`b > 0` the value is computed using
   a continued fraction expansion.  For all other values it is computed using 
   the relation
   
   .. only:: not texinfo

      .. math:: I_x(a,b,x) = (1/a) x^a {}_2F_1(a,1-b,a+1,x)/B(a,b)

   .. only:: texinfo

      ::

         I_x(a,b,x) = (1/a) x^a 2F1(a,1-b,a+1,x) / B(a,b)
