.. index::
   single: exponential function
   single: exp

The functions described in this section are declared in the header file
:file:`gsl_sf_exp.h`.

Exponential Function
--------------------

.. function:: double gsl_sf_exp (double x)
              int gsl_sf_exp_e (double x, gsl_sf_result * result)

   These routines provide an exponential function :math:`\exp(x)` using GSL
   semantics and error checking.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: int gsl_sf_exp_e10_e (double x, gsl_sf_result_e10 * result)

   This function computes the exponential :math:`\exp(x)` using the
   :type:`gsl_sf_result_e10` type to return a result with extended range.
   This function may be useful if the value of :math:`\exp(x)` would
   overflow the  numeric range of :code:`double`.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_exp_mult (double x, double y)
              int gsl_sf_exp_mult_e (double x, double y, gsl_sf_result * result)

   These routines exponentiate :data:`x` and multiply by the factor :data:`y`
   to return the product :math:`y \exp(x)`.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: int gsl_sf_exp_mult_e10_e (const double x, const double y, gsl_sf_result_e10 * result)

   This function computes the product :math:`y \exp(x)` using the
   :type:`gsl_sf_result_e10` type to return a result with extended numeric
   range.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

Relative Exponential Functions
------------------------------

.. function:: double gsl_sf_expm1 (double x)
              int gsl_sf_expm1_e (double x, gsl_sf_result * result)

   These routines compute the quantity :math:`\exp(x)-1` using an algorithm
   that is accurate for small :math:`x`.
.. Exceptional Return Values:  GSL_EOVRFLW

.. function:: double gsl_sf_exprel (double x)
              int gsl_sf_exprel_e (double x, gsl_sf_result * result)

   These routines compute the quantity :math:`(\exp(x)-1)/x` using an
   algorithm that is accurate for small :data:`x`.  For small :data:`x` the
   algorithm is based on the expansion
   :math:`(\exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + \dots`.
.. Exceptional Return Values:  GSL_EOVRFLW

.. function:: double gsl_sf_exprel_2 (double x)
              int gsl_sf_exprel_2_e (double x, gsl_sf_result * result)

   These routines compute the quantity :math:`2(\exp(x)-1-x)/x^2` using an
   algorithm that is accurate for small :data:`x`.  For small :data:`x` the
   algorithm is based on the expansion
   :math:`2(\exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + \dots`.
.. Exceptional Return Values:  GSL_EOVRFLW

.. function:: double gsl_sf_exprel_n (int n, double x)
              int gsl_sf_exprel_n_e (int n, double x, gsl_sf_result * result)

   These routines compute the :math:`N`-relative exponential, which is the
   :data:`n`-th generalization of the functions :func:`gsl_sf_exprel` and
   :func:`gsl_sf_exprel_2`.  The :math:`N`-relative exponential is given by,

   .. only:: not texinfo

      .. math::

         \hbox{exprel}_N(x)
                     &= N!/x^N \left(\exp(x) - \sum_{k=0}^{N-1} x^k/k!\right)\cr
                     &= 1 + x/(N+1) + x^2/((N+1)(N+2)) + \dots\cr
                     &= {}_1F_1(1,1+N,x)\cr

   .. only:: texinfo

      ::

         exprel_N(x) = N!/x^N (\exp(x) - \sum_{k=0}^{N-1} x^k/k!)
                     = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
                     = 1F1 (1,1+N,x)
.. Exceptional Return Values: 

Exponentiation With Error Estimate
----------------------------------

.. function:: int gsl_sf_exp_err_e (double x, double dx, gsl_sf_result * result)

   This function exponentiates :data:`x` with an associated absolute error
   :data:`dx`.
.. Exceptional Return Values: 

.. function:: int gsl_sf_exp_err_e10_e (double x, double dx, gsl_sf_result_e10 * result)

   This function exponentiates a quantity :data:`x` with an associated absolute 
   error :data:`dx` using the :type:`gsl_sf_result_e10` type to return a result with
   extended range.
.. Exceptional Return Values: 

.. function:: int gsl_sf_exp_mult_err_e (double x, double dx, double y, double dy, gsl_sf_result * result)

   This routine computes the product :math:`y \exp(x)` for the quantities
   :data:`x`, :data:`y` with associated absolute errors :data:`dx`, :data:`dy`.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: int gsl_sf_exp_mult_err_e10_e (double x, double dx, double y, double dy, gsl_sf_result_e10 * result)

   This routine computes the product :math:`y \exp(x)` for the quantities
   :data:`x`, :data:`y` with associated absolute errors :data:`dx`, :data:`dy` using the
   :type:`gsl_sf_result_e10` type to return a result with extended range.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW
