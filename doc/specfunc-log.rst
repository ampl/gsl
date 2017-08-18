.. index:: logarithm and related functions

Information on the properties of the Logarithm function can be found in
Abramowitz & Stegun, Chapter 4.  The functions described in this section
are declared in the header file :file:`gsl_sf_log.h`.

.. function:: double gsl_sf_log (double x)
              int gsl_sf_log_e (double x, gsl_sf_result * result)

   These routines compute the logarithm of :data:`x`, :math:`\log(x)`, for
   :math:`x > 0`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_log_abs (double x)
              int gsl_sf_log_abs_e (double x, gsl_sf_result * result)

   These routines compute the logarithm of the magnitude of :data:`x`,
   :math:`\log(|x|)`, for :math:`x \ne 0`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_complex_log_e (double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * theta)

   This routine computes the complex logarithm of :math:`z = z_r + i z_i`.
   The results are returned as :data:`lnr`, :data:`theta` such that
   :math:`\exp(lnr + i \theta) = z_r + i z_i`, where :math:`\theta` lies in
   the range :math:`[-\pi,\pi]`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_log_1plusx (double x)
              int gsl_sf_log_1plusx_e (double x, gsl_sf_result * result)

   These routines compute :math:`\log(1 + x)` for :math:`x > -1` using an
   algorithm that is accurate for small :data:`x`.
.. Domain: x > -1.0
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_log_1plusx_mx (double x)
              int gsl_sf_log_1plusx_mx_e (double x, gsl_sf_result * result)

   These routines compute :math:`\log(1 + x) - x` for :math:`x > -1` using an
   algorithm that is accurate for small :data:`x`.
.. Domain: x > -1.0 
.. Exceptional Return Values: GSL_EDOM
