.. index:: synchrotron functions

The functions described in this section are declared in the header file
:file:`gsl_sf_synchrotron.h`.

.. function:: double gsl_sf_synchrotron_1 (double x)
              int gsl_sf_synchrotron_1_e (double x, gsl_sf_result * result)

   These routines compute the first synchrotron function 
   :math:`x \int_x^\infty dt K_{5/3}(t)`
   for :math:`x \ge 0`.
.. Domain: x >= 0.0
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW
 
.. function:: double gsl_sf_synchrotron_2 (double x)
              int gsl_sf_synchrotron_2_e (double x, gsl_sf_result * result)

   These routines compute the second synchrotron function 
   :math:`x K_{2/3}(x)` for :math:`x \ge 0`.
.. Domain: x >= 0.0
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW
