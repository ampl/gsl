.. index:: Dawson function

The Dawson integral is defined by

.. math:: \exp(-x^2) \int_0^x dt \exp(t^2)

A table of Dawson's integral can be found in Abramowitz &
Stegun, Table 7.5.  The Dawson functions are declared in the header file
:file:`gsl_sf_dawson.h`.

.. function:: double gsl_sf_dawson (double x)
              int gsl_sf_dawson_e (double x, gsl_sf_result * result)

   These routines compute the value of Dawson's integral for :data:`x`.
.. Exceptional Return Values: GSL_EUNDRFLW
