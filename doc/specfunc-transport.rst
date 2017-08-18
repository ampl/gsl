.. index:: transport functions

The transport functions :math:`J(n,x)` are defined by the integral 
representations

.. math:: J(n,x) = \int_0^x t^n e^t /(e^t - 1)^2 dt

They are declared in the header file :file:`gsl_sf_transport.h`.

.. function:: double gsl_sf_transport_2 (double x)
              int gsl_sf_transport_2_e (double x, gsl_sf_result * result)

   These routines compute the transport function :math:`J(2,x)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_transport_3 (double x)
              int gsl_sf_transport_3_e (double x, gsl_sf_result * result)

   These routines compute the transport function :math:`J(3,x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_transport_4 (double x)
              int gsl_sf_transport_4_e (double x, gsl_sf_result * result)

   These routines compute the transport function :math:`J(4,x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_transport_5 (double x)
              int gsl_sf_transport_5_e (double x, gsl_sf_result * result)

   These routines compute the transport function :math:`J(5,x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW
