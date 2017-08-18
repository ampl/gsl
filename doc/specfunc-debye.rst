.. index:: Debye functions

The Debye functions :math:`D_n(x)` are defined by the following integral,

.. only:: not texinfo

   .. math:: D_n(x) = {n \over x^n} \int_0^x dt {t^n \over e^t - 1}

.. only:: texinfo

   .. math:: D_n(x) = n/x^n \int_0^x dt (t^n/(e^t - 1))

For further information see Abramowitz &
Stegun, Section 27.1.  The Debye functions are declared in the header
file :file:`gsl_sf_debye.h`.

.. function:: double gsl_sf_debye_1 (double x)
              int gsl_sf_debye_1_e (double x, gsl_sf_result * result)

   These routines compute the first-order Debye function :math:`D_1(x)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_debye_2 (double x)
              int gsl_sf_debye_2_e (double x, gsl_sf_result * result)

   These routines compute the second-order Debye function :math:`D_2(x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_debye_3 (double x)
              int gsl_sf_debye_3_e (double x, gsl_sf_result * result)

   These routines compute the third-order Debye function :math:`D_3(x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_debye_4 (double x)
              int gsl_sf_debye_4_e (double x, gsl_sf_result * result)

   These routines compute the fourth-order Debye function :math:`D_4(x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_debye_5 (double x)
              int gsl_sf_debye_5_e (double x, gsl_sf_result * result)

   These routines compute the fifth-order Debye function :math:`D_5(x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_debye_6 (double x)
              int gsl_sf_debye_6_e (double x, gsl_sf_result * result)

   These routines compute the sixth-order Debye function :math:`D_6(x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW
