.. index:: Fermi-Dirac function

The functions described in this section are declared in the header file
:file:`gsl_sf_fermi_dirac.h`.

Complete Fermi-Dirac Integrals
------------------------------
.. index::
   single: complete Fermi-Dirac integrals
   single: Fj(x), Fermi-Dirac integral

The complete Fermi-Dirac integral :math:`F_j(x)` is given by,

.. only:: not texinfo

   .. math:: F_j(x) := {1\over\Gamma(j+1)} \int_0^\infty dt {t^j  \over (\exp(t-x) + 1)}

.. only:: texinfo

   ::

      F_j(x) := (1/\Gamma(j+1)) \int_0^\infty dt (t^j / (\exp(t-x) + 1))

Note that the Fermi-Dirac integral is sometimes defined without the
normalisation factor in other texts.

.. function:: double gsl_sf_fermi_dirac_m1 (double x)
              int gsl_sf_fermi_dirac_m1_e (double x, gsl_sf_result * result)

   These routines compute the complete Fermi-Dirac integral with an index of :math:`-1`. 
   This integral is given by 
   :math:`F_{-1}(x) = e^x / (1 + e^x)`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_fermi_dirac_0 (double x)
              int gsl_sf_fermi_dirac_0_e (double x, gsl_sf_result * result)

   These routines compute the complete Fermi-Dirac integral with an index of :math:`0`. 
   This integral is given by :math:`F_0(x) = \ln(1 + e^x)`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_fermi_dirac_1 (double x)
              int gsl_sf_fermi_dirac_1_e (double x, gsl_sf_result * result)

   These routines compute the complete Fermi-Dirac integral with an index of :math:`1`,
   :math:`F_1(x) = \int_0^\infty dt (t /(\exp(t-x)+1))`.
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: double gsl_sf_fermi_dirac_2 (double x)
              int gsl_sf_fermi_dirac_2_e (double x, gsl_sf_result * result)

   These routines compute the complete Fermi-Dirac integral with an index
   of :math:`2`,
   :math:`F_2(x) = (1/2) \int_0^\infty dt (t^2 /(\exp(t-x)+1))`.
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: double gsl_sf_fermi_dirac_int (int j, double x)
              int gsl_sf_fermi_dirac_int_e (int j, double x, gsl_sf_result * result)

   These routines compute the complete Fermi-Dirac integral with an integer
   index of :math:`j`,
   :math:`F_j(x) = (1/\Gamma(j+1)) \int_0^\infty dt (t^j /(\exp(t-x)+1))`.
.. Complete integral F_j(x) for integer j
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: double gsl_sf_fermi_dirac_mhalf (double x)
              int gsl_sf_fermi_dirac_mhalf_e (double x, gsl_sf_result * result)

   These routines compute the complete Fermi-Dirac integral 
   :math:`F_{-1/2}(x)`.
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: double gsl_sf_fermi_dirac_half (double x)
              int gsl_sf_fermi_dirac_half_e (double x, gsl_sf_result * result)

   These routines compute the complete Fermi-Dirac integral 
   :math:`F_{1/2}(x)`.
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: double gsl_sf_fermi_dirac_3half (double x)
              int gsl_sf_fermi_dirac_3half_e (double x, gsl_sf_result * result)

   These routines compute the complete Fermi-Dirac integral 
   :math:`F_{3/2}(x)`.
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EOVRFLW

Incomplete Fermi-Dirac Integrals
--------------------------------
.. index::
   single: incomplete Fermi-Dirac integral
   single:  Fj(x,b), incomplete Fermi-Dirac integral

The incomplete Fermi-Dirac integral :math:`F_j(x,b)` is given by,

.. only:: not texinfo

   .. math:: F_j(x,b) := {1\over\Gamma(j+1)} \int_b^\infty dt {t^j  \over (\exp(t-x) + 1)}

.. only:: texinfo

   ::

      F_j(x,b) := (1/\Gamma(j+1)) \int_b^\infty dt (t^j / (\Exp(t-x) + 1))

.. function:: double gsl_sf_fermi_dirac_inc_0 (double x, double b)
              int gsl_sf_fermi_dirac_inc_0_e (double x, double b, gsl_sf_result * result)

   These routines compute the incomplete Fermi-Dirac integral with an index
   of zero,
   :math:`F_0(x,b) = \ln(1 + e^{b-x}) - (b-x)`
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EDOM
