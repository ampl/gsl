.. index::
   single: W function
   single: Lambert function

Lambert's W functions, :math:`W(x)`, are defined to be solutions
of the equation :math:`W(x) \exp(W(x)) = x`. This function has
multiple branches for :math:`x < 0`; however, it has only
two real-valued branches. We define :math:`W_0(x)` to be the
principal branch, where :math:`W > -1` for :math:`x < 0`, and 
:math:`W_{-1}(x)`
to be the other real branch, where
:math:`W < -1` for :math:`x < 0`.  The Lambert functions are
declared in the header file :file:`gsl_sf_lambert.h`.

.. function:: double gsl_sf_lambert_W0 (double x)
              int gsl_sf_lambert_W0_e (double x, gsl_sf_result * result)

   These compute the principal branch of the Lambert W function, :math:`W_0(x)`.
.. exceptions: GSL_EDOM, GSL_EMAXITER

.. function:: double gsl_sf_lambert_Wm1 (double x)
              int gsl_sf_lambert_Wm1_e (double x, gsl_sf_result * result)

   These compute the secondary real-valued branch of the Lambert W function, 
   :math:`W_{-1}(x)`.
.. exceptions: GSL_EDOM, GSL_EMAXITER
