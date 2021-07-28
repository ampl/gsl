.. index:: trigonometric functions

The library includes its own trigonometric functions in order to provide
consistency across platforms and reliable error estimates.  These
functions are declared in the header file :file:`gsl_sf_trig.h`.

Circular Trigonometric Functions
--------------------------------

.. index::
   single: sine function, special functions

.. function:: double gsl_sf_sin (double x)
              int gsl_sf_sin_e (double x, gsl_sf_result * result)

   These routines compute the sine function :math:`\sin(x)`.
.. Exceptional Return Values:

.. index::
   single: cosine function, special functions

.. function:: double gsl_sf_cos (double x)
              int gsl_sf_cos_e (double x, gsl_sf_result * result)

   These routines compute the cosine function :math:`\cos(x)`.
.. Exceptional Return Values:

.. index::
   single: hypot function, special functions

.. function:: double gsl_sf_hypot (double x, double y)
              int gsl_sf_hypot_e (double x, double y, gsl_sf_result * result)

   These routines compute the hypotenuse function :math:`\sqrt{x^2 + y^2}`
   avoiding overflow and underflow.
.. Exceptional Return Values:

.. index::
   single: complex sinc function, special functions

.. function:: double gsl_sf_sinc (double x)
              int gsl_sf_sinc_e (double x, gsl_sf_result * result)

   These routines compute :math:`\sinc(x) = \sin(\pi x) / (\pi x)` for any
   value of :data:`x`.
.. Exceptional Return Values: none

Trigonometric Functions for Complex Arguments
---------------------------------------------

.. index::
   single: complex sine function, special functions

.. function:: int gsl_sf_complex_sin_e (double zr, double zi, gsl_sf_result * szr, gsl_sf_result * szi)

   This function computes the complex sine, :math:`\sin(z_r + i z_i)` storing
   the real and imaginary parts in :data:`szr`, :data:`szi`.
.. Exceptional Return Values: GSL_EOVRFLW

.. index::
   single: complex cosine function, special functions

.. function:: int gsl_sf_complex_cos_e (double zr, double zi, gsl_sf_result * czr, gsl_sf_result * czi)

   This function computes the complex cosine, :math:`\cos(z_r + i z_i)` storing
   the real and imaginary parts in :data:`czr`, :data:`czi`.
.. Exceptional Return Values: GSL_EOVRFLW

.. index::
   single: complex log sine function, special functions

.. function:: int gsl_sf_complex_logsin_e (double zr, double zi, gsl_sf_result * lszr, gsl_sf_result * lszi)

   This function computes the logarithm of the complex sine,
   :math:`\log(\sin(z_r + i z_i))` storing the real and imaginary parts in
   :data:`lszr`, :data:`lszi`.
.. Exceptional Return Values: GSL_EDOM, GSL_ELOSS

Hyperbolic Trigonometric Functions
----------------------------------

.. index:: logarithm of sinh function, special functions

.. function:: double gsl_sf_lnsinh (double x)
              int gsl_sf_lnsinh_e (double x, gsl_sf_result * result)

   These routines compute :math:`\log(\sinh(x))` for :math:`x > 0`.
.. Domain: x > 0 
.. Exceptional Return Values: GSL_EDOM

.. index::
   single: logarithm of cosh function, special functions

.. function:: double gsl_sf_lncosh (double x)
              int gsl_sf_lncosh_e (double x, gsl_sf_result * result)

   These routines compute :math:`\log(\cosh(x))` for any :data:`x`.
.. Exceptional Return Values: none

Conversion Functions
--------------------
.. index::
   single: polar to rectangular conversion
   single: rectangular to polar conversion

.. function:: int gsl_sf_polar_to_rect (double r, double theta, gsl_sf_result * x, gsl_sf_result * y)

   This function converts the polar coordinates (:data:`r`, :data:`theta`) to
   rectilinear coordinates (:data:`x`, :data:`y`), :math:`x = r\cos(\theta)`,
   :math:`y = r\sin(\theta)`.
.. Exceptional Return Values: GSL_ELOSS

.. function:: int gsl_sf_rect_to_polar (double x, double y, gsl_sf_result * r, gsl_sf_result * theta)

   This function converts the rectilinear coordinates (:data:`x`, :data:`y`) to
   polar coordinates (:data:`r`, :data:`theta`), such that :math:`x = r\cos(\theta)`,
   :math:`y = r\sin(\theta)`.  The argument :data:`theta`
   lies in the range :math:`[-\pi, \pi]`.
.. Exceptional Return Values: GSL_EDOM

Restriction Functions
---------------------
.. index::
   single: angular reduction
   single: reduction of angular variables

.. function:: double gsl_sf_angle_restrict_symm (double theta)
              int gsl_sf_angle_restrict_symm_e (double * theta)

   These routines force the angle :data:`theta` to lie in the range
   :math:`(-\pi,\pi]`.  

   Note that the mathematical value of :math:`\pi` is slightly greater
   than :macro:`M_PI`, so the machine numbers :macro:`M_PI` and :code:`-M_PI`
   are included in the range.
.. Exceptional Return Values: GSL_ELOSS

.. function:: double gsl_sf_angle_restrict_pos (double theta)
              int gsl_sf_angle_restrict_pos_e (double * theta)

   These routines force the angle :data:`theta` to lie in the range :math:`[0, 2\pi)`.

   Note that the mathematical value of :math:`2\pi` is slightly greater
   than :code:`2*M_PI`, so the machine number :code:`2*M_PI` is included in
   the range.

.. Exceptional Return Values: GSL_ELOSS

Trigonometric Functions With Error Estimates
--------------------------------------------

.. function:: int gsl_sf_sin_err_e (double x, double dx, gsl_sf_result * result)

   This routine computes the sine of an angle :data:`x` with an associated 
   absolute error :data:`dx`,
   :math:`\sin(x \pm dx)`.  Note that this function is provided in the error-handling form only since
   its purpose is to compute the propagated error.

.. function:: int gsl_sf_cos_err_e (double x, double dx, gsl_sf_result * result)

   This routine computes the cosine of an angle :data:`x` with an associated
   absolute error :data:`dx`, 
   :math:`\cos(x \pm dx)`.  Note that this function is provided in the error-handling form only since
   its purpose is to compute the propagated error.
