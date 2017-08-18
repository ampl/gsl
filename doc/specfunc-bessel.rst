.. index:: Bessel functions

The routines described in this section compute the Cylindrical Bessel
functions :math:`J_n(x)`, :math:`Y_n(x)`, Modified cylindrical Bessel
functions :math:`I_n(x)`, :math:`K_n(x)`, Spherical Bessel functions
:math:`j_l(x)`, :math:`y_l(x)`, and Modified Spherical Bessel functions
:math:`i_l(x)`, :math:`k_l(x)`.  For more information see Abramowitz & Stegun,
Chapters 9 and 10.  The Bessel functions are defined in the header file
:file:`gsl_sf_bessel.h`.

Regular Cylindrical Bessel Functions
------------------------------------
.. index:: Cylindrical Bessel Functions
.. index:: Regular Cylindrical Bessel Functions
.. index::
   single: J(x), Bessel Functions

.. function:: double gsl_sf_bessel_J0 (double x)
              int gsl_sf_bessel_J0_e (double x, gsl_sf_result * result)

   These routines compute the regular cylindrical Bessel function of zeroth
   order, :math:`J_0(x)`.

.. function:: double gsl_sf_bessel_J1 (double x)
              int gsl_sf_bessel_J1_e (double x, gsl_sf_result * result)

   These routines compute the regular cylindrical Bessel function of first
   order, :math:`J_1(x)`.

.. function:: double gsl_sf_bessel_Jn (int n, double x)
              int gsl_sf_bessel_Jn_e (int n, double x, gsl_sf_result * result)

   These routines compute the regular cylindrical Bessel function of 
   order :data:`n`, :math:`J_n(x)`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_Jn_array (int nmin, int nmax, double x, double result_array[])

   This routine computes the values of the regular cylindrical Bessel
   functions :math:`J_n(x)` for :math:`n` from :data:`nmin` to :data:`nmax`
   inclusive, storing the results in the array :data:`result_array`.  The
   values are computed using recurrence relations for efficiency, and
   therefore may differ slightly from the exact values.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW


Irregular Cylindrical Bessel Functions
--------------------------------------
.. index:: Irregular Cylindrical Bessel Functions
.. index::
   single: Y(x), Bessel Functions

.. function:: double gsl_sf_bessel_Y0 (double x)
              int gsl_sf_bessel_Y0_e (double x, gsl_sf_result * result)

   These routines compute the irregular cylindrical Bessel function of zeroth
   order, :math:`Y_0(x)`, for :math:`x>0`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_Y1 (double x)
              int gsl_sf_bessel_Y1_e (double x, gsl_sf_result * result)

   These routines compute the irregular cylindrical Bessel function of first
   order, :math:`Y_1(x)`, for :math:`x>0`.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_Yn (int n, double x)
              int gsl_sf_bessel_Yn_e (int n, double x, gsl_sf_result * result)

   These routines compute the irregular cylindrical Bessel function of 
   order :data:`n`, :math:`Y_n(x)`, for :math:`x>0`.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_Yn_array (int nmin, int nmax, double x, double result_array[])

   This routine computes the values of the irregular cylindrical Bessel
   functions :math:`Y_n(x)` for :math:`n` from :data:`nmin` to :data:`nmax`
   inclusive, storing the results in the array :data:`result_array`.  The
   domain of the function is :math:`x>0`.  The values are computed using
   recurrence relations for efficiency, and therefore may differ slightly
   from the exact values.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

Regular Modified Cylindrical Bessel Functions
---------------------------------------------
.. index:: Modified Cylindrical Bessel Functions
.. index:: Regular Modified Cylindrical Bessel Functions
.. index::
   single: I(x), Bessel Functions

.. function:: double gsl_sf_bessel_I0 (double x)
              int gsl_sf_bessel_I0_e (double x, gsl_sf_result * result)

   These routines compute the regular modified cylindrical Bessel function
   of zeroth order, :math:`I_0(x)`.
.. Exceptional Return Values: GSL_EOVRFLW

.. function:: double gsl_sf_bessel_I1 (double x)
              int gsl_sf_bessel_I1_e (double x, gsl_sf_result * result)

   These routines compute the regular modified cylindrical Bessel function
   of first order, :math:`I_1(x)`.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_In (int n, double x)
              int gsl_sf_bessel_In_e (int n, double x, gsl_sf_result * result)

   These routines compute the regular modified cylindrical Bessel function
   of order :data:`n`, :math:`I_n(x)`.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_In_array (int nmin, int nmax, double x, double result_array[])

   This routine computes the values of the regular modified cylindrical
   Bessel functions :math:`I_n(x)` for :math:`n` from :data:`nmin` to
   :data:`nmax` inclusive, storing the results in the array
   :data:`result_array`.  The start of the range :data:`nmin` must be positive
   or zero.  The values are computed using recurrence relations for
   efficiency, and therefore may differ slightly from the exact values.
.. Domain: nmin >=0, nmax >= nmin 
.. Conditions: n=nmin,...,nmax, nmin >=0, nmax >= nmin 
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_I0_scaled (double x)
              int gsl_sf_bessel_I0_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled regular modified cylindrical Bessel
   function of zeroth order :math:`\exp(-|x|) I_0(x)`.
.. Exceptional Return Values: none

.. function:: double gsl_sf_bessel_I1_scaled (double x)
              int gsl_sf_bessel_I1_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled regular modified cylindrical Bessel
   function of first order :math:`\exp(-|x|) I_1(x)`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_In_scaled (int n, double x)
              int gsl_sf_bessel_In_scaled_e (int n, double x, gsl_sf_result * result)

   These routines compute the scaled regular modified cylindrical Bessel
   function of order :data:`n`, :math:`\exp(-|x|) I_n(x)` 
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_In_scaled_array (int nmin, int nmax, double x, double result_array[])

   This routine computes the values of the scaled regular cylindrical
   Bessel functions :math:`\exp(-|x|) I_n(x)` for :math:`n` from
   :data:`nmin` to :data:`nmax` inclusive, storing the results in the array
   :data:`result_array`. The start of the range :data:`nmin` must be positive
   or zero.  The values are computed using recurrence relations for
   efficiency, and therefore may differ slightly from the exact values.
.. Domain: nmin >=0, nmax >= nmin 
.. Conditions:  n=nmin,...,nmax 
.. Exceptional Return Values: GSL_EUNDRFLW

Irregular Modified Cylindrical Bessel Functions
-----------------------------------------------
.. index:: Irregular Modified Cylindrical Bessel Functions
.. index::
   single: K(x), Bessel Functions

.. function:: double gsl_sf_bessel_K0 (double x)
              int gsl_sf_bessel_K0_e (double x, gsl_sf_result * result)

   These routines compute the irregular modified cylindrical Bessel
   function of zeroth order, :math:`K_0(x)`, for :math:`x > 0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_K1 (double x)
              int gsl_sf_bessel_K1_e (double x, gsl_sf_result * result)

   These routines compute the irregular modified cylindrical Bessel
   function of first order, :math:`K_1(x)`, for :math:`x > 0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_Kn (int n, double x)
              int gsl_sf_bessel_Kn_e (int n, double x, gsl_sf_result * result)

   These routines compute the irregular modified cylindrical Bessel
   function of order :data:`n`, :math:`K_n(x)`, for :math:`x > 0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_Kn_array (int nmin, int nmax, double x, double result_array[])

   This routine computes the values of the irregular modified cylindrical
   Bessel functions :math:`K_n(x)` for :math:`n` from :data:`nmin` to
   :data:`nmax` inclusive, storing the results in the array
   :data:`result_array`. The start of the range :data:`nmin` must be positive
   or zero. The domain of the function is :math:`x>0`. The values are
   computed using recurrence relations for efficiency, and therefore
   may differ slightly from the exact values.
.. Conditions: n=nmin,...,nmax 
.. Domain: x > 0.0, nmin>=0, nmax >= nmin
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_K0_scaled (double x)
              int gsl_sf_bessel_K0_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled irregular modified cylindrical Bessel
   function of zeroth order :math:`\exp(x) K_0(x)` for :math:`x>0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_bessel_K1_scaled (double x) 
              int gsl_sf_bessel_K1_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled irregular modified cylindrical Bessel
   function of first order :math:`\exp(x) K_1(x)` for :math:`x>0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_Kn_scaled (int n, double x)
              int gsl_sf_bessel_Kn_scaled_e (int n, double x, gsl_sf_result * result)

   These routines compute the scaled irregular modified cylindrical Bessel
   function of order :data:`n`, :math:`\exp(x) K_n(x)`, for :math:`x>0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_Kn_scaled_array (int nmin, int nmax, double x, double result_array[])

   This routine computes the values of the scaled irregular cylindrical
   Bessel functions :math:`\exp(x) K_n(x)` for :math:`n` from :data:`nmin` to
   :data:`nmax` inclusive, storing the results in the array
   :data:`result_array`. The start of the range :data:`nmin` must be positive
   or zero.  The domain of the function is :math:`x>0`. The values are
   computed using recurrence relations for efficiency, and therefore
   may differ slightly from the exact values.
.. Domain: x > 0.0, nmin >=0, nmax >= nmin 
.. Conditions: n=nmin,...,nmax 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

Regular Spherical Bessel Functions
----------------------------------
.. index:: Spherical Bessel Functions
.. index:: Regular Spherical Bessel Functions
.. index::
   single: j(x), Bessel Functions

.. function:: double gsl_sf_bessel_j0 (double x)
              int gsl_sf_bessel_j0_e (double x, gsl_sf_result * result)

   These routines compute the regular spherical Bessel function of zeroth
   order, :math:`j_0(x) = \sin(x)/x`.
.. Exceptional Return Values: none

.. function:: double gsl_sf_bessel_j1 (double x)
              int gsl_sf_bessel_j1_e (double x, gsl_sf_result * result)

   These routines compute the regular spherical Bessel function of first
   order, :math:`j_1(x) = (\sin(x)/x - \cos(x))/x`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_j2 (double x)
              int gsl_sf_bessel_j2_e (double x, gsl_sf_result * result)

   These routines compute the regular spherical Bessel function of second
   order, :math:`j_2(x) = ((3/x^2 - 1)\sin(x) - 3\cos(x)/x)/x`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_jl (int l, double x)
              int gsl_sf_bessel_jl_e (int l, double x, gsl_sf_result * result)

   These routines compute the regular spherical Bessel function of 
   order :data:`l`, :math:`j_l(x)`, for
   :math:`l \geq 0` and :math:`x \geq 0`.
.. Domain: l >= 0, x >= 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_jl_array (int lmax, double x, double result_array[])

   This routine computes the values of the regular spherical Bessel
   functions :math:`j_l(x)` for :math:`l` from 0 to :data:`lmax`
   inclusive  for
   :math:`lmax \geq 0` and
   :math:`x \geq 0`, storing the results in the array :data:`result_array`.
   The values are computed using recurrence relations for
   efficiency, and therefore may differ slightly from the exact values.
.. Domain: lmax >= 0 
.. Conditions: l=0,1,...,lmax 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_jl_steed_array (int lmax, double x, double * result_array)

   This routine uses Steed's method to compute the values of the regular
   spherical Bessel functions :math:`j_l(x)` for :math:`l` from 0 to
   :data:`lmax` inclusive for
   :math:`lmax \geq 0` and
   :math:`x \geq 0`, storing the results in the array
   :data:`result_array`.
   The Steed/Barnett algorithm is described in Comp. Phys. Comm. 21,
   297 (1981).  Steed's method is more stable than the
   recurrence used in the other functions but is also slower.
.. Domain: lmax >= 0 
.. Conditions: l=0,1,...,lmax 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

Irregular Spherical Bessel Functions
------------------------------------
.. index:: Irregular Spherical Bessel Functions
.. index::
   single: y(x), Bessel Functions

.. function:: double gsl_sf_bessel_y0 (double x)
              int gsl_sf_bessel_y0_e (double x, gsl_sf_result * result)

   These routines compute the irregular spherical Bessel function of zeroth
   order, :math:`y_0(x) = -\cos(x)/x`.
.. Exceptional Return Values: none

.. function:: double gsl_sf_bessel_y1 (double x)
              int gsl_sf_bessel_y1_e (double x, gsl_sf_result * result)

   These routines compute the irregular spherical Bessel function of first
   order, :math:`y_1(x) = -(\cos(x)/x + \sin(x))/x`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_y2 (double x)
              int gsl_sf_bessel_y2_e (double x, gsl_sf_result * result)

   These routines compute the irregular spherical Bessel function of second
   order, :math:`y_2(x) = (-3/x^3 + 1/x)\cos(x) - (3/x^2)\sin(x)`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_yl (int l, double x)
              int gsl_sf_bessel_yl_e (int l, double x, gsl_sf_result * result)

   These routines compute the irregular spherical Bessel function of 
   order :data:`l`, :math:`y_l(x)`, for
   :math:`l \geq 0`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_yl_array (int lmax, double x, double result_array[])

   This routine computes the values of the irregular spherical Bessel
   functions :math:`y_l(x)` for :math:`l` from 0 to :data:`lmax`
   inclusive for
   :math:`lmax \geq 0`, storing the results in the array :data:`result_array`.
   The values are computed using recurrence relations for
   efficiency, and therefore may differ slightly from the exact values.
.. Domain: lmax >= 0 
.. Conditions: l=0,1,...,lmax 
.. Exceptional Return Values: GSL_EUNDRFLW

Regular Modified Spherical Bessel Functions
-------------------------------------------
.. index:: Modified Spherical Bessel Functions
.. index:: Regular Modified Spherical Bessel Functions
.. index::
   single: i(x), Bessel Functions

The regular modified spherical Bessel functions :math:`i_l(x)` 
are related to the modified Bessel functions of fractional order,
:math:`i_l(x) = \sqrt{\pi/(2x)} I_{l+1/2}(x)`

.. function:: double gsl_sf_bessel_i0_scaled (double x)
              int gsl_sf_bessel_i0_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled regular modified spherical Bessel
   function of zeroth order, :math:`\exp(-|x|) i_0(x)`.
.. Exceptional Return Values: none

.. function:: double gsl_sf_bessel_i1_scaled (double x)
              int gsl_sf_bessel_i1_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled regular modified spherical Bessel
   function of first order, :math:`\exp(-|x|) i_1(x)`.
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_i2_scaled (double x)
              int gsl_sf_bessel_i2_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled regular modified spherical Bessel
   function of second order, :math:`\exp(-|x|) i_2(x)` 
.. Exceptional Return Values: GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_il_scaled (int l, double x)
              int gsl_sf_bessel_il_scaled_e (int l, double x, gsl_sf_result * result)

   These routines compute the scaled regular modified spherical Bessel
   function of order :data:`l`, :math:`\exp(-|x|) i_l(x)`
.. Domain: l >= 0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_il_scaled_array (int lmax, double x, double result_array[])

   This routine computes the values of the scaled regular modified
   spherical Bessel functions :math:`\exp(-|x|) i_l(x)` for :math:`l` from
   0 to :data:`lmax` inclusive for
   :math:`lmax \geq 0`, storing the results in
   the array :data:`result_array`. 
   The values are computed using recurrence relations for
   efficiency, and therefore may differ slightly from the exact values.
.. Domain: lmax >= 0 
.. Conditions: l=0,1,...,lmax 
.. Exceptional Return Values: GSL_EUNDRFLW

Irregular Modified Spherical Bessel Functions
---------------------------------------------
.. index:: Irregular Modified Spherical Bessel Functions
.. index::
   single: k(x), Bessel Functions

The irregular modified spherical Bessel functions :math:`k_l(x)`
are related to the irregular modified Bessel functions of fractional order,
:math:`k_l(x) = \sqrt{\pi/(2x)} K_{l+1/2}(x)`.

.. function:: double gsl_sf_bessel_k0_scaled (double x)
              int gsl_sf_bessel_k0_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled irregular modified spherical Bessel
   function of zeroth order, :math:`\exp(x) k_0(x)`, for :math:`x>0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_k1_scaled (double x)
              int gsl_sf_bessel_k1_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled irregular modified spherical Bessel
   function of first order, :math:`\exp(x) k_1(x)`, for :math:`x>0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: double gsl_sf_bessel_k2_scaled (double x)
              int gsl_sf_bessel_k2_scaled_e (double x, gsl_sf_result * result)

   These routines compute the scaled irregular modified spherical Bessel
   function of second order, :math:`\exp(x) k_2(x)`, for :math:`x>0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: double gsl_sf_bessel_kl_scaled (int l, double x)
              int gsl_sf_bessel_kl_scaled_e (int l, double x, gsl_sf_result * result)

   These routines compute the scaled irregular modified spherical Bessel
   function of order :data:`l`, :math:`\exp(x) k_l(x)`, for :math:`x>0`.
.. Domain: x > 0.0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_kl_scaled_array (int lmax, double x, double result_array[])

   This routine computes the values of the scaled irregular modified
   spherical Bessel functions :math:`\exp(x) k_l(x)` for :math:`l` from
   0 to :data:`lmax` inclusive for
   :math:`lmax \geq 0` and :math:`x>0`, storing the results in
   the array :data:`result_array`. 
   The values are computed using recurrence relations for
   efficiency, and therefore may differ slightly from the exact values.
.. Domain: lmax >= 0 
.. Conditions: l=0,1,...,lmax 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

Regular Bessel Function---Fractional Order
------------------------------------------
.. index::
   single: Fractional Order Bessel Functions
   single: Bessel Functions, Fractional Order
   single: Regular Bessel Functions, Fractional Order

.. function:: double gsl_sf_bessel_Jnu (double nu, double x)
              int gsl_sf_bessel_Jnu_e (double nu, double x, gsl_sf_result * result)

   These routines compute the regular cylindrical Bessel function of
   fractional order :math:`\nu`, :math:`J_\nu(x)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: int gsl_sf_bessel_sequence_Jnu_e (double nu, gsl_mode_t mode, size_t size, double v[])

   This function computes the regular cylindrical Bessel function of
   fractional order :math:`\nu`, :math:`J_\nu(x)`, evaluated at a series of
   :math:`x` values.  The array :data:`v` of length :data:`size` contains the
   :math:`x` values.  They are assumed to be strictly ordered and positive.
   The array is over-written with the values of :math:`J_\nu(x_i)`.
.. Exceptional Return Values: GSL_EDOM, GSL_EINVAL

Irregular Bessel Functions---Fractional Order
---------------------------------------------

.. function:: double gsl_sf_bessel_Ynu (double nu, double x)
              int gsl_sf_bessel_Ynu_e (double nu, double x, gsl_sf_result * result)

   These routines compute the irregular cylindrical Bessel function of
   fractional order :math:`\nu`, :math:`Y_\nu(x)`.
.. Exceptional Return Values: 

Regular Modified Bessel Functions---Fractional Order
----------------------------------------------------
.. index::
   single: Modified Bessel Functions, Fractional Order
   single: Regular Modified Bessel Functions, Fractional Order

.. function:: double gsl_sf_bessel_Inu (double nu, double x)
              int gsl_sf_bessel_Inu_e (double nu, double x, gsl_sf_result * result)

   These routines compute the regular modified Bessel function of
   fractional order :math:`\nu`, :math:`I_\nu(x)` for :math:`x>0`,
   :math:`\nu>0`.
.. Domain: x >= 0, nu >= 0 
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

.. function:: double gsl_sf_bessel_Inu_scaled (double nu, double x)
              int gsl_sf_bessel_Inu_scaled_e (double nu, double x, gsl_sf_result * result)

   These routines compute the scaled regular modified Bessel function of
   fractional order :math:`\nu`, :math:`\exp(-|x|)I_\nu(x)` for :math:`x>0`,
   :math:`\nu>0`.
.. Domain: x >= 0, nu >= 0 
.. Exceptional Return Values: GSL_EDOM

Irregular Modified Bessel Functions---Fractional Order
------------------------------------------------------
.. index::
   single: Irregular Modified Bessel Functions, Fractional Order

.. function:: double gsl_sf_bessel_Knu (double nu, double x)
              int gsl_sf_bessel_Knu_e (double nu, double x, gsl_sf_result * result)

   These routines compute the irregular modified Bessel function of
   fractional order :math:`\nu`, :math:`K_\nu(x)` for :math:`x>0`,
   :math:`\nu>0`.
.. Domain: x > 0, nu >= 0 
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW

.. function:: double gsl_sf_bessel_lnKnu (double nu, double x)
              int gsl_sf_bessel_lnKnu_e (double nu, double x, gsl_sf_result * result)

   These routines compute the logarithm of the irregular modified Bessel
   function of fractional order :math:`\nu`, :math:`\ln(K_\nu(x))` for
   :math:`x>0`, :math:`\nu>0`. 
.. Domain: x > 0, nu >= 0 
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_bessel_Knu_scaled (double nu, double x)
              int gsl_sf_bessel_Knu_scaled_e (double nu, double x, gsl_sf_result * result)

   These routines compute the scaled irregular modified Bessel function of
   fractional order :math:`\nu`, :math:`\exp(+|x|) K_\nu(x)` for :math:`x>0`,
   :math:`\nu>0`.
.. Domain: x > 0, nu >= 0 
.. Exceptional Return Values: GSL_EDOM

Zeros of Regular Bessel Functions
---------------------------------
.. index::
   single: Zeros of Regular Bessel Functions
   single: Regular Bessel Functions, Zeros of 

.. function:: double gsl_sf_bessel_zero_J0 (unsigned int s)
              int gsl_sf_bessel_zero_J0_e (unsigned int s, gsl_sf_result * result)

   These routines compute the location of the :data:`s`-th positive zero of
   the Bessel function :math:`J_0(x)`.
.. Exceptional Return Values: 

.. function:: double gsl_sf_bessel_zero_J1 (unsigned int s)
              int gsl_sf_bessel_zero_J1_e (unsigned int s, gsl_sf_result * result)

   These routines compute the location of the :data:`s`-th positive zero of
   the Bessel function :math:`J_1(x)`.
.. Exceptional Return Values: 

.. function:: double gsl_sf_bessel_zero_Jnu (double nu, unsigned int s)
              int gsl_sf_bessel_zero_Jnu_e (double nu, unsigned int s, gsl_sf_result * result)

   These routines compute the location of the :data:`s`-th positive zero of
   the Bessel function :math:`J_\nu(x)`.  The current implementation does not
   support negative values of :data:`nu`. 
.. Exceptional Return Values: 
