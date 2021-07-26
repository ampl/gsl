.. index::
   single: Legendre polynomials
   single: Legendre functions
   single: spherical harmonics
   single: conical functions
   single: hyperbolic space

The Legendre Functions and Legendre Polynomials are described in
Abramowitz & Stegun, Chapter 8.  These functions are declared in 
the header file :file:`gsl_sf_legendre.h`.

Legendre Polynomials
--------------------

.. function:: double gsl_sf_legendre_P1 (double x)
              double gsl_sf_legendre_P2 (double x)
              double gsl_sf_legendre_P3 (double x)
              int gsl_sf_legendre_P1_e (double x, gsl_sf_result * result)
              int gsl_sf_legendre_P2_e (double x, gsl_sf_result * result)
              int gsl_sf_legendre_P3_e (double x, gsl_sf_result * result)

   These functions evaluate the Legendre polynomials
   :math:`P_l(x)` using explicit representations for :math:`l = 1, 2, 3`.
.. Exceptional Return Values: none

.. function:: double gsl_sf_legendre_Pl (int l, double x)
              int gsl_sf_legendre_Pl_e (int l, double x, gsl_sf_result * result)

   These functions evaluate the Legendre polynomial :math:`P_l(x)`
   for a specific value of :data:`l`, :data:`x` subject to :math:`l \ge 0` and
   :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_legendre_Pl_array (int lmax, double x, double result_array[])
              int gsl_sf_legendre_Pl_deriv_array (int lmax, double x, double result_array[], double result_deriv_array[])

   These functions compute arrays of Legendre polynomials
   :math:`P_l(x)` and derivatives :math:`dP_l(x)/dx`
   for :math:`l = 0, \dots, lmax` and :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Q0 (double x)
              int gsl_sf_legendre_Q0_e (double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_0(x)` for
   :math:`x > -1` and :math:`x \ne 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Q1 (double x)
              int gsl_sf_legendre_Q1_e (double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_1(x)` for
   :math:`x > -1` and :math:`x \ne 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Ql (int l, double x)
              int gsl_sf_legendre_Ql_e (int l, double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_l(x)` for
   :math:`x > -1`, :math:`x \ne 1` and :math:`l \ge 0`.
.. Exceptional Return Values: GSL_EDOM

Associated Legendre Polynomials and Spherical Harmonics
-------------------------------------------------------

The following functions compute the associated Legendre polynomials
:math:`P_l^m(x)` which are solutions of the differential equation

.. only:: not texinfo

   .. math:: (1 - x^2) {d^2 \over dx^2} P_l^m(x) - 2x {d \over dx} P_l^m(x) +
             \left( l(l+1) - {m^2 \over 1 - x^2} \right) P_l^m(x) = 0

.. only:: texinfo

   ::

      (1 - x^2) d^2 P_l^m(x) / dx^2 P_l^m(x) - 2x d/dx P_l^m(x) +
      ( l(l+1) - m^2 / (1 - x^2) ) P_l^m(x) = 0

where the degree :math:`l` and order :math:`m` satisfy :math:`0 \le l` and
:math:`0 \le m \le l`.
The functions :math:`P_l^m(x)` grow combinatorially with
:math:`l` and can overflow for :math:`l` larger than about 150.
Alternatively, one may calculate normalized associated Legendre
polynomials. There are a number of different normalization conventions,
and these
functions can be stably computed up to degree and order 2700. The
following normalizations are provided:

* Schmidt semi-normalization

  Schmidt semi-normalized associated Legendre polynomials are often
  used in the magnetics community and are defined as

  .. only:: not texinfo

     .. math::

        S_l^0(x) &= P_l^0(x) \\
        S_l^m(x) &= (-1)^m \sqrt{2 {(l-m)! \over (l+m)!}} P_l^m(x), m > 0 

  .. only:: texinfo

     ::

        S_l^0(x) = P_l^0(x)
        S_l^m(x) = (-1)^m \sqrt((2(l-m)! / (l+m)!)) P_l^m(x), m > 0 

  The factor of :math:`(-1)^m` is called the Condon-Shortley phase
  factor and can be excluded if desired by setting the parameter
  :code:`csphase = 1` in the functions below.

* Spherical Harmonic Normalization

  The associated Legendre polynomials suitable for calculating spherical
  harmonics are defined as

  .. only:: not texinfo

     .. math:: Y_l^m(x) = (-1)^m \sqrt{{2l + 1 \over 4 \pi} {(l-m)! \over (l+m)!}} P_l^m(x)

  .. only:: texinfo

     ::

        Y_l^m(x) = (-1)^m \sqrt((2l + 1) * (l-m)! / (4 \pi) / (l+m)!) P_l^m(x)

  where again the phase factor :math:`(-1)^m` can be included or excluded
  if desired.

* Full Normalization

  The fully normalized associated Legendre polynomials are defined as

  .. only:: not texinfo

     .. math:: N_l^m(x) = (-1)^m \sqrt{(l + {1 \over 2}) {(l-m)! \over (l+m)!}} P_l^m(x)

  .. only:: texinfo
  
     ::
     
        N_l^m(x) = (-1)^m \sqrt((l + 1/2) (l-m)! / (l+m)!) P_l^m(x)

  and have the property

  .. math:: \int_{-1}^1 N_l^m(x)^2 dx = 1

The normalized associated Legendre routines below use a recurrence
relation which is stable up to a degree and order of about 2700.
Beyond this, the computed functions could suffer from underflow
leading to incorrect results. Routines are provided to compute
first and second derivatives
:math:`dP_l^m(x)/dx` and :math:`d^2 P_l^m(x)/dx^2` as well as their alternate
versions :math:`d P_l^m(\cos{\theta})/d\theta` and
:math:`d^2 P_l^m(\cos{\theta})/d\theta^2`. While there is a simple
scaling relationship between the two forms, the derivatives
involving :math:`\theta` are heavily used in spherical harmonic
expansions and so these routines are also provided.

In the functions below, a parameter of type :type:`gsl_sf_legendre_t`
specifies the type of normalization to use. The possible values are

.. type:: gsl_sf_legendre_t

   ================================== ===============================================================================
   Value                              Description
   ================================== ===============================================================================
   :code:`GSL_SF_LEGENDRE_NONE`       The unnormalized associated Legendre polynomials :math:`P_l^m(x)`
   :code:`GSL_SF_LEGENDRE_SCHMIDT`    The Schmidt semi-normalized associated Legendre polynomials :math:`S_l^m(x)`
   :code:`GSL_SF_LEGENDRE_SPHARM`     The spherical harmonic associated Legendre polynomials :math:`Y_l^m(x)`
   :code:`GSL_SF_LEGENDRE_FULL`       The fully normalized associated Legendre polynomials :math:`N_l^m(x)`
   ================================== ===============================================================================

.. function:: int gsl_sf_legendre_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[])
              int gsl_sf_legendre_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[])

   These functions calculate all normalized associated Legendre
   polynomials for :math:`0 \le l \le lmax` and
   :math:`0 \le m \le l` for :math:`|x| \le 1`.
   The :data:`norm` parameter specifies which normalization is used.
   The normalized :math:`P_l^m(x)` values are stored in :data:`result_array`, whose
   minimum size can be obtained from calling :func:`gsl_sf_legendre_array_n`.
   The array index of :math:`P_l^m(x)` is obtained from calling
   :code:`gsl_sf_legendre_array_index(l, m)`. To include or exclude
   the Condon-Shortley phase factor of :math:`(-1)^m`, set the parameter
   :data:`csphase` to either :math:`-1` or :math:`1` respectively in the
   :code:`_e` function. This factor is excluded by default.

.. function:: int gsl_sf_legendre_deriv_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[])
              int gsl_sf_legendre_deriv_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[], double result_deriv_array[])

   These functions calculate all normalized associated Legendre
   functions and their first derivatives up to degree :data:`lmax` for
   :math:`|x| < 1`.
   The parameter :data:`norm` specifies the normalization used. The
   normalized :math:`P_l^m(x)` values and their derivatives
   :math:`dP_l^m(x)/dx` are stored in :data:`result_array` and
   :data:`result_deriv_array` respectively.
   To include or exclude
   the Condon-Shortley phase factor of :math:`(-1)^m`, set the parameter
   :data:`csphase` to either :math:`-1` or :math:`1` respectively in the
   :code:`_e` function. This factor is excluded by default.

.. function:: int gsl_sf_legendre_deriv_alt_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[])
              int gsl_sf_legendre_deriv_alt_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[], double result_deriv_array[])

   These functions calculate all normalized associated Legendre
   functions and their (alternate) first derivatives up to degree :data:`lmax` for
   :math:`|x| < 1`.
   The normalized :math:`P_l^m(x)` values and their derivatives
   :math:`dP_l^m(\cos{\theta})/d\theta` are stored in :data:`result_array` and
   :data:`result_deriv_array` respectively.
   To include or exclude
   the Condon-Shortley phase factor of :math:`(-1)^m`, set the parameter
   :data:`csphase` to either :math:`-1` or :math:`1` respectively in the
   :code:`_e` function. This factor is excluded by default.

.. function:: int gsl_sf_legendre_deriv2_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[], double result_deriv2_array[])
              int gsl_sf_legendre_deriv2_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[], double result_deriv_array[], double result_deriv2_array[])

   These functions calculate all normalized associated Legendre
   functions and their first and second derivatives up to degree :data:`lmax` for
   :math:`|x| < 1`.
   The parameter :data:`norm` specifies the normalization used. The
   normalized :math:`P_l^m(x)`, their first derivatives
   :math:`dP_l^m(x)/dx`, and their second derivatives
   :math:`d^2 P_l^m(x)/dx^2` are stored in :data:`result_array`,
   :data:`result_deriv_array`, and :data:`result_deriv2_array` respectively.
   To include or exclude
   the Condon-Shortley phase factor of :math:`(-1)^m`, set the parameter
   :data:`csphase` to either :math:`-1` or :math:`1` respectively in the
   :code:`_e` function. This factor is excluded by default.

.. function:: int gsl_sf_legendre_deriv2_alt_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[], double result_deriv2_array[])
              int gsl_sf_legendre_deriv2_alt_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[], double result_deriv_array[], double result_deriv2_array[])

   These functions calculate all normalized associated Legendre
   functions and their (alternate) first and second derivatives up to degree
   :data:`lmax` for
   :math:`|x| < 1`.
   The parameter :data:`norm` specifies the normalization used. The
   normalized :math:`P_l^m(x)`, their first derivatives
   :math:`dP_l^m(\cos{\theta})/d\theta`, and their second derivatives
   :math:`d^2 P_l^m(\cos{\theta})/d\theta^2` are stored in :data:`result_array`,
   :data:`result_deriv_array`, and :data:`result_deriv2_array` respectively.
   To include or exclude
   the Condon-Shortley phase factor of :math:`(-1)^m`, set the parameter
   :data:`csphase` to either :math:`-1` or :math:`1` respectively in the
   :code:`_e` function. This factor is excluded by default.

.. function:: size_t gsl_sf_legendre_nlm(const size_t lmax)

   This function returns the total number of associated Legendre
   functions :math:`P_l^m(x)` for a given :data:`lmax`. The number is
   :code:`(lmax+1) * (lmax+2) / 2`.

.. function:: size_t gsl_sf_legendre_array_n (const size_t lmax)

   This function returns the minimum array size for maximum degree :data:`lmax`
   needed for the array versions of the associated Legendre functions.
   Size is calculated as the total number of :math:`P_l^m(x)` functions
   (see :func:`gsl_sf_legendre_nlm`),
   plus extra space for precomputing multiplicative factors used in the
   recurrence relations.

.. function:: size_t gsl_sf_legendre_array_index (const size_t l, const size_t m)

   This function returns the index into :data:`result_array`,
   :data:`result_deriv_array`, or :data:`result_deriv2_array` corresponding
   to :math:`P_l^m(x)`, :math:`P_l^{'m}(x)`, or :math:`P_l^{''m}(x)`. The
   index is given by :math:`l(l+1)/2 + m`.

   An inline version of this function is used if :macro:`HAVE_INLINE` is
   defined.

.. function:: double gsl_sf_legendre_Plm (int l, int m, double x)
              int gsl_sf_legendre_Plm_e (int l, int m, double x, gsl_sf_result * result)

   These routines compute the associated Legendre polynomial
   :math:`P_l^m(x)` for :math:`m \ge 0`,
   :math:`l \ge m`, and :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

.. function:: double gsl_sf_legendre_sphPlm (int l, int m, double x)
              int gsl_sf_legendre_sphPlm_e (int l, int m, double x, gsl_sf_result * result)

   These routines compute the normalized associated Legendre polynomial
   :math:`\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x)` suitable
   for use in spherical harmonics.  The parameters must satisfy :math:`m \ge 0`,
   :math:`l \ge m`, and :math:`|x| \le 1`.
   These routines avoid the overflows
   that occur for the standard normalization of :math:`P_l^m(x)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_legendre_Plm_array (int lmax, int m, double x, double result_array[])
              int gsl_sf_legendre_Plm_deriv_array (int lmax, int m, double x, double result_array[], double result_deriv_array[])

   These functions are now deprecated and will be removed in a future
   release; see :func:`gsl_sf_legendre_array` and
   :func:`gsl_sf_legendre_deriv_array`.

.. function:: int gsl_sf_legendre_sphPlm_array (int lmax, int m, double x, double result_array[])
              int gsl_sf_legendre_sphPlm_deriv_array (int lmax, int m, double x, double result_array[], double result_deriv_array[])

   These functions are now deprecated and will be removed in a future
   release; see :func:`gsl_sf_legendre_array` and
   :func:`gsl_sf_legendre_deriv_array`.

.. function:: int gsl_sf_legendre_array_size (const int lmax, const int m)

   This function is now deprecated and will be removed in a future
   release.

Conical Functions
-----------------

The Conical Functions :math:`P^\mu_{-(1/2)+i\lambda}(x)`
and :math:`Q^\mu_{-(1/2)+i\lambda}`
are described in Abramowitz & Stegun, Section 8.12.

.. function:: double gsl_sf_conicalP_half (double lambda, double x)
              int gsl_sf_conicalP_half_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the irregular Spherical Conical Function
   :math:`P^{1/2}_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_mhalf (double lambda, double x)
              int gsl_sf_conicalP_mhalf_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the regular Spherical Conical Function
   :math:`P^{-1/2}_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_0 (double lambda, double x)
              int gsl_sf_conicalP_0_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the conical function
   :math:`P^0_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_1 (double lambda, double x)
              int gsl_sf_conicalP_1_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the conical function 
   :math:`P^1_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_sph_reg (int l, double lambda, double x)
              int gsl_sf_conicalP_sph_reg_e (int l, double lambda, double x, gsl_sf_result * result)

   These routines compute the Regular Spherical Conical Function
   :math:`P^{-1/2-l}_{-1/2 + i \lambda}(x)`
   for :math:`x > -1` and :math:`l \ge -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_cyl_reg (int m, double lambda, double x)
              int gsl_sf_conicalP_cyl_reg_e (int m, double lambda, double x, gsl_sf_result * result)

   These routines compute the Regular Cylindrical Conical Function
   :math:`P^{-m}_{-1/2 + i \lambda}(x)`
   for :math:`x > -1` and :math:`m \ge -1`.
.. Exceptional Return Values: GSL_EDOM

Radial Functions for Hyperbolic Space
-------------------------------------

The following spherical functions are specializations of Legendre
functions which give the regular eigenfunctions of the Laplacian on a
3-dimensional hyperbolic space :math:`H^3`.  Of particular interest is
the flat limit, :math:`\lambda \to \infty`, :math:`\eta \to 0`,
:math:`\lambda\eta` fixed.
  
.. function:: double gsl_sf_legendre_H3d_0 (double lambda, double eta)
              int gsl_sf_legendre_H3d_0_e (double lambda, double eta, gsl_sf_result * result)

   These routines compute the zeroth radial eigenfunction of the Laplacian on the
   3-dimensional hyperbolic space,

   .. math:: L^{H3d}_0(\lambda,\eta) := {\sin(\lambda\eta) \over \lambda\sinh(\eta)}

   for :math:`\eta \ge 0`.
   In the flat limit this takes the form
   :math:`L^{H3d}_0(\lambda,\eta) = j_0(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_H3d_1 (double lambda, double eta)
              int gsl_sf_legendre_H3d_1_e (double lambda, double eta, gsl_sf_result * result)

   These routines compute the first radial eigenfunction of the Laplacian on
   the 3-dimensional hyperbolic space,

   .. math:: L^{H3d}_1(\lambda,\eta) := {1\over\sqrt{\lambda^2 + 1}} {\left(\sin(\lambda \eta)\over \lambda \sinh(\eta)\right)} \left(\coth(\eta) - \lambda \cot(\lambda\eta)\right)

   for :math:`\eta \ge 0`
   In the flat limit this takes the form 
   :math:`L^{H3d}_1(\lambda,\eta) = j_1(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_H3d (int l, double lambda, double eta)
              int gsl_sf_legendre_H3d_e (int l, double lambda, double eta, gsl_sf_result * result)

   These routines compute the :data:`l`-th radial eigenfunction of the
   Laplacian on the 3-dimensional hyperbolic space :math:`\eta \ge 0` and
   :math:`l \ge 0`.
   In the flat limit this takes the form
   :math:`L^{H3d}_l(\lambda,\eta) = j_l(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_legendre_H3d_array (int lmax, double lambda, double eta, double result_array[])

   This function computes an array of radial eigenfunctions
   :math:`L^{H3d}_l( \lambda, \eta)`
   for :math:`0 \le l \le lmax`.
.. Exceptional Return Values:
