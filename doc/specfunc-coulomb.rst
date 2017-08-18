.. index::
   single: Coulomb wave functions
   single: hydrogen atom

The prototypes of the Coulomb functions are declared in the header file
:file:`gsl_sf_coulomb.h`.  Both bound state and scattering solutions are
available.

Normalized Hydrogenic Bound States
----------------------------------

.. function:: double gsl_sf_hydrogenicR_1 (double Z, double r)
              int gsl_sf_hydrogenicR_1_e (double Z, double r, gsl_sf_result * result)

   These routines compute the lowest-order normalized hydrogenic bound
   state radial wavefunction
   :math:`R_1 := 2Z \sqrt{Z} \exp(-Z r)`.

.. function:: double gsl_sf_hydrogenicR (int n, int l, double Z, double r)
              int gsl_sf_hydrogenicR_e (int n, int l, double Z, double r, gsl_sf_result * result)

   These routines compute the :data:`n`-th normalized hydrogenic bound state
   radial wavefunction,

   .. only:: not texinfo

      .. math:: R_n := {2 Z^{3/2} \over n^2}  \left({2Z r \over n}\right)^l  \sqrt{(n-l-1)! \over (n+l)!} \exp(-Z r/n) L^{2l+1}_{n-l-1}(2Z r / n).

   .. only:: texinfo

      | R_n := 2 (Z^{3/2}/n^2) \sqrt{(n-l-1)!/(n+l)!} \exp(-Z r/n) (2Zr/n)^l
      |           L^{2l+1}_{n-l-1}(2Zr/n).  

   where :math:`L^a_b(x)` is the :ref:`generalized Laguerre polynomial <laguerre-functions>`.
   The normalization is chosen such that the wavefunction :math:`\psi` is
   given by :math:`\psi(n,l,r) = R_n Y_{lm}`.   

Coulomb Wave Functions
----------------------

The Coulomb wave functions :math:`F_L(\eta,x)`, :math:`G_L(\eta,x)` are
described in Abramowitz & Stegun, Chapter 14.  Because there can be a
large dynamic range of values for these functions, overflows are handled
gracefully.  If an overflow occurs, :code:`GSL_EOVRFLW` is signalled and
exponent(s) are returned through the modifiable parameters :data:`exp_F`,
:data:`exp_G`. The full solution can be reconstructed from the following
relations,

.. only:: not texinfo

   .. math::

      F_L(\eta,x) &= fc[k_L] * \exp(exp_F) \\
      G_L(\eta,x) &= gc[k_L] * \exp(exp_G)

   .. math::

      F_L'(\eta,x) &= fcp[k_L] * \exp(exp_F) \\
      G_L'(\eta,x) &= gcp[k_L] * \exp(exp_G)

.. only:: texinfo

   | F_L(\eta,x) = fc[k_L] * \exp(exp_F)
   | G_L(\eta,x) = gc[k_L] * \exp(exp_G)
   |
   | F_L'(\eta,x) = fcp[k_L] * \exp(exp_F)
   | G_L'(\eta,x) = gcp[k_L] * \exp(exp_G)

.. function:: int gsl_sf_coulomb_wave_FG_e (double eta, double x, double L_F, int k, gsl_sf_result * F, gsl_sf_result * Fp, gsl_sf_result * G, gsl_sf_result * Gp, double * exp_F, double * exp_G)

   This function computes the Coulomb wave functions :math:`F_L(\eta,x)`,
   :math:`G_{L-k}(\eta,x)` and their derivatives 
   :math:`F'_L(\eta,x)`, 
   :math:`G'_{L-k}(\eta,x)`
   with respect to :math:`x`.  The parameters are restricted to :math:`L, L-k > -1/2`,
   :math:`x > 0` and integer :math:`k`.  Note that :math:`L`
   itself is not restricted to being an integer. The results are stored in
   the parameters F, G for the function values and :data:`Fp`,
   :data:`Gp` for the derivative values.  If an overflow occurs,
   :code:`GSL_EOVRFLW` is returned and scaling exponents are stored in
   the modifiable parameters :data:`exp_F`, :data:`exp_G`.

.. function:: int gsl_sf_coulomb_wave_F_array (double L_min, int kmax, double eta, double x, double fc_array[], double * F_exponent)

   This function computes the Coulomb wave function :math:`F_L(\eta,x)` for
   :math:`L = Lmin \dots Lmin + kmax`, storing the results in :data:`fc_array`.
   In the case of overflow the exponent is stored in :data:`F_exponent`.

.. function:: int gsl_sf_coulomb_wave_FG_array (double L_min, int kmax, double eta, double x, double fc_array[], double gc_array[], double * F_exponent, double * G_exponent)

   This function computes the functions :math:`F_L(\eta,x)`,
   :math:`G_L(\eta,x)` for :math:`L = Lmin \dots Lmin + kmax` storing the
   results in :data:`fc_array` and :data:`gc_array`.  In the case of overflow the
   exponents are stored in :data:`F_exponent` and :data:`G_exponent`.

.. function:: int gsl_sf_coulomb_wave_FGp_array (double L_min, int kmax, double eta, double x, double fc_array[], double fcp_array[], double gc_array[], double gcp_array[], double * F_exponent, double * G_exponent)

   This function computes the functions :math:`F_L(\eta,x)`,
   :math:`G_L(\eta,x)` and their derivatives :math:`F'_L(\eta,x)`,
   :math:`G'_L(\eta,x)` for :math:`L = Lmin \dots Lmin + kmax` storing the
   results in :data:`fc_array`, :data:`gc_array`, :data:`fcp_array` and :data:`gcp_array`.
   In the case of overflow the exponents are stored in :data:`F_exponent` 
   and :data:`G_exponent`.

.. function:: int gsl_sf_coulomb_wave_sphF_array (double L_min, int kmax, double eta, double x, double fc_array[], double F_exponent[])

   This function computes the Coulomb wave function divided by the argument
   :math:`F_L(\eta, x)/x` for :math:`L = Lmin \dots Lmin + kmax`, storing the
   results in :data:`fc_array`.  In the case of overflow the exponent is
   stored in :data:`F_exponent`. This function reduces to spherical Bessel
   functions in the limit :math:`\eta \to 0`.

Coulomb Wave Function Normalization Constant
--------------------------------------------

The Coulomb wave function normalization constant is defined in
Abramowitz 14.1.7.

.. function:: int gsl_sf_coulomb_CL_e (double L, double eta, gsl_sf_result * result)

   This function computes the Coulomb wave function normalization constant
   :math:`C_L(\eta)` for :math:`L > -1`.

.. function:: int gsl_sf_coulomb_CL_array (double Lmin, int kmax, double eta, double cl[])

   This function computes the Coulomb wave function normalization constant
   :math:`C_L(\eta)` for :math:`L = Lmin \dots Lmin + kmax`, :math:`Lmin > -1`.
