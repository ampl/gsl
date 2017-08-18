.. index:: elliptic integrals

The functions described in this section are declared in the header
file :file:`gsl_sf_ellint.h`.  Further information about the elliptic
integrals can be found in Abramowitz & Stegun, Chapter 17.

Definition of Legendre Forms
----------------------------
.. index:: Legendre forms of elliptic integrals

The Legendre forms of elliptic integrals :math:`F(\phi,k)`,
:math:`E(\phi,k)` and :math:`\Pi(\phi,k,n)` are defined by,

.. only:: not texinfo

   .. math::

      F(\phi,k)   &= \int_0^\phi dt {1 \over \sqrt{(1 - k^2 \sin^2(t))}} \\
      E(\phi,k)   &= \int_0^\phi dt   \sqrt{(1 - k^2 \sin^2(t))} \\
      \Pi(\phi,k,n) &= \int_0^\phi dt {1 \over (1 + n \sin^2(t)) \sqrt{1 - k^2 \sin^2(t)}}

.. only:: texinfo

   ::

      F(\phi,k) = \int_0^\phi dt 1/\sqrt((1 - k^2 \sin^2(t)))
      E(\phi,k) = \int_0^\phi dt   \sqrt((1 - k^2 \sin^2(t)))
      Pi(\phi,k,n) = \int_0^\phi dt 1/((1 + n \sin^2(t))\sqrt(1 - k^2 \sin^2(t)))

The complete Legendre forms are denoted by :math:`K(k) = F(\pi/2, k)` and
:math:`E(k) = E(\pi/2, k)`.  

The notation used here is based on Carlson, "Numerische
Mathematik" 33 (1979) 1 and differs slightly from that used by
Abramowitz & Stegun, where the functions are given in terms of the
parameter :math:`m = k^2` and :math:`n` is replaced by :math:`-n`.

Definition of Carlson Forms
---------------------------
.. index:: Carlson forms of Elliptic integrals

The Carlson symmetric forms of elliptical integrals :math:`RC(x,y)`,
:math:`RD(x,y,z)`, :math:`RF(x,y,z)` and :math:`RJ(x,y,z,p)` are defined
by,

.. only:: not texinfo

   .. math::

      RC(x,y)   &= 1/2 \int_0^\infty dt (t+x)^{-1/2} (t+y)^{-1} \\
      RD(x,y,z) &= 3/2 \int_0^\infty dt (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-3/2} \\
      RF(x,y,z) &= 1/2 \int_0^\infty dt (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-1/2} \\
      RJ(x,y,z,p) &= 3/2 \int_0^\infty dt (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-1/2} (t+p)^{-1}

.. only:: texinfo

   ::

     RC(x,y) = 1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1)
     RD(x,y,z) = 3/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2)
     RF(x,y,z) = 1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2)
     RJ(x,y,z,p) = 3/2 \int_0^\infty dt 
                      (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1)

Legendre Form of Complete Elliptic Integrals
--------------------------------------------

.. function:: double gsl_sf_ellint_Kcomp (double k, gsl_mode_t mode)
              int gsl_sf_ellint_Kcomp_e (double k, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the complete elliptic integral :math:`K(k)` to
   the accuracy specified by the mode variable :data:`mode`.  
   Note that Abramowitz & Stegun define this function in terms of the
   parameter :math:`m = k^2`.
.. Exceptional Return Values:  GSL_EDOM

.. function:: double gsl_sf_ellint_Ecomp (double k, gsl_mode_t mode)
              int gsl_sf_ellint_Ecomp_e (double k, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the complete elliptic integral :math:`E(k)` to the
   accuracy specified by the mode variable :data:`mode`.
   Note that Abramowitz & Stegun define this function in terms of the
   parameter :math:`m = k^2`.
.. Exceptional Return Values:  GSL_EDOM

.. function:: double gsl_sf_ellint_Pcomp (double k, double n, gsl_mode_t mode)
              int gsl_sf_ellint_Pcomp_e (double k, double n,  gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the complete elliptic integral :math:`\Pi(k,n)` to the
   accuracy specified by the mode variable :data:`mode`.
   Note that Abramowitz & Stegun define this function in terms of the
   parameters :math:`m = k^2` and :math:`\sin^2(\alpha) = k^2`, with the
   change of sign :math:`n \to -n`.
.. Exceptional Return Values:  GSL_EDOM

Legendre Form of Incomplete Elliptic Integrals
----------------------------------------------

.. function:: double gsl_sf_ellint_F (double phi, double k, gsl_mode_t mode)
              int gsl_sf_ellint_F_e (double phi, double k, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the incomplete elliptic integral :math:`F(\phi,k)`
   to the accuracy specified by the mode variable :data:`mode`.
   Note that Abramowitz & Stegun define this function in terms of the
   parameter :math:`m = k^2`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_ellint_E (double phi, double k, gsl_mode_t mode)
              int gsl_sf_ellint_E_e (double phi, double k, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the incomplete elliptic integral :math:`E(\phi,k)`
   to the accuracy specified by the mode variable :data:`mode`.
   Note that Abramowitz & Stegun define this function in terms of the
   parameter :math:`m = k^2`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_ellint_P (double phi, double k, double n, gsl_mode_t mode)
              int gsl_sf_ellint_P_e (double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the incomplete elliptic integral :math:`\Pi(\phi,k,n)`
   to the accuracy specified by the mode variable :data:`mode`.
   Note that Abramowitz & Stegun define this function in terms of the
   parameters :math:`m = k^2` and :math:`\sin^2(\alpha) = k^2`, with the
   change of sign :math:`n \to -n`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_ellint_D (double phi, double k, gsl_mode_t mode)
              int gsl_sf_ellint_D_e (double phi, double k, gsl_mode_t mode, gsl_sf_result * result)

   These functions compute the incomplete elliptic integral
   :math:`D(\phi,k)` which is defined through the Carlson form :math:`RD(x,y,z)`
   by the following relation, 

   .. only:: not texinfo

      .. math:: D(\phi,k) = {1 \over 3} (\sin \phi)^3 RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1)


   .. only:: texinfo

      ::

        D(\phi,k) = (1/3)(\sin(\phi))^3 RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1).

.. Exceptional Return Values: GSL_EDOM

Carlson Forms
-------------

.. function:: double gsl_sf_ellint_RC (double x, double y, gsl_mode_t mode)
              int gsl_sf_ellint_RC_e (double x, double y, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the incomplete elliptic integral :math:`RC(x,y)`
   to the accuracy specified by the mode variable :data:`mode`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_ellint_RD (double x, double y, double z, gsl_mode_t mode)
              int gsl_sf_ellint_RD_e (double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the incomplete elliptic integral :math:`RD(x,y,z)`
   to the accuracy specified by the mode variable :data:`mode`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_ellint_RF (double x, double y, double z, gsl_mode_t mode)
              int gsl_sf_ellint_RF_e (double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the incomplete elliptic integral :math:`RF(x,y,z)`
   to the accuracy specified by the mode variable :data:`mode`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_ellint_RJ (double x, double y, double z, double p, gsl_mode_t mode)
              int gsl_sf_ellint_RJ_e (double x, double y, double z, double p, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the incomplete elliptic integral :math:`RJ(x,y,z,p)`
   to the accuracy specified by the mode variable :data:`mode`.
.. Exceptional Return Values: GSL_EDOM
