.. index::
   single: coupling coefficients
   single: 3-j symbols
   single: 6-j symbols
   single: 9-j symbols
   single: Wigner coefficients
   single: Racah coefficients

The Wigner 3-j, 6-j and 9-j symbols give the coupling coefficients for
combined angular momentum vectors.  Since the arguments of the standard
coupling coefficient functions are integer or half-integer, the
arguments of the following functions are, by convention, integers equal
to twice the actual spin value.  For information on the 3-j coefficients
see Abramowitz & Stegun, Section 27.9.  The functions described in this
section are declared in the header file :file:`gsl_sf_coupling.h`.

3-j Symbols
-----------

.. function:: double gsl_sf_coupling_3j (int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc)
              int gsl_sf_coupling_3j_e (int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc, gsl_sf_result * result)

   These routines compute the Wigner 3-j coefficient, 

   .. only:: not texinfo

      .. math::

         \left(
         \begin{array}{ccc}
           ja & jb & jc \\
           ma & mb & mc
         \end{array}
         \right)

   .. only:: texinfo

      | ( ja jb jc )
      | ( ma mb mc )

   where the arguments are given in half-integer units, :math:`ja` =
   :data:`two_ja`/2, :math:`ma` = :data:`two_ma`/2, etc.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

6-j Symbols
-----------

.. function:: double gsl_sf_coupling_6j (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf)
              int gsl_sf_coupling_6j_e (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, gsl_sf_result * result) 

   These routines compute the Wigner 6-j coefficient, 

   .. only:: not texinfo

      .. math::

         \left\{
         \begin{array}{ccc}
           ja & jb & jc \\
           jd & je & jf
         \end{array}
         \right\}

   .. only:: texinfo

      | { ja jb jc }
      | { jd je jf }

   where the arguments are given in half-integer units, :math:`ja` =
   :data:`two_ja`/2, :math:`ma` = :data:`two_ma`/2, etc.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

9-j Symbols
-----------

.. function:: double gsl_sf_coupling_9j (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji)
              int gsl_sf_coupling_9j_e (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji, gsl_sf_result * result) 

   These routines compute the Wigner 9-j coefficient, 

   .. only:: not texinfo

      .. math::

         \left\{
         \begin{array}{ccc}
           ja & jb & jc \\
           jd & je & jf \\
           jg & jh & ji
         \end{array}
         \right\}

   .. only:: texinfo
   
      | { ja jb jc }
      | { jd je jf }
      | { jg jh ji }

   where the arguments are given in half-integer units, :math:`ja` =
   :data:`two_ja`/2, :math:`ma` = :data:`two_ma`/2, etc.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW
