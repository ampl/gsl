.. index::
   single: hypergeometric functions
   single: confluent hypergeometric functions

Hypergeometric functions are described in Abramowitz & Stegun, Chapters
13 and 15.  These functions are declared in the header file
:file:`gsl_sf_hyperg.h`.

.. function:: double gsl_sf_hyperg_0F1 (double c, double x)
              int gsl_sf_hyperg_0F1_e (double c, double x, gsl_sf_result * result)

   These routines compute the hypergeometric function
   
   .. only:: not texinfo
      
      .. math:: {}_0F_1(c,x)

   .. only:: texinfo

      .. math:: 0F1(c,x)

.. It is related to Bessel functions
.. 0F1[c,x] =
..   Gamma[c]    x^(1/2(1-c)) I_(c-1)(2 Sqrt[x])
..   Gamma[c] (-x)^(1/2(1-c)) J_(c-1)(2 Sqrt[-x])
.. exceptions: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_hyperg_1F1_int (int m, int n, double x)
              int gsl_sf_hyperg_1F1_int_e (int m, int n, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function

   .. only:: not texinfo

      .. math:: {}_1F_1(m,n,x) = M(m,n,x)

   .. only:: texinfo

      .. math:: 1F1(m,n,x) = M(m,n,x)

   for integer parameters :data:`m`, :data:`n`.
.. exceptions: 

.. function:: double gsl_sf_hyperg_1F1 (double a, double b, double x)
              int gsl_sf_hyperg_1F1_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function

   .. only:: not texinfo

      .. math:: {}_1F_1(a,b,x) = M(a,b,x)

   .. only:: texinfo

      .. math:: 1F1(a,b,x) = M(a,b,x)

   for general parameters :data:`a`, :data:`b`.
.. exceptions:

.. function:: double gsl_sf_hyperg_U_int (int m, int n, double x)
              int gsl_sf_hyperg_U_int_e (int m, int n, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function
   :math:`U(m,n,x)` for integer parameters :data:`m`, :data:`n`.
.. exceptions:

.. function:: int gsl_sf_hyperg_U_int_e10_e (int m, int n, double x, gsl_sf_result_e10 * result)

   This routine computes the confluent hypergeometric function
   :math:`U(m,n,x)` for integer parameters :data:`m`, :data:`n` using the
   :type:`gsl_sf_result_e10` type to return a result with extended range.

.. function:: double gsl_sf_hyperg_U (double a, double b, double x)
              int gsl_sf_hyperg_U_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function :math:`U(a,b,x)`.
.. exceptions:

.. function:: int gsl_sf_hyperg_U_e10_e (double a, double b, double x, gsl_sf_result_e10 * result)

   This routine computes the confluent hypergeometric function
   :math:`U(a,b,x)` using the :type:`gsl_sf_result_e10` type to return a
   result with extended range. 
.. exceptions:

.. function:: double gsl_sf_hyperg_2F1 (double a, double b, double c, double x)
              int gsl_sf_hyperg_2F1_e (double a, double b, double c, double x, gsl_sf_result * result)

   These routines compute the Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a,b,c,x) = F(a,b,c,x)

   .. only:: texinfo

      .. math:: 2F1(a,b,c,x) = F(a,b,c,x)
         
   for :math:`|x| < 1`. If the arguments :math:`(a,b,c,x)` are too close to a singularity then
   the function can return the error code :macro:`GSL_EMAXITER` when the
   series approximation converges too slowly.  This occurs in the region of
   :math:`x = 1`, :math:`c - a - b = m` for integer m.
.. exceptions:

.. function:: double gsl_sf_hyperg_2F1_conj (double aR, double aI, double c, double x)
              int gsl_sf_hyperg_2F1_conj_e (double aR, double aI, double c, double x, gsl_sf_result * result)

   These routines compute the Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a_R + i a_I, aR - i aI, c, x)

   .. only:: texinfo

      .. math:: 2F1(a_R + i a_I, aR - i aI, c, x)

   with complex parameters for :math:`|x| < 1`.

.. function:: double gsl_sf_hyperg_2F1_renorm (double a, double b, double c, double x)
              int gsl_sf_hyperg_2F1_renorm_e (double a, double b, double c, double x, gsl_sf_result * result)

   These routines compute the renormalized Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a,b,c,x) / \Gamma(c)

   .. only:: texinfo

      .. math:: 2F1(a,b,c,x) / \Gamma(c)

   for :math:`|x| < 1`.
.. exceptions:

.. function:: double gsl_sf_hyperg_2F1_conj_renorm (double aR, double aI, double c, double x)
              int gsl_sf_hyperg_2F1_conj_renorm_e (double aR, double aI, double c, double x, gsl_sf_result * result)

   These routines compute the renormalized Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c)

   .. only:: texinfo

      .. math:: 2F1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c)

   for :math:`|x| < 1`.
.. exceptions:

.. function:: double gsl_sf_hyperg_2F0 (double a, double b, double x)
              int gsl_sf_hyperg_2F0_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the hypergeometric function
   
   .. only:: not texinfo
      
      .. math:: {}_2F_0(a,b,x)

   .. only:: texinfo

      .. math:: 2F0(a,b,x)

   The series representation is a divergent hypergeometric series.
   However, for :math:`x < 0` we have 

   .. only:: not texinfo

      .. math:: {}_2F_0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)

   .. only:: texinfo

      .. math:: 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)

.. exceptions: GSL_EDOM
