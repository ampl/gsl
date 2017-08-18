.. index::
   single: psi function
   single: digamma function
   single: polygamma functions

The polygamma functions of order :math:`n` are defined by

.. only:: not texinfo

   .. math:: \psi^{(n)}(x) = \left(d \over dx\right)^n \psi(x) = \left(d \over dx\right)^{n+1} \log(\Gamma(x))

.. only:: texinfo

   ::

      \psi^{(n)}(x) = (d/dx)^n \psi(x) = (d/dx)^{n+1} \log(\Gamma(x))

where :math:`\psi(x) = \Gamma'(x)/\Gamma(x)` is known as the digamma function.
These functions are declared in the header file :file:`gsl_sf_psi.h`.

Digamma Function
----------------

.. function:: double gsl_sf_psi_int (int n)
              int gsl_sf_psi_int_e (int n, gsl_sf_result * result)

   These routines compute the digamma function :math:`\psi(n)` for positive
   integer :data:`n`.  The digamma function is also called the Psi function.
.. Domain: n integer, n > 0
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_psi (double x)
              int gsl_sf_psi_e (double x, gsl_sf_result * result)

   These routines compute the digamma function :math:`\psi(x)` for general
   :data:`x`, :math:`x \ne 0`.
.. Domain: x != 0.0, -1.0, -2.0, ...
.. Exceptional Return Values: GSL_EDOM, GSL_ELOSS

.. function:: double gsl_sf_psi_1piy (double y)
              int gsl_sf_psi_1piy_e (double y, gsl_sf_result * result)

   These routines compute the real part of the digamma function on the line
   :math:`1 + i y`, :math:`\Re[\psi(1 + i y)]`.
.. exceptions: none
.. Exceptional Return Values: none

Trigamma Function
-----------------

.. function:: double gsl_sf_psi_1_int (int n)
              int gsl_sf_psi_1_int_e (int n, gsl_sf_result * result)

   These routines compute the Trigamma function :math:`\psi'(n)` for
   positive integer :math:`n`.
.. Domain: n integer, n > 0 
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_psi_1 (double x)
              int gsl_sf_psi_1_e (double x, gsl_sf_result * result)

   These routines compute the Trigamma function :math:`\psi'(x)` for
   general :data:`x`.
.. Domain: x != 0.0, -1.0, -2.0, ...
.. Exceptional Return Values: GSL_EDOM, GSL_ELOSS

Polygamma Function
------------------

.. function:: double gsl_sf_psi_n (int n, double x)
              int gsl_sf_psi_n_e (int n, double x, gsl_sf_result * result)

   These routines compute the polygamma function :math:`\psi^{(n)}(x)` for
   :math:`n \ge 0`, :math:`x > 0`.
.. Domain: n >= 0, x > 0.0
.. Exceptional Return Values: GSL_EDOM
