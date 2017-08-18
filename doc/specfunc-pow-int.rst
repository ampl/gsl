.. index::
   single: power function
   single: integer powers

The following functions are equivalent to the function :func:`gsl_pow_int`
with an error estimate.  These functions are
declared in the header file :file:`gsl_sf_pow_int.h`.

.. function:: double gsl_sf_pow_int (double x, int n)
              int gsl_sf_pow_int_e (double x, int n, gsl_sf_result * result)

   These routines compute the power :math:`x^n` for integer :data:`n`.  The
   power is computed using the minimum number of multiplications. For
   example, :math:`x^8` is computed as :math:`((x^2)^2)^2`, requiring only 3
   multiplications.  For reasons of efficiency, these functions do not
   check for overflow or underflow conditions. The following is a simple example::

      #include <gsl/gsl_sf_pow_int.h>
      /* compute 3.0**12 */
      double y = gsl_sf_pow_int(3.0, 12); 
