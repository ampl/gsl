.. index:: elementary operations
.. index:: multiplication

The following functions allow for the propagation of errors when
combining quantities by multiplication.  The functions are declared in
the header file :file:`gsl_sf_elementary.h`.

.. function:: double gsl_sf_multiply (double x, double y)
              int gsl_sf_multiply_e (double x, double y, gsl_sf_result * result)

   This function multiplies :data:`x` and :data:`y` storing the product and its
   associated error in :data:`result`.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: int gsl_sf_multiply_err_e (double x, double dx, double y, double dy, gsl_sf_result * result)

   This function multiplies :data:`x` and :data:`y` with associated absolute
   errors :data:`dx` and :data:`dy`.  The product 
   :math:`xy \pm xy \sqrt{(dx/x)^2 +(dy/y)^2}`
   is stored in :data:`result`.
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW
