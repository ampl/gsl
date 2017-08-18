.. index::
   single: Jacobi elliptic functions
   single: elliptic functions (Jacobi)

The Jacobian Elliptic functions are defined in Abramowitz & Stegun,
Chapter 16.  The functions are declared in the header file
:file:`gsl_sf_elljac.h`.

.. function:: int gsl_sf_elljac_e (double u, double m, double * sn, double * cn, double * dn)

   This function computes the Jacobian elliptic functions :math:`sn(u|m)`,
   :math:`cn(u|m)`, :math:`dn(u|m)` by descending Landen
   transformations.
.. Exceptional Return Values: GSL_EDOM
