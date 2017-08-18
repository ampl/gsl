.. index:: Clausen functions

The Clausen function is defined by the following integral,

.. only:: not texinfo

   .. math:: Cl_2(x) = - \int_0^x dt \log{\left( 2 \sin{(t/2)} \right)}

.. only:: texinfo

   Cl_2(x) = - \int_0^x dt log( 2 \sin(t/2) )

It is related to the :ref:`dilogarithm <dilog-function>` by 
:math:`Cl_2(\theta) = \Im Li_2(\exp(i\theta))`.
The Clausen functions are declared in the header file
:file:`gsl_sf_clausen.h`.

.. function:: double gsl_sf_clausen (double x)
              int gsl_sf_clausen_e (double x, gsl_sf_result * result)

   These routines compute the Clausen integral :math:`Cl_2(x)`.
