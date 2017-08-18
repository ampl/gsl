.. index::
   single: discrete Hankel transforms
   single: Hankel transforms, discrete
   single: transforms, Hankel

**************************
Discrete Hankel Transforms
**************************

This chapter describes functions for performing Discrete Hankel
Transforms (DHTs). The functions are declared in the header file
:file:`gsl_dht.h`.

Definitions
===========

The discrete Hankel transform acts on a vector of sampled data, where
the samples are assumed to have been taken at points related to the
zeros of a Bessel function of fixed order; compare this to the case of
the discrete Fourier transform, where samples are taken at points
related to the zeroes of the sine or cosine function.

Starting with its definition, the Hankel transform (or Bessel transform) of
order :math:`\nu` of a function :math:`f` with :math:`\nu > -1/2` is defined as
(see Johnson, 1987 and Lemoine, 1994)

.. math:: F_\nu(u) = \int_0^\infty f(t) J_\nu(u t) t dt

If the integral exists, :math:`F_\nu` is called the Hankel transformation
of :math:`f`. The reverse transform is given by

.. math:: f(t) = \int_0^\infty F_\nu(u) J_\nu(u t) u du

where :math:`\int_0^\infty f(t) t^{1/2} dt` must exist and be
absolutely convergent, and where :math:`f(t)` satisfies Dirichlet's
conditions (of limited total fluctuations) in the interval
:math:`[0,\infty]`.

Now the discrete Hankel transform works on a discrete function
:math:`f`, which is sampled on points :math:`n=1...M` located at
positions :math:`t_n=(j_{\nu,n}/j_{\nu,M}) X` in real space and
at :math:`u_n=j_{\nu,n}/X` in reciprocal space. Here,
:math:`j_{\nu,m}` are the m-th zeros of the Bessel function
:math:`J_\nu(x)` arranged in ascending order. Moreover, the
discrete functions are assumed to be band limited, so
:math:`f(t_n)=0` and :math:`F(u_n)=0` for :math:`n>M`. Accordingly,
the function :math:`f` is defined on the interval :math:`[0,X]`.

Following the work of Johnson, 1987 and
Lemoine, 1994, the discrete Hankel transform is given by

.. only:: not texinfo

   .. math::

      F_\nu(u_m) = {{2 X^2}\over{j_{\nu,M}^2}}
            \sum_{k=1}^{M-1} f\left({{j_{\nu,k} X}\over{j_{\nu,M}}}\right)
                {{J_\nu(j_{\nu,m} j_{\nu,k} / j_{\nu,M})}\over{J_{\nu+1}(j_{\nu,k})^2}}.

.. only:: texinfo

   ::

      F_\nu(u_m) = (2 X^2 / j_(\nu,M)^2)
            \sum_{k=1}^{M-1} f(j_(\nu,k) X/j_(\nu,M))
                (J_\nu(j_(\nu,m) j_(\nu,k) / j_(\nu,M)) / J_(\nu+1)(j_(\nu,k))^2).

It is this discrete expression which defines the discrete Hankel
transform calculated by GSL. In GSL, forward and backward transforms
are defined equally and calculate :math:`F_\nu(u_m)`.
Following Johnson, the backward transform reads

.. only:: not texinfo

   .. math::

      f(t_k) = {{2}\over{X^2}}
            \sum_{m=1}^{M-1} F\left({{j_{\nu,m}}\over{X}}\right)
                {{J_\nu(j_{\nu,m} j_{\nu,k} / j_{\nu,M})}\over{J_{\nu+1}(j_{\nu,m})^2}}.

.. only:: texinfo

   ::

      f(t_k) = (2 / X^2)
            \sum_{m=1}^{M-1} F(j_(\nu,m)/X)
                (J_\nu(j_(\nu,m) j_(\nu,k) / j_(\nu,M)) / J_(\nu+1)(j_(\nu,m))^2).

Obviously, using the forward transform instead of the backward transform gives an
additional factor :math:`X^4/j_{\nu,M}^2=t_m^2/u_m^2`.

The kernel in the summation above defines the matrix of the
:math:`\nu`-Hankel transform of size :math:`M-1`. The coefficients of
this matrix, being dependent on :math:`\nu` and :math:`M`, must be
precomputed and stored; the :type:`gsl_dht` object encapsulates this
data. The allocation function :func:`gsl_dht_alloc` returns a
:type:`gsl_dht` object which must be properly initialized with
:func:`gsl_dht_init` before it can be used to perform transforms on data
sample vectors, for fixed :math:`\nu` and :math:`M`, using the
:func:`gsl_dht_apply` function. The implementation allows to define the
length :math:`X` of the fundamental interval, for convenience, while
discrete Hankel transforms are often defined on the unit interval
instead of :math:`[0,X]`.

Notice that by assumption :math:`f(t)` vanishes at the endpoints
of the interval, consistent with the inversion formula
and the sampling formula given above. Therefore, this transform
corresponds to an orthogonal expansion in eigenfunctions
of the Dirichlet problem for the Bessel differential equation.

Functions
=========

.. type:: gsl_dht

   Workspace for computing discrete Hankel transforms

.. function:: gsl_dht * gsl_dht_alloc (size_t size)

   This function allocates a Discrete Hankel transform object of size
   :data:`size`.

.. function:: int gsl_dht_init (gsl_dht * t, double nu, double xmax)

   This function initializes the transform :data:`t` for the given values of
   :data:`nu` and :data:`xmax`.

.. function:: gsl_dht * gsl_dht_new (size_t size, double nu, double xmax)

   This function allocates a Discrete Hankel transform object of size
   :data:`size` and initializes it for the given values of :data:`nu` and
   :data:`xmax`.

.. function:: void gsl_dht_free (gsl_dht * t)

   This function frees the transform :data:`t`.

.. function:: int gsl_dht_apply (const gsl_dht * t, double * f_in, double * f_out)

   This function applies the transform :data:`t` to the array :data:`f_in`
   whose size is equal to the size of the transform.  The result is stored
   in the array :data:`f_out` which must be of the same length.   

   Applying this function to its output gives the original data
   multiplied by :math:`(X^2/j_{\nu,M})^2`,
   up to numerical errors.

.. function:: double gsl_dht_x_sample (const gsl_dht * t, int n)

   This function returns the value of the :data:`n`-th sample point in the unit interval,
   :math:`{({j_{\nu,n+1}} / {j_{\nu,M}}}) X`.
   These are the points where the function :math:`f(t)` is assumed to be sampled.

.. function:: double gsl_dht_k_sample (const gsl_dht * t, int n)

   This function returns the value of the :data:`n`-th sample point in "k-space",
   :math:`{{j_{\nu,n+1}} / X}`.

References and Further Reading
==============================

The algorithms used by these functions are described in the following papers,

* H. Fisk Johnson, Comp.: Phys.: Comm.: 43, 181 (1987).

* D. Lemoine, J. Chem.: Phys.: 101, 3936 (1994).
