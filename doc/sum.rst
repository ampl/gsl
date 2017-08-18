.. index::
   single: acceleration of series
   single: summation, acceleration
   single: series, acceleration
   single: u-transform for series
   single: Levin u-transform
   single: convergence, accelerating a series

*******************
Series Acceleration
*******************

The functions described in this chapter accelerate the convergence of a
series using the Levin :math:`u`-transform.  This method takes a small number of
terms from the start of a series and uses a systematic approximation to
compute an extrapolated value and an estimate of its error.  The
:math:`u`-transform works for both convergent and divergent series, including
asymptotic series.

These functions are declared in the header file :file:`gsl_sum.h`.

Acceleration functions
======================

The following functions compute the full Levin :math:`u`-transform of a series
with its error estimate.  The error estimate is computed by propagating
rounding errors from each term through to the final extrapolation. 

These functions are intended for summing analytic series where each term
is known to high accuracy, and the rounding errors are assumed to
originate from finite precision. They are taken to be relative errors of
order :macro:`GSL_DBL_EPSILON` for each term.

The calculation of the error in the extrapolated value is an
:math:`O(N^2)` process, which is expensive in time and memory.  A faster
but less reliable method which estimates the error from the convergence
of the extrapolated value is described in the next section.  For the
method described here a full table of intermediate values and
derivatives through to :math:`O(N)` must be computed and stored, but this
does give a reliable error estimate.

.. type:: gsl_sum_levin_u_workspace

   Workspace for a Leven :math:`u`-transform.

.. function:: gsl_sum_levin_u_workspace * gsl_sum_levin_u_alloc (size_t n)

   This function allocates a workspace for a Levin :math:`u`-transform of :data:`n`
   terms.  The size of the workspace is :math:`O(2n^2 + 3n)`.

.. function:: void gsl_sum_levin_u_free (gsl_sum_levin_u_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_sum_levin_u_accel (const double * array, size_t array_size, gsl_sum_levin_u_workspace * w, double * sum_accel, double * abserr)

   This function takes the terms of a series in :data:`array` of size
   :data:`array_size` and computes the extrapolated limit of the series using
   a Levin :math:`u`-transform.  Additional working space must be provided in
   :data:`w`.  The extrapolated sum is stored in :data:`sum_accel`, with an
   estimate of the absolute error stored in :data:`abserr`.  The actual
   term-by-term sum is returned in :code:`w->sum_plain`. The algorithm
   calculates the truncation error (the difference between two successive
   extrapolations) and round-off error (propagated from the individual
   terms) to choose an optimal number of terms for the extrapolation.  
   All the terms of the series passed in through :data:`array` should be non-zero.

Acceleration functions without error estimation
===============================================

The functions described in this section compute the Levin :math:`u`-transform of
series and attempt to estimate the error from the "truncation error" in
the extrapolation, the difference between the final two approximations.
Using this method avoids the need to compute an intermediate table of
derivatives because the error is estimated from the behavior of the
extrapolated value itself. Consequently this algorithm is an :math:`O(N)`
process and only requires :math:`O(N)` terms of storage.  If the series
converges sufficiently fast then this procedure can be acceptable.  It
is appropriate to use this method when there is a need to compute many
extrapolations of series with similar convergence properties at high-speed.
For example, when numerically integrating a function defined by a
parameterized series where the parameter varies only slightly. A
reliable error estimate should be computed first using the full
algorithm described above in order to verify the consistency of the
results.

.. type:: gsl_sum_levin_utrunc_workspace

   Workspace for a Levin :math:`u`-transform without error estimation

.. function:: gsl_sum_levin_utrunc_workspace * gsl_sum_levin_utrunc_alloc (size_t n)

   This function allocates a workspace for a Levin :math:`u`-transform of :data:`n`
   terms, without error estimation.  The size of the workspace is
   :math:`O(3n)`.

.. function:: void gsl_sum_levin_utrunc_free (gsl_sum_levin_utrunc_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_sum_levin_utrunc_accel (const double * array, size_t array_size, gsl_sum_levin_utrunc_workspace * w, double * sum_accel, double * abserr_trunc)

   This function takes the terms of a series in :data:`array` of size
   :data:`array_size` and computes the extrapolated limit of the series using
   a Levin :math:`u`-transform.  Additional working space must be provided in
   :data:`w`.  The extrapolated sum is stored in :data:`sum_accel`.  The actual
   term-by-term sum is returned in :code:`w->sum_plain`. The algorithm
   terminates when the difference between two successive extrapolations
   reaches a minimum or is sufficiently small. The difference between these
   two values is used as estimate of the error and is stored in
   :data:`abserr_trunc`.  To improve the reliability of the algorithm the
   extrapolated values are replaced by moving averages when calculating the
   truncation error, smoothing out any fluctuations.

Examples
========

The following code calculates an estimate of :math:`\zeta(2) = \pi^2 / 6`
using the series,

.. math:: \zeta(2) = 1 + 1/2^2 + 1/3^2 + 1/4^2 + \dots

After :data:`N` terms the error in the sum is :math:`O(1/N)`, making direct
summation of the series converge slowly.

.. include:: examples/sum.c
   :code:

The output below shows that the Levin :math:`u`-transform is able to obtain an 
estimate of the sum to 1 part in 
:math:`10^{10}`
using the first eleven terms of the series.  The
error estimate returned by the function is also accurate, giving
the correct number of significant digits. 

.. include:: examples/sum.txt
   :code:

Note that a direct summation of this series would require 
:math:`10^{10}`
terms to achieve the same precision as the accelerated 
sum does in 13 terms.

References and Further Reading
==============================

The algorithms used by these functions are described in the following papers,

* T. Fessler, W.F. Ford, D.A. Smith,
  HURRY: An acceleration algorithm for scalar sequences and series
  *ACM Transactions on Mathematical Software*, 9(3):346--354, 1983.
  and Algorithm 602 9(3):355--357, 1983.

The theory of the :math:`u`-transform was presented by Levin,

* D. Levin,
  Development of Non-Linear Transformations for Improving Convergence of
  Sequences, *Intern.: J.: Computer Math.* B3:371--388, 1973.

A review paper on the Levin Transform is available online,

* Herbert H. H. Homeier, Scalar Levin-Type Sequence Transformations,
  http://arxiv.org/abs/math/0005209
