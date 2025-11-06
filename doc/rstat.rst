.. index::
   single: running statistics
   single: online statistics

******************
Running Statistics
******************

This chapter describes routines for computing running statistics,
also known as online statistics, of data. These routines are
suitable for handling large datasets for which it may be
inconvenient or impractical to store in memory all at once.
The data can be processed in a single pass, one point at a time.
Each time a data point is added to the accumulator, internal
parameters are updated in order to compute the current mean, variance,
standard deviation, skewness, and kurtosis. These statistics are
exact, and are updated with numerically stable single-pass algorithms.
The median and arbitrary quantiles are also available, however these
calculations use algorithms which provide approximations, and grow
more accurate as more data is added to the accumulator.

The functions described in this chapter are declared in the header file
:file:`gsl_rstat.h`.

Initializing the Accumulator
============================

.. type:: gsl_rstat_workspace

   This workspace contains parameters used to calculate various statistics
   and are updated after each data point is added to the accumulator.

.. function:: gsl_rstat_workspace * gsl_rstat_alloc (void)

   This function allocates a workspace for computing running statistics.
   The size of the workspace is :math:`O(1)`.

.. function:: void gsl_rstat_free (gsl_rstat_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_rstat_reset (gsl_rstat_workspace * w)

   This function resets the workspace :data:`w` to its initial state,
   so it can begin working on a new set of data.

Adding Data to the Accumulator
==============================

.. function:: int gsl_rstat_add (const double x, gsl_rstat_workspace * w)

   This function adds the data point :data:`x` to the statistical
   accumulator, updating calculations of the mean, variance,
   standard deviation, skewness, kurtosis, and median.

.. function:: size_t gsl_rstat_n (const gsl_rstat_workspace * w)

   This function returns the number of data so far added to the accumulator.

Current Statistics
==================

.. function:: double gsl_rstat_min (const gsl_rstat_workspace * w)

   This function returns the minimum value added to the accumulator.

.. function:: double gsl_rstat_max (const gsl_rstat_workspace * w)

   This function returns the maximum value added to the accumulator.

.. function:: double gsl_rstat_mean (const gsl_rstat_workspace * w)

   This function returns the mean of all data added to the accumulator,
   defined as

   .. only:: not texinfo

      .. math:: {\Hat\mu} = {1 \over N} \sum x_i

   .. only:: texinfo

      ::

         \Hat\mu = (1/N) \sum x_i

.. function:: double gsl_rstat_variance (const gsl_rstat_workspace * w)

   This function returns the variance of all data added to the accumulator,
   defined as

   .. only:: not texinfo

      .. math:: {\Hat\sigma}^2 = {1 \over (N-1)} \sum (x_i - {\Hat\mu})^2

   .. only:: texinfo

      ::

         \Hat\sigma^2 = (1/(N-1)) \sum (x_i - \Hat\mu)^2

.. function:: double gsl_rstat_sd (const gsl_rstat_workspace * w)

   This function returns the standard deviation of all data added to the
   accumulator, defined as the square root of the variance given above.

.. function:: double gsl_rstat_sd_mean (const gsl_rstat_workspace * w)

   This function returns the standard deviation of the mean, defined as

   .. only:: not texinfo
   
      .. math:: \Hat\sigma_{\Hat\mu} = {\Hat\sigma \over \sqrt{N}}

   .. only:: texinfo

      ::

         sd_mean = \Hat\sigma / \sqrt{N}

.. function:: double gsl_rstat_rms (const gsl_rstat_workspace * w)

   This function returns the root mean square of all data added to the
   accumulator, defined as

   .. math:: rms = \sqrt{{1 \over N} \sum x_i^2}

.. function:: double gsl_rstat_skew (const gsl_rstat_workspace * w)

   This function returns the skewness of all data added to the accumulator,
   defined as

   .. only:: not texinfo

      .. math::

         skew = {1 \over N} \sum 
          {\left( x_i - {\Hat\mu} \over {\Hat\sigma} \right)}^3

   .. only:: texinfo

      ::

         skew = (1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^3

.. function:: double gsl_rstat_kurtosis (const gsl_rstat_workspace * w)

   This function returns the kurtosis of all data added to the accumulator,
   defined as

   .. only:: not texinfo

      .. math::

         kurtosis = \left( {1 \over N} \sum 
          {\left(x_i - {\Hat\mu} \over {\Hat\sigma} \right)}^4 
          \right) 
          - 3

   .. only:: texinfo

      ::

         kurtosis = ((1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^4)  - 3

.. function:: double gsl_rstat_median (gsl_rstat_workspace * w)

   This function returns an estimate of the median of the data added to
   the accumulator.

.. function:: double gsl_rstat_norm (const gsl_rstat_workspace * w)

   This function returns the Euclidean norm of all data added to the
   accumulator, defined as

   .. math:: \left|\left| x \right|\right| = \sqrt{\sum x_i^2}

Quantiles
=========

The functions in this section estimate quantiles dynamically without
storing the entire dataset, using the algorithm of Jain and Chlamtec, 1985.
Only five points (markers) are stored which represent the minimum
and maximum of the data, as well as current estimates of the
:math:`p/2`-, :math:`p`-, and :math:`(1+p)/2`-quantiles. Each time
a new data point is added, the marker positions and heights are
updated.

.. type:: gsl_rstat_quantile_workspace

   This workspace contains parameters for estimating quantiles of the
   current dataset

.. function:: gsl_rstat_quantile_workspace * gsl_rstat_quantile_alloc (const double p)

   This function allocates a workspace for the dynamic estimation of
   :data:`p`-quantiles, where :data:`p` is between :math:`0` and :math:`1`.
   The median corresponds to :math:`p = 0.5`. The size of the workspace
   is :math:`O(1)`.

.. function:: void gsl_rstat_quantile_free (gsl_rstat_quantile_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_rstat_quantile_reset (gsl_rstat_quantile_workspace * w)

   This function resets the workspace :data:`w` to its initial state,
   so it can begin working on a new set of data.

.. function:: int gsl_rstat_quantile_add (const double x, gsl_rstat_quantile_workspace * w)

   This function updates the estimate of the :math:`p`-quantile with
   the new data point :data:`x`.

.. function:: double gsl_rstat_quantile_get (gsl_rstat_quantile_workspace * w)

   This function returns the current estimate of the :math:`p`-quantile.

Examples
========

Here is a basic example of how to use the statistical functions:

.. include:: examples/rstat.c
   :code:

The program should produce the following output,

.. include:: examples/rstat.txt
   :code:

This next program estimates the lower quartile, median and upper
quartile from 10,000 samples of a random Rayleigh distribution,
using the :math:`P^2` algorithm of Jain and Chlamtec. For
comparison, the exact values are also computed from the sorted
dataset.

.. include:: examples/rquantile.c
   :code:

The program should produce the following output,

.. include:: examples/rquantile.txt
   :code:

References and Further Reading
==============================

The algorithm used to dynamically estimate :math:`p`-quantiles is described
in the paper,

* R. Jain and I. Chlamtac.
  *The P^2 algorithm for dynamic calculation of quantiles and histograms without storing observations*,
  Communications of the ACM, Volume 28 (October), Number 10, 1985,
  p. 1076-1085.
