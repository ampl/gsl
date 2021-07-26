.. index::
   single: moving window statistics
   single: statistics, moving window
   single: online statistics

************************
Moving Window Statistics
************************

This chapter describes routines for computing *moving window
statistics* (also called rolling statistics and running statistics),
using a window around a sample which is used to calculate various
local statistical properties of an input data stream. The window is
then slid forward by one sample to process the next data point and so on.

The functions described in this chapter are declared in the header file
:file:`gsl_movstat.h`.

Introduction
============

This chapter is concerned with calculating various statistics from
subsets of a given dataset. The main idea is to compute statistics
in the vicinity of a given data sample by defining a *window* which
includes the sample itself as well as some specified number of samples
before and after the sample in question. For a sample :math:`x_i`, we
define a window :math:`W_i^{H,J}` as

.. only:: not texinfo

   .. math:: W_i^{H,J} = \left\{ x_{i-H}, \dots, x_i, \dots, x_{i+J} \right\}

.. only:: texinfo

   ::

      W_i^{H,J} = {x_{i-H},...,x_i,...,x_{i+J}}

The parameters :math:`H` and :math:`J` are non-negative integers specifying
the number of samples to include before and after the sample :math:`x_i`.
Statistics such as the mean and standard deviation of the window :math:`W_i^{H,J}`
may be computed, and then the window is shifted forward by one sample to
focus on :math:`x_{i+1}`. The total number of samples in the window is
:math:`K = H + J + 1`. To define a symmetric window centered on :math:`x_i`,
one would set :math:`H = J = \left\lfloor K / 2 \right\rfloor`.

Handling Endpoints
==================

When processing samples near the ends of the input signal, there will not
be enough samples to fill the window :math:`W_i^{H,J}` defined above.
Therefore the user must specify how to construct the windows near the end points.
This is done by passing an input argument of type :type:`gsl_movstat_end_t`:

.. type:: gsl_movstat_end_t

   This data type specifies how to construct windows near end points and can
   be selected from the following choices:

   .. macro:: GSL_MOVSTAT_END_PADZERO

      With this option, a full window of length :math:`K` will be constructed
      by inserting zeros into the window near the signal end points. Effectively,
      the input signal is modified to

      .. only:: not texinfo

         .. math:: \tilde{x} = \{ \underbrace{0, \dots, 0}_{H \textrm{ zeros}}, x_1, x_2, \dots, x_{n-1}, x_n, \underbrace{0, \dots, 0}_{J \textrm{ zeros} } \}

      .. only:: texinfo

         ::

            x~ = {0, ..., 0, x_1, x_2, ..., x_{n-1}, x_n, 0, ..., 0}

      to ensure a well-defined window for all :math:`x_i`.

   .. macro:: GSL_MOVSTAT_END_PADVALUE

      With this option, a full window of length :math:`K` will be constructed
      by padding the window with the first and last sample in the input signal.
      Effectively, the input signal is modified to

      .. only:: not texinfo

         .. math:: \tilde{x} = \{ \underbrace{x_1, \dots, x_1}_{H}, x_1, x_2, \dots, x_{n-1}, x_n, \underbrace{x_n, \dots, x_n}_{J} \}

      .. only:: texinfo

         ::

            x~ = {x_1, ..., x_1, x_1, x_2, ..., x_{n-1}, x_n, x_n, ..., x_n}

   .. macro:: GSL_MOVSTAT_END_TRUNCATE

      With this option, no padding is performed, and the windows are simply truncated
      as the end points are approached.

.. index::
   single: moving window, allocation

Allocation for Moving Window Statistics
=======================================

.. type:: gsl_movstat_workspace

   The moving window statistical routines use a common workspace.

.. function:: gsl_movstat_workspace * gsl_movstat_alloc(const size_t K)

   This function allocates a workspace for computing symmetric, centered moving statistics with a window
   length of :math:`K` samples. In this case, :math:`H = J = \left\lfloor K/2 \right\rfloor`. The size of the workspace
   is :math:`O(7K)`.

.. function:: gsl_movstat_workspace * gsl_movstat_alloc2(const size_t H, const size_t J)

   This function allocates a workspace for computing moving statistics using a window with :math:`H`
   samples prior to the current sample, and :math:`J` samples after the current sample. The
   total window size is :math:`K = H + J + 1`. The size of the workspace is :math:`O(7K)`.

.. function:: void * gsl_movstat_free(gsl_movstat_workspace * w)

   This function frees the memory associated with :data:`w`.

.. index::
   single: moving mean
   single: rolling mean

Moving Mean
===========

The moving window mean calculates the mean of the values of each window :math:`W_i^{H,J}`.

.. only:: not texinfo

   .. math:: \hat{\mu}_i = \frac{1}{\left| W_i^{H,J} \right|} \sum_{x_m \in W_i^{H,J}} x_m

.. only:: texinfo

   ::

      \hat{\mu}_i = 1/| W_i^{H,J} | \sum_{x_m \in W_i^{H,J}} x_m

Here, :math:`\left| W_i^{H,J} \right|` represents the number of elements in the window
:math:`W_i^{H,J}`. This will normally be :math:`K`, unless the :macro:`GSL_MOVSTAT_END_TRUNCATE`
option is selected, in which case it could be less than :math:`K` near the signal end points.

.. function:: int gsl_movstat_mean(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)

   This function computes the moving window mean of the input vector :data:`x`, storing
   the output in :data:`y`. The parameter :data:`endtype` specifies how windows near
   the ends of the input should be handled. It is allowed to have :data:`x` = :data:`y`
   for an in-place moving mean.

.. index::
   single: moving variance
   single: moving standard deviation
   single: rolling variance
   single: rolling standard deviation

Moving Variance and Standard Deviation
======================================

The moving window variance calculates the *sample variance* of the values of each window :math:`W_i^{H,J}`,
defined by

.. only:: not texinfo

   .. math:: \hat{\sigma}_i^2 = \frac{1}{\left( \left| W_i^{H,J} \right| - 1 \right)} \sum_{x_m \in W_i^{H,J}} \left( x_m - \hat{\mu}_i \right)^2

.. only:: texinfo

   ::

      \hat{\sigma}_i^2 = 1/(|W_i^{H,J}| - 1) \sum_{x_m \in W_i^{H,J}} ( x_m - \hat{\mu}_i )^2

where :math:`\hat{\mu}_i` is the mean of :math:`W_i^{H,J}` defined above. The standard deviation :math:`\hat{\sigma}_i`
is the square root of the variance.

.. function:: int gsl_movstat_variance(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)

   This function computes the moving window variance of the input vector :data:`x`, storing
   the output in :data:`y`. The parameter :data:`endtype` specifies how windows near
   the ends of the input should be handled. It is allowed to have :data:`x` = :data:`y`
   for an in-place moving variance.

.. function:: int gsl_movstat_sd(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)

   This function computes the moving window standard deviation of the input vector :data:`x`, storing
   the output in :data:`y`. The parameter :data:`endtype` specifies how windows near
   the ends of the input should be handled. It is allowed to have :data:`x` = :data:`y`
   for an in-place moving standard deviation.

.. index::
   single: moving minimum
   single: moving maximum
   single: rolling minimum
   single: rolling maximum

Moving Minimum and Maximum
==========================

The moving minimum/maximum calculates the minimum and maximum values of
each window :math:`W_i^{H,J}`.

.. only:: not texinfo

    .. math::

       y_i^{min} &= \min \left( W_i^{H,J} \right) \\
       y_i^{max} &= \max \left( W_i^{H,J} \right)

.. only:: texinfo

   ::

      y_i^{min} = \min W_i^{H,J}
      y_i^{max} = \max W_i^{H,J}

.. function:: int gsl_movstat_min(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)

   This function computes the moving minimum of the input vector :data:`x`, storing
   the result in :data:`y`. The parameter :data:`endtype` specifies how windows near
   the ends of the input should be handled.
   It is allowed to have :data:`x` = :data:`y` for an in-place moving minimum.

.. function:: int gsl_movstat_max(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)

   This function computes the moving maximum of the input vector :data:`x`, storing
   the result in :data:`y`. The parameter :data:`endtype` specifies how windows near
   the ends of the input should be handled.
   It is allowed to have :data:`x` = :data:`y` for an in-place moving maximum.

.. function:: int gsl_movstat_minmax(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y_min, gsl_vector * y_max, gsl_movstat_workspace * w)

   This function computes the moving minimum and maximum of the input vector :data:`x`, storing
   the window minimums in :data:`y_min` and the window maximums in :data:`y_max`.
   The parameter :data:`endtype` specifies how windows near the ends of the input should be handled.

.. index::
   single: moving sum
   single: rolling sum

Moving Sum
==========

The moving window sum calculates the sum of the values of each window :math:`W_i^{H,J}`.

.. math:: y_i = \sum_{x_m \in W_i^{H,J}} x_m

.. function:: int gsl_movstat_sum(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)

   This function computes the moving window sum of the input vector :data:`x`, storing
   the output in :data:`y`. The parameter :data:`endtype` specifies how windows near
   the ends of the input should be handled. It is allowed to have :data:`x` = :data:`y`
   for an in-place moving sum.

.. index::
   single: moving median
   single: rolling median

Moving Median
=============

The moving median calculates the median of the window :math:`W_i^{H,J}` for
each sample :math:`x_i`:

.. only:: not texinfo

    .. math:: y_i = \textrm{median} \left( W_i^{H,J} \right)

.. only:: texinfo

   ::

      y_i = median(W_i^{H,J})

.. function:: int gsl_movstat_median(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)

   This function computes the moving median of the input vector :data:`x`, storing
   the output in :data:`y`. The parameter :data:`endtype` specifies how windows near
   the ends of the input should be handled. It is allowed for
   :data:`x` = :data:`y` for an in-place moving window median.

Robust Scale Estimation
=======================

A common problem in statistics is to quantify the dispersion (also known as the variability, scatter, and spread) of
a set of data. Often this is done by calculating the variance or standard deviation. However these statistics
are strongly influenced by outliers, and can often provide erroneous results when even a small number of outliers
are present.

Several useful statistics have emerged to provide robust estimates of scale which are not as susceptible to data outliers.
A few of these statistical scale estimators are described below.

.. index::
   single: moving median absolute deviation
   single: rolling median absolute deviation

Moving MAD
----------

The median absolute deviation (MAD) for the window :math:`W_i^{H,J}` is defined
to be the median of the absolute deviations from the window's median:

.. only:: not texinfo

    .. math:: MAD_i = 1.4826 \times \textrm{median} \left( \left| W_i^{H,J} - \textrm{median} \left( W_i^{H,J} \right) \right| \right)

.. only:: texinfo

   ::

      MAD_i = 1.4826 * median[ |W_i^{H,J} - median(W_i^{H,J})| ]

The factor of :math:`1.4826` makes the MAD an unbiased estimator of the standard deviation
for Gaussian data. The MAD has an efficiency of 37%.  See :ref:`here <sec_mad-statistic>` for more information.

.. function:: int gsl_movstat_mad0(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * xmedian, gsl_vector * xmad, gsl_movstat_workspace * w)
              int gsl_movstat_mad(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * xmedian, gsl_vector * xmad, gsl_movstat_workspace * w)

   These functions compute the moving MAD of the input vector :data:`x` and store the result
   in :data:`xmad`. The medians of each window :math:`W_i^{H,J}` are stored in :data:`xmedian`
   on output. The inputs :data:`x`, :data:`xmedian`, and :data:`xmad` must all be the same length.
   The parameter :data:`endtype` specifies how windows near the ends of the input should be handled.
   The function :code:`mad0` does not include the scale factor of :math:`1.4826`, while the
   function :code:`mad` does include this factor.

.. index::
   single: moving quantile range
   single: rolling quantile range

Moving QQR
----------

The q-quantile range (QQR) is the difference between the :math:`(1-q)` and :math:`q` quantiles
of a set of data,

.. math:: QQR = Q_{1-q} - Q_q

The case :math:`q = 0.25` corresponds to the well-known *interquartile range (IQR)*, which
is the difference between the 75th and 25th percentiles of a set of data. The QQR is
a *trimmed estimator*, the main idea being to discard the largest and smallest values in
a data window and compute a scale estimate from the remaining middle values. In the case
of the IQR, the largest and smallest 25% of the data are discarded and the scale is
estimated from the remaining (middle) 50%.

.. function:: int gsl_movstat_qqr(const gsl_movstat_end_t endtype, const gsl_vector * x, const double q, gsl_vector * xqqr, gsl_movstat_workspace * w)

   This function computes the moving QQR of the input vector :data:`x` and stores the q-quantile ranges
   of each window :math:`W_i^{H,J}` in :data:`xqqr`. The quantile parameter :data:`q` must be between
   :math:`0` and :math:`0.5`. The input :math:`q = 0.25` corresponds to the IQR.
   The inputs :data:`x` and :data:`xqqr` must be the same length.
   The parameter :data:`endtype` specifies how windows near the ends of the input should be handled.

Moving :math:`S_n`
------------------

The :math:`S_n` statistic proposed by Croux and Rousseeuw is based on pairwise differences between
all samples in the window. It has an efficiency of 58%, significantly higher than the MAD.
See :ref:`here <sec_Sn-statistic>` for more information.

.. function:: int gsl_movstat_Sn(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * xscale, gsl_movstat_workspace * w)

   This function computes the moving :math:`S_n` of the input vector :data:`x` and stores the output
   in :data:`xscale`. The inputs :data:`x` and :data:`xscale` must be the same length.
   The parameter :data:`endtype` specifies how windows near the ends of the input should be handled.
   It is allowed for :data:`x` = :data:`xscale` for an in-place moving window :math:`S_n`.

Moving :math:`Q_n`
------------------

The :math:`Q_n` statistic proposed by Croux and Rousseeuw is loosely based on the Hodges-Lehmann location
estimator. It has a relatively high efficiency of 82%. See :ref:`here <sec_Qn-statistic>` for more information.

.. function:: int gsl_movstat_Qn(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * xscale, gsl_movstat_workspace * w)

   This function computes the moving :math:`Q_n` of the input vector :data:`x` and stores the output
   in :data:`xscale`. The inputs :data:`x` and :data:`xscale` must be the same length.
   The parameter :data:`endtype` specifies how windows near the ends of the input should be handled.
   It is allowed for :data:`x` = :data:`xscale` for an in-place moving window :math:`Q_n`.

.. index::
   single: moving window accumulators
   single: rolling window accumulators

User-defined Moving Statistics
==============================

GSL offers an interface for users to define their own moving window statistics
functions, without needing to implement the edge-handling and accumulator
machinery. This can be done by explicitly constructing the windows
:math:`W_i^{H,J}` for a given input signal (:func:`gsl_movstat_fill`), or by calculating a user-defined
function for each window automatically. In order to apply a user-defined
function to each window, users must define a variable of type
:type:`gsl_movstat_function` to pass into :func:`gsl_movstat_apply`.
This structure is defined as follows.

.. type:: gsl_movstat_function

   Structure specifying user-defined moving window statistical function::

     typedef struct
     {
       double (* function) (const size_t n, double x[], void * params);
       void * params;
     } gsl_movstat_function;

   This structure contains a pointer to the user-defined function as well
   as possible parameters to pass to the function.

   .. member:: double (* function) (const size_t n, double x[], void * params)

      This function returns the user-defined statistic of the array :data:`x`
      of length :data:`n`. User-specified parameters are passed in via :data:`params`.
      It is allowed to modify the array :data:`x`.

   .. member:: void * params

      User-specified parameters to be passed into the function.

.. function:: int gsl_movstat_apply(const gsl_movstat_end_t endtype, const gsl_movstat_function * F, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)

   This function applies the user-defined moving window statistic specified in :data:`F`
   to the input vector :data:`x`, storing the output in :data:`y`.
   The parameter :data:`endtype` specifies how windows near the ends of the input should be handled.
   It is allowed for :data:`x` = :data:`y` for an in-place moving window calculation.

.. function:: size_t gsl_movstat_fill(const gsl_movstat_end_t endtype, const gsl_vector * x, const size_t idx, const size_t H, const size_t J, double * window)

   This function explicitly constructs the sliding window for the input vector :data:`x` which
   is centered on the sample :data:`idx`. On output, the array :data:`window` will contain
   :math:`W_{idx}^{H,J}`. The number of samples to the left and right
   of the sample :data:`idx` are specified by :data:`H` and :data:`J` respectively.
   The parameter :data:`endtype` specifies how windows near the ends of the input should be handled.
   The function returns the size of the window.

Accumulators
============

Many of the algorithms of this chapter are based on an accumulator design, which
process the input vector one sample at a time, updating calculations of the
desired statistic for the current window. Each accumulator is stored in the
following structure:

.. type:: gsl_movstat_accum

   Structure specifying accumulator for moving window statistics::

     typedef struct
     {
       size_t (* size) (const size_t n);
       int (* init) (const size_t n, void * vstate);
       int (* insert) (const double x, void * vstate);
       int (* delete) (void * vstate);
       int (* get) (void * params, double * result, const void * vstate);
     } gsl_movstat_accum;

   The structure contains function pointers responsible for performing
   different tasks for the accumulator.

   .. member:: size_t (* size) (const size_t n)

        This function returns the size of the workspace (in bytes) needed by the accumulator
        for a moving window of length :data:`n`.

   .. member:: int (* init) (const size_t n, void * vstate)

        This function initializes the workspace :data:`vstate` for a moving window of length :data:`n`.

   .. member:: int (* insert) (const double x, void * vstate)

        This function inserts a single sample :data:`x` into the accumulator, updating internal
        calculations of the desired statistic. If the accumulator is full (i.e. :math:`n` samples
        have already been inserted), then the oldest sample is deleted from the accumulator.

   .. member:: int (* delete) (void * vstate)

        This function deletes the oldest sample from the accumulator, updating internal
        calculations of the desired statistic.

   .. member:: int (* get) (void * params, double * result, const void * vstate)

        This function stores the desired statistic for the current window in
        :data:`result`. The input :data:`params` specifies optional parameters
        for calculating the statistic.

The following accumulators of type :type:`gsl_movstat_accum` are defined by GSL to perform moving window statistics
calculations.

.. var:: gsl_movstat_accum * gsl_movstat_accum_min
         gsl_movstat_accum * gsl_movstat_accum_max
         gsl_movstat_accum * gsl_movstat_accum_minmax

   These accumulators calculate moving window minimum/maximums efficiently, using
   the algorithm of D. Lemire.

.. var:: gsl_movstat_accum * gsl_movstat_accum_mean
         gsl_movstat_accum * gsl_movstat_accum_sd
         gsl_movstat_accum * gsl_movstat_accum_variance

   These accumulators calculate the moving window mean, standard deviation, and variance,
   using the algorithm of B. P. Welford.

.. var:: gsl_movstat_accum * gsl_movstat_accum_median

   This accumulator calculates the moving window median using the min/max heap algorithm
   of Härdle and Steiger.

.. var:: gsl_movstat_accum * gsl_movstat_accum_Sn
         gsl_movstat_accum * gsl_movstat_accum_Qn

   These accumulators calculate the moving window :math:`S_n` and :math:`Q_n` statistics
   developed by Croux and Rousseeuw.

.. var:: gsl_movstat_accum * gsl_movstat_accum_sum

   This accumulator calculates the moving window sum.

.. var:: gsl_movstat_accum * gsl_movstat_accum_qqr

   This accumulator calculates the moving window q-quantile range.

Examples
========

Example 1
---------

The following example program computes the moving mean, minimum and maximum of a noisy
sinusoid signal of length :math:`N = 500` with a symmetric moving window of size :math:`K = 11`.

.. _fig_movstat1:

.. figure:: /images/movstat1.png
   :scale: 60%

   Original signal time series (gray) with moving mean (green), moving minimum (blue),
   and moving maximum (orange).

The program is given below.

.. include:: examples/movstat1.c
   :code:

Example 2: Robust Scale
-----------------------

The following example program analyzes a time series of length :math:`N = 1000` composed
of Gaussian random variates with zero mean whose standard deviation changes in a piecewise constant fashion
as shown in the table below.

============ ==============
Sample Range :math:`\sigma`
============ ==============
1-200        1.0
201-450      5.0
451-600      1.0
601-850      3.0
851-1000     5.0
============ ==============

Additionally, about 1% of the samples are perturbed to represent outliers by adding
:math:`\pm 15` to the random Gaussian variate.
The program calculates the moving statistics MAD, IQR, :math:`S_n`, :math:`Q_n`, and
the standard deviation using a symmetric moving window of length :math:`K = 41`. The results are shown in
:numref:`fig_movstat2`.

.. _fig_movstat2:

.. figure:: /images/movstat2.png
   :scale: 60%

   Top: time series of piecewise constant variance. Bottom: scale estimates using a moving
   window; the true sigma value is in light blue, MAD in green, IQR in red, :math:`S_n` in yellow, and
   :math:`Q_n` in dark blue. The moving standard deviation is shown in gray.

The robust statistics follow the true standard deviation piecewise changes well, without being
influenced by the outliers. The moving standard deviation (gray curve) is heavily influenced by
the presence of the outliers. The program is given below.

.. include:: examples/movstat2.c
   :code:

Example 3: User-defined Moving Window
-------------------------------------

This example program illustrates how a user can define their own moving window function to apply
to an input vector. It constructs a random noisy time series of length :math:`N = 1000` with
some outliers added. Then it applies a moving window trimmed mean to the time series with
trim parameter :math:`\alpha = 0.1`. The length of the moving window is :math:`K = 11`, so
the smallest and largest sample of each window is discarded prior to computing the mean.
The results are shown in :numref:`fig_movstat3`.

.. _fig_movstat3:

.. figure:: /images/movstat3.png
   :scale: 60%

   Noisy time series data (black) with moving window trimmed mean (red)

The program is given below.

.. include:: examples/movstat3.c
   :code:

References and Further Reading
==============================

The following publications are relevant to the algorithms described
in this chapter,

* W.Hardle and W. Steiger, *Optimal Median Smoothing*, Appl. Statist., 44 (2), 1995.

* D. Lemire, *Streaming Maximum-Minimum Filter Using No More than Three Comparisons per Element*,
  Nordic Journal of Computing, 13 (4), 2006 (https://arxiv.org/abs/cs/0610046).

* B. P. Welford, *Note on a method for calculating corrected sums of squares and products*,
  Technometrics, 4 (3), 1962.
