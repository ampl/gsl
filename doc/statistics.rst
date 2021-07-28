.. index::
   single: statistics
   single: mean
   single: standard deviation
   single: variance
   single: estimated standard deviation
   single: estimated variance
   single: t-test
   single: range
   single: min
   single: max

**********
Statistics
**********

This chapter describes the statistical functions in the library.  The
basic statistical functions include routines to compute the mean,
variance and standard deviation.  More advanced functions allow you to
calculate absolute deviations, skewness, and kurtosis as well as the
median and arbitrary percentiles.  The algorithms use recurrence
relations to compute average quantities in a stable way, without large
intermediate values that might overflow. 

The functions are available in versions for datasets in the standard
floating-point and integer types.  The versions for double precision
floating-point data have the prefix :code:`gsl_stats` and are declared in
the header file :file:`gsl_statistics_double.h`.  The versions for integer
data have the prefix :code:`gsl_stats_int` and are declared in the header
file :file:`gsl_statistics_int.h`.   All the functions operate on C 
arrays with a :code:`stride` parameter specifying the spacing between 
elements.  

Mean, Standard Deviation and Variance
=====================================

.. function:: double gsl_stats_mean (const double data[], size_t stride, size_t n)

   This function returns the arithmetic mean of :data:`data`, a dataset of
   length :data:`n` with stride :data:`stride`.  The arithmetic mean, or
   *sample mean*, is denoted by :math:`\Hat\mu` and defined as,

   .. math:: \Hat\mu = {1 \over N} \sum x_i

   where :math:`x_i` are the elements of the dataset :data:`data`.  For
   samples drawn from a gaussian distribution the variance of
   :math:`\Hat\mu` is :math:`\sigma^2 / N`.

.. function:: double gsl_stats_variance (const double data[], size_t stride, size_t n)

   This function returns the estimated, or *sample*, variance of
   :data:`data`, a dataset of length :data:`n` with stride :data:`stride`.  The
   estimated variance is denoted by :math:`\Hat\sigma^2` and is defined by,

   .. only:: not texinfo

      .. math:: {\Hat\sigma}^2 = {1 \over (N-1)} \sum (x_i - {\Hat\mu})^2

   .. only:: texinfo

      ::

         \Hat\sigma^2 = (1/(N-1)) \sum (x_i - \Hat\mu)^2

   where :math:`x_i` are the elements of the dataset :data:`data`.  Note that
   the normalization factor of :math:`1/(N-1)` results from the derivation
   of :math:`\Hat\sigma^2` as an unbiased estimator of the population
   variance :math:`\sigma^2`.  For samples drawn from a Gaussian distribution
   the variance of :math:`\Hat\sigma^2` itself is :math:`2 \sigma^4 / N`.

   This function computes the mean via a call to :func:`gsl_stats_mean`.  If
   you have already computed the mean then you can pass it directly to
   :func:`gsl_stats_variance_m`.

.. function:: double gsl_stats_variance_m (const double data[], size_t stride, size_t n, double mean)

   This function returns the sample variance of :data:`data` relative to the
   given value of :data:`mean`.  The function is computed with :math:`\Hat\mu`
   replaced by the value of :data:`mean` that you supply,

   .. only:: not texinfo

      .. math:: {\Hat\sigma}^2 = {1 \over (N-1)} \sum (x_i - mean)^2

   .. only:: texinfo

      ::

         \Hat\sigma^2 = (1/(N-1)) \sum (x_i - mean)^2

.. function:: double gsl_stats_sd (const double data[], size_t stride, size_t n)
              double gsl_stats_sd_m (const double data[], size_t stride, size_t n, double mean)

   The standard deviation is defined as the square root of the variance.
   These functions return the square root of the corresponding variance
   functions above.

.. function:: double gsl_stats_tss (const double data[], size_t stride, size_t n)
              double gsl_stats_tss_m (const double data[], size_t stride, size_t n, double mean)

   These functions return the total sum of squares (TSS) of :data:`data` about
   the mean.  For :func:`gsl_stats_tss_m` the user-supplied value of
   :data:`mean` is used, and for :func:`gsl_stats_tss` it is computed using
   :func:`gsl_stats_mean`.

   .. only:: not texinfo

      .. math:: {\rm TSS} = \sum (x_i - mean)^2

   .. only:: texinfo

      ::

         TSS =  \sum (x_i - mean)^2

.. function:: double gsl_stats_variance_with_fixed_mean (const double data[], size_t stride, size_t n, double mean)

   This function computes an unbiased estimate of the variance of
   :data:`data` when the population mean :data:`mean` of the underlying
   distribution is known *a priori*.  In this case the estimator for
   the variance uses the factor :math:`1/N` and the sample mean
   :math:`\Hat\mu` is replaced by the known population mean :math:`\mu`,

   .. only:: not texinfo

      .. math:: {\Hat\sigma}^2 = {1 \over N} \sum (x_i - \mu)^2

   .. only:: texinfo

      ::

         \Hat\sigma^2 = (1/N) \sum (x_i - \mu)^2

.. function:: double gsl_stats_sd_with_fixed_mean (const double data[], size_t stride, size_t n, double mean)

   This function calculates the standard deviation of :data:`data` for a
   fixed population mean :data:`mean`.  The result is the square root of the
   corresponding variance function.

Absolute deviation
==================

.. function:: double gsl_stats_absdev (const double data[], size_t stride, size_t n)

   This function computes the absolute deviation from the mean of
   :data:`data`, a dataset of length :data:`n` with stride :data:`stride`.  The
   absolute deviation from the mean is defined as,

   .. only:: not texinfo

      .. math:: absdev  = {1 \over N} \sum |x_i - {\Hat\mu}|

   .. only:: texinfo

      ::

         absdev  = (1/N) \sum |x_i - \Hat\mu|

   where :math:`x_i` are the elements of the dataset :data:`data`.  The
   absolute deviation from the mean provides a more robust measure of the
   width of a distribution than the variance.  This function computes the
   mean of :data:`data` via a call to :func:`gsl_stats_mean`.

.. function:: double gsl_stats_absdev_m (const double data[], size_t stride, size_t n, double mean)

   This function computes the absolute deviation of the dataset :data:`data`
   relative to the given value of :data:`mean`,

   .. only:: not texinfo

      .. math:: absdev  = {1 \over N} \sum |x_i - mean|

   .. only:: texinfo

      ::

         absdev  = (1/N) \sum |x_i - mean|

   This function is useful if you have already computed the mean of
   :data:`data` (and want to avoid recomputing it), or wish to calculate the
   absolute deviation relative to another value (such as zero, or the
   median).

.. index:: skewness, kurtosis

Higher moments (skewness and kurtosis)
======================================

.. function:: double gsl_stats_skew (const double data[], size_t stride, size_t n)

   This function computes the skewness of :data:`data`, a dataset of length
   :data:`n` with stride :data:`stride`.  The skewness is defined as,

   .. only:: not texinfo

      .. math::

         skew = {1 \over N} \sum 
          {\left( x_i - {\Hat\mu} \over {\Hat\sigma} \right)}^3

   .. only:: texinfo

      ::

         skew = (1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^3

   where :math:`x_i` are the elements of the dataset :data:`data`.  The skewness
   measures the asymmetry of the tails of a distribution.

   The function computes the mean and estimated standard deviation of
   :data:`data` via calls to :func:`gsl_stats_mean` and :func:`gsl_stats_sd`.

.. function:: double gsl_stats_skew_m_sd (const double data[], size_t stride, size_t n, double mean, double sd)

   This function computes the skewness of the dataset :data:`data` using the
   given values of the mean :data:`mean` and standard deviation :data:`sd`,

   .. only:: not texinfo

      .. math:: skew = {1 \over N} \sum {\left( x_i - mean \over sd \right)}^3

   .. only:: texinfo

      ::

         skew = (1/N) \sum ((x_i - mean)/sd)^3

   These functions are useful if you have already computed the mean and
   standard deviation of :data:`data` and want to avoid recomputing them.

.. function:: double gsl_stats_kurtosis (const double data[], size_t stride, size_t n)

   This function computes the kurtosis of :data:`data`, a dataset of length
   :data:`n` with stride :data:`stride`.  The kurtosis is defined as,

   .. only:: not texinfo

      .. math::

         kurtosis = \left( {1 \over N} \sum 
          {\left(x_i - {\Hat\mu} \over {\Hat\sigma} \right)}^4 
          \right) 
          - 3

   .. only:: texinfo

      ::

         kurtosis = ((1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^4)  - 3

   The kurtosis measures how sharply peaked a distribution is, relative to
   its width.  The kurtosis is normalized to zero for a Gaussian
   distribution.

.. function:: double gsl_stats_kurtosis_m_sd (const double data[], size_t stride, size_t n, double mean, double sd)

   This function computes the kurtosis of the dataset :data:`data` using the
   given values of the mean :data:`mean` and standard deviation :data:`sd`,

   .. only:: not texinfo

      .. math::

         kurtosis = {1 \over N}
           \left( \sum {\left(x_i - mean \over sd \right)}^4 \right) 
           - 3

   .. only:: texinfo

      ::

         kurtosis = ((1/N) \sum ((x_i - mean)/sd)^4) - 3

   This function is useful if you have already computed the mean and
   standard deviation of :data:`data` and want to avoid recomputing them.

Autocorrelation
===============

.. function:: double gsl_stats_lag1_autocorrelation (const double data[], const size_t stride, const size_t n)

   This function computes the lag-1 autocorrelation of the dataset :data:`data`.

   .. only:: not texinfo

      .. math::

         a_1 = {\sum_{i = 2}^{n} (x_{i} - \Hat\mu) (x_{i-1} - \Hat\mu)
         \over
         \sum_{i = 1}^{n} (x_{i} - \Hat\mu) (x_{i} - \Hat\mu)}

   .. only:: texinfo

      ::

         a_1 = {\sum_{i = 2}^{n} (x_{i} - \Hat\mu) (x_{i-1} - \Hat\mu)
                \over
                \sum_{i = 1}^{n} (x_{i} - \Hat\mu) (x_{i} - \Hat\mu)}

.. function:: double gsl_stats_lag1_autocorrelation_m (const double data[], const size_t stride, const size_t n, const double mean)

   This function computes the lag-1 autocorrelation of the dataset
   :data:`data` using the given value of the mean :data:`mean`.

.. index::
   single: covariance, of two datasets

Covariance
==========

.. function:: double gsl_stats_covariance (const double data1[], const size_t stride1, const double data2[], const size_t stride2, const size_t n)

   This function computes the covariance of the datasets :data:`data1` and
   :data:`data2` which must both be of the same length :data:`n`.

   .. only:: not texinfo

      .. math:: covar = {1 \over (n - 1)} \sum_{i = 1}^{n} (x_{i} - \Hat x) (y_{i} - \Hat y)

   .. only:: texinfo

      ::

         covar = (1/(n - 1)) \sum_{i = 1}^{n} (x_i - \Hat x) (y_i - \Hat y)

.. function:: double gsl_stats_covariance_m (const double data1[], const size_t stride1, const double data2[], const size_t stride2, const size_t n, const double mean1, const double mean2)

   This function computes the covariance of the datasets :data:`data1` and
   :data:`data2` using the given values of the means, :data:`mean1` and
   :data:`mean2`.  This is useful if you have already computed the means of
   :data:`data1` and :data:`data2` and want to avoid recomputing them.

.. index::
   single: correlation, of two datasets

Correlation
===========

.. function:: double gsl_stats_correlation (const double data1[], const size_t stride1, const double data2[], const size_t stride2, const size_t n)

   This function efficiently computes the Pearson correlation coefficient
   between the datasets :data:`data1` and :data:`data2` which must both be of
   the same length :data:`n`.

   .. only:: not texinfo

      .. math::

         r = {cov(x, y) \over \Hat\sigma_x \Hat\sigma_y} =
         {{1 \over n-1} \sum (x_i - \Hat x) (y_i - \Hat y)
         \over
         \sqrt{{1 \over n-1} \sum (x_i - {\Hat x})^2}
         \sqrt{{1 \over n-1} \sum (y_i - {\Hat y})^2}
         }

   .. only:: texinfo

      ::

         r = cov(x, y) / (\Hat\sigma_x \Hat\sigma_y)
           = {1/(n-1) \sum (x_i - \Hat x) (y_i - \Hat y)
              \over
              \sqrt{1/(n-1) \sum (x_i - \Hat x)^2} \sqrt{1/(n-1) \sum (y_i - \Hat y)^2}
             }

.. function:: double gsl_stats_spearman (const double data1[], const size_t stride1, const double data2[], const size_t stride2, const size_t n, double work[])

   This function computes the Spearman rank correlation coefficient between
   the datasets :data:`data1` and :data:`data2` which must both be of the same
   length :data:`n`. Additional workspace of size 2 * :data:`n` is required in
   :data:`work`. The Spearman rank correlation between vectors :math:`x` and
   :math:`y` is equivalent to the Pearson correlation between the ranked
   vectors :math:`x_R` and :math:`y_R`, where ranks are defined to be the
   average of the positions of an element in the ascending order of the values.

Weighted Samples
================

The functions described in this section allow the computation of
statistics for weighted samples.  The functions accept an array of
samples, :math:`x_i`, with associated weights, :math:`w_i`.  Each sample
:math:`x_i` is considered as having been drawn from a Gaussian
distribution with variance :math:`\sigma_i^2`.  The sample weight
:math:`w_i` is defined as the reciprocal of this variance, :math:`w_i = 1/\sigma_i^2`.
Setting a weight to zero corresponds to removing a sample from a dataset.

.. function:: double gsl_stats_wmean (const double w[], size_t wstride, const double data[], size_t stride, size_t n)

   This function returns the weighted mean of the dataset :data:`data` with
   stride :data:`stride` and length :data:`n`, using the set of weights :data:`w`
   with stride :data:`wstride` and length :data:`n`.  The weighted mean is defined as,

   .. only:: not texinfo

      .. math:: {\Hat\mu} = {{\sum w_i x_i} \over {\sum w_i}}

   .. only:: texinfo

      ::

         \Hat\mu = (\sum w_i x_i) / (\sum w_i)

.. function:: double gsl_stats_wvariance (const double w[], size_t wstride, const double data[], size_t stride, size_t n)

   This function returns the estimated variance of the dataset :data:`data`
   with stride :data:`stride` and length :data:`n`, using the set of weights
   :data:`w` with stride :data:`wstride` and length :data:`n`.  The estimated
   variance of a weighted dataset is calculated as,

   .. only:: not texinfo

      .. math::

         \Hat\sigma^2 = {{\sum w_i} \over {(\sum w_i)^2 - \sum (w_i^2)}} 
                         \sum w_i (x_i - \Hat\mu)^2

   .. only:: texinfo

      ::

         \Hat\sigma^2 = ((\sum w_i)/((\sum w_i)^2 - \sum (w_i^2))) 
                         \sum w_i (x_i - \Hat\mu)^2

   Note that this expression reduces to an unweighted variance with the
   familiar :math:`1/(N-1)` factor when there are :math:`N` equal non-zero
   weights.

.. function:: double gsl_stats_wvariance_m (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean)

   This function returns the estimated variance of the weighted dataset
   :data:`data` using the given weighted mean :data:`wmean`.

.. function:: double gsl_stats_wsd (const double w[], size_t wstride, const double data[], size_t stride, size_t n)

   The standard deviation is defined as the square root of the variance.
   This function returns the square root of the corresponding variance
   function :func:`gsl_stats_wvariance` above.

.. function:: double gsl_stats_wsd_m (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean)

   This function returns the square root of the corresponding variance
   function :func:`gsl_stats_wvariance_m` above.

.. function:: double gsl_stats_wvariance_with_fixed_mean (const double w[], size_t wstride, const double data[], size_t stride, size_t n, const double mean)

   This function computes an unbiased estimate of the variance of the weighted
   dataset :data:`data` when the population mean :data:`mean` of the underlying
   distribution is known *a priori*.  In this case the estimator for
   the variance replaces the sample mean :math:`\Hat\mu` by the known
   population mean :math:`\mu`,

   .. only:: not texinfo

      .. math:: \Hat\sigma^2 = {{\sum w_i (x_i - \mu)^2} \over {\sum w_i}}

   .. only:: texinfo

      ::

         \Hat\sigma^2 = (\sum w_i (x_i - \mu)^2) / (\sum w_i)

.. function:: double gsl_stats_wsd_with_fixed_mean (const double w[], size_t wstride, const double data[], size_t stride, size_t n, const double mean)

   The standard deviation is defined as the square root of the variance.
   This function returns the square root of the corresponding variance
   function above.

.. function:: double gsl_stats_wtss (const double w[], const size_t wstride, const double data[], size_t stride, size_t n)
              double gsl_stats_wtss_m (const double w[], const size_t wstride, const double data[], size_t stride, size_t n, double wmean)

   These functions return the weighted total sum of squares (TSS) of
   :data:`data` about the weighted mean.  For :func:`gsl_stats_wtss_m` the
   user-supplied value of :data:`wmean` is used, and for :func:`gsl_stats_wtss`
   it is computed using :func:`gsl_stats_wmean`.

   .. only:: not texinfo

      .. math:: {\rm TSS} = \sum w_i (x_i - wmean)^2

   .. only:: texinfo

      ::

         TSS =  \sum w_i (x_i - wmean)^2

.. function:: double gsl_stats_wabsdev (const double w[], size_t wstride, const double data[], size_t stride, size_t n)

   This function computes the weighted absolute deviation from the weighted
   mean of :data:`data`.  The absolute deviation from the mean is defined as,

   .. only:: not texinfo

      .. math:: absdev = {{\sum w_i |x_i - \Hat\mu|} \over {\sum w_i}}

   .. only:: texinfo

      ::

         absdev = (\sum w_i |x_i - \Hat\mu|) / (\sum w_i)

.. function:: double gsl_stats_wabsdev_m (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean)

   This function computes the absolute deviation of the weighted dataset
   :data:`data` about the given weighted mean :data:`wmean`.

.. function:: double gsl_stats_wskew (const double w[], size_t wstride, const double data[], size_t stride, size_t n)

   This function computes the weighted skewness of the dataset :data:`data`.

   .. only:: not texinfo

      .. math:: skew = {{\sum w_i ((x_i - {\Hat x})/{\Hat \sigma})^3} \over {\sum w_i}}

   .. only:: texinfo

      ::

         skew = (\sum w_i ((x_i - \Hat x)/\Hat \sigma)^3) / (\sum w_i)

.. function:: double gsl_stats_wskew_m_sd (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean, double wsd)

   This function computes the weighted skewness of the dataset :data:`data`
   using the given values of the weighted mean and weighted standard
   deviation, :data:`wmean` and :data:`wsd`.

.. function:: double gsl_stats_wkurtosis (const double w[], size_t wstride, const double data[], size_t stride, size_t n)

   This function computes the weighted kurtosis of the dataset :data:`data`.

   .. only:: not texinfo

      .. math:: kurtosis = {{\sum w_i ((x_i - {\Hat x})/{\Hat \sigma})^4} \over {\sum w_i}} - 3

   .. only:: texinfo

      ::

         kurtosis = ((\sum w_i ((x_i - \Hat x)/\Hat \sigma)^4) / (\sum w_i)) - 3

.. function:: double gsl_stats_wkurtosis_m_sd (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean, double wsd)

   This function computes the weighted kurtosis of the dataset :data:`data`
   using the given values of the weighted mean and weighted standard
   deviation, :data:`wmean` and :data:`wsd`.

Maximum and Minimum values
==========================

The following functions find the maximum and minimum values of a
dataset (or their indices).  If the data contains :code:`NaN`-s then a
:code:`NaN` will be returned, since the maximum or minimum value is
undefined.  For functions which return an index, the location of the
first :code:`NaN` in the array is returned.

.. function:: double gsl_stats_max (const double data[], size_t stride, size_t n)

   This function returns the maximum value in :data:`data`, a dataset of
   length :data:`n` with stride :data:`stride`.  The maximum value is defined
   as the value of the element :math:`x_i` which satisfies :math:`x_i \ge x_j`
   for all :math:`j`.

   If you want instead to find the element with the largest absolute
   magnitude you will need to apply :func:`fabs` or :func:`abs` to your data
   before calling this function.

.. function:: double gsl_stats_min (const double data[], size_t stride, size_t n)

   This function returns the minimum value in :data:`data`, a dataset of
   length :data:`n` with stride :data:`stride`.  The minimum value is defined
   as the value of the element :math:`x_i` which satisfies :math:`x_i \le x_j`
   for all :math:`j`.

   If you want instead to find the element with the smallest absolute
   magnitude you will need to apply :func:`fabs` or :func:`abs` to your data
   before calling this function.

.. function:: void gsl_stats_minmax (double * min, double * max, const double data[], size_t stride, size_t n)

   This function finds both the minimum and maximum values :data:`min`,
   :data:`max` in :data:`data` in a single pass.

.. function:: size_t gsl_stats_max_index (const double data[], size_t stride, size_t n)

   This function returns the index of the maximum value in :data:`data`, a
   dataset of length :data:`n` with stride :data:`stride`.  The maximum value is
   defined as the value of the element :math:`x_i` which satisfies 
   :math:`x_i \ge x_j`
   for all :math:`j`.  When there are several equal maximum
   elements then the first one is chosen.

.. function:: size_t gsl_stats_min_index (const double data[], size_t stride, size_t n)

   This function returns the index of the minimum value in :data:`data`, a
   dataset of length :data:`n` with stride :data:`stride`.  The minimum value
   is defined as the value of the element :math:`x_i` which satisfies
   :math:`x_i \ge x_j`
   for all :math:`j`.  When there are several equal
   minimum elements then the first one is chosen.

.. function:: void gsl_stats_minmax_index (size_t * min_index, size_t * max_index, const double data[], size_t stride, size_t n)

   This function returns the indexes :data:`min_index`, :data:`max_index` of
   the minimum and maximum values in :data:`data` in a single pass.

Median and Percentiles
======================

The median and percentile functions described in this section operate on
sorted data in :math:`O(1)` time. There is also a routine for computing
the median of an unsorted input array in average :math:`O(n)` time using
the quickselect algorithm. For convenience we use *quantiles*, measured on a scale
of 0 to 1, instead of percentiles (which use a scale of 0 to 100).

.. function:: double gsl_stats_median_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n)

   This function returns the median value of :data:`sorted_data`, a dataset
   of length :data:`n` with stride :data:`stride`.  The elements of the array
   must be in ascending numerical order.  There are no checks to see
   whether the data are sorted, so the function :func:`gsl_sort` should
   always be used first.

   When the dataset has an odd number of elements the median is the value
   of element :math:`(n-1)/2`.  When the dataset has an even number of
   elements the median is the mean of the two nearest middle values,
   elements :math:`(n-1)/2` and :math:`n/2`.  Since the algorithm for
   computing the median involves interpolation this function always returns
   a floating-point number, even for integer data types.

.. function:: double gsl_stats_median (double data[], const size_t stride, const size_t n)

   This function returns the median value of :data:`data`, a dataset
   of length :data:`n` with stride :data:`stride`. The median is found
   using the quickselect algorithm. The input array does not need to be
   sorted, but note that the algorithm rearranges the array and so the input
   is not preserved on output.

.. function:: double gsl_stats_quantile_from_sorted_data (const double sorted_data[], size_t stride, size_t n, double f)

   This function returns a quantile value of :data:`sorted_data`, a
   double-precision array of length :data:`n` with stride :data:`stride`.  The
   elements of the array must be in ascending numerical order.  The
   quantile is determined by the :data:`f`, a fraction between 0 and 1.  For
   example, to compute the value of the 75th percentile :data:`f` should have
   the value 0.75.

   There are no checks to see whether the data are sorted, so the function
   :func:`gsl_sort` should always be used first.

   The quantile is found by interpolation, using the formula

   .. only:: not texinfo

      .. math:: \hbox{quantile} = (1 - \delta) x_i + \delta x_{i+1}

   .. only:: texinfo

      ::

         quantile = (1 - \delta) x_i + \delta x_{i+1}

   where :math:`i` is :code:`floor((n - 1)f)` and :math:`\delta` is
   :math:`(n-1)f - i`.

   Thus the minimum value of the array (:code:`data[0*stride]`) is given by
   :data:`f` equal to zero, the maximum value (:code:`data[(n-1)*stride]`) is
   given by :data:`f` equal to one and the median value is given by :data:`f`
   equal to 0.5.  Since the algorithm for computing quantiles involves
   interpolation this function always returns a floating-point number, even
   for integer data types.

.. @node Statistical tests
.. @section Statistical tests

.. FIXME, do more work on the statistical tests

.. -@deftypefun double gsl_stats_ttest (const double data1[], double data2[], size_t n1, size_t n2)
.. -@deftypefunx Statistics double gsl_stats_int_ttest (const double data1[], double data2[], size_t n1, size_t n2)

.. The function :func:`gsl_stats_ttest` computes the t-test statistic for
.. the two arrays :data:`data1`[] and :data:`data2`[], of lengths :data:`n1` and
.. -:data:`n2` respectively.

.. The t-test statistic measures the difference between the means of two
.. datasets.

Order Statistics
================

The :math:`k`-th *order statistic* of a sample is equal to its :math:`k`-th smallest value.
The :math:`k`-th order statistic of a set of :math:`n` values :math:`x = \left\{ x_i \right\}, 1 \le i \le n` is
denoted :math:`x_{(k)}`. The median of the set :math:`x` is equal to :math:`x_{\left( \frac{n}{2} \right)}` if
:math:`n` is odd, or the average of :math:`x_{\left( \frac{n}{2} \right)}` and :math:`x_{\left( \frac{n}{2} + 1 \right)}`
if :math:`n` is even. The :math:`k`-th smallest element of a length :math:`n` vector can be found
in average :math:`O(n)` time using the quickselect algorithm.

.. function:: double gsl_stats_select(double data[], const size_t stride, const size_t n, const size_t k)

   This function finds the :data:`k`-th smallest element of the input array :data:`data`
   of length :data:`n` and stride :data:`stride` using the quickselect method. The
   algorithm rearranges the elements of :data:`data` and so the input array is not preserved
   on output.

.. index::
   single: robust location estimators
   single: location estimation
   single: estimation, location

Robust Location Estimates
=========================

A *location estimate* refers to a typical or central value which best describes a given
dataset. The mean and median are both examples of location estimators. However, the
mean has a severe sensitivity to data outliers and can give erroneous values when
even a small number of outliers are present. The median on the other hand, has
a strong insensitivity to data outliers, but due to its non-smoothness it can
behave unexpectedly in certain situations. GSL offers the following alternative
location estimators, which are robust to the presence of outliers.

.. index::
   single: trimmed mean
   single: truncated mean
   single: mean, trimmed
   single: mean, truncated

Trimmed Mean
------------

The trimmed mean, or *truncated mean*, discards a certain number of smallest and largest
samples from the input vector before computing the mean of the remaining samples. The
amount of trimming is specified by a factor :math:`\alpha \in [0,0.5]`. Then the
number of samples discarded from both ends of the input vector is
:math:`\left\lfloor \alpha n \right\rfloor`, where :math:`n` is the length of the input.
So to discard 25% of the samples from each end, one would set :math:`\alpha = 0.25`.

.. function:: double gsl_stats_trmean_from_sorted_data (const double alpha, const double sorted_data[], const size_t stride, const size_t n)

   This function returns the trimmed mean of :data:`sorted_data`, a dataset
   of length :data:`n` with stride :data:`stride`. The elements of the array
   must be in ascending numerical order.  There are no checks to see
   whether the data are sorted, so the function :func:`gsl_sort` should
   always be used first. The trimming factor :math:`\alpha` is given in :data:`alpha`.
   If :math:`\alpha \ge 0.5`, then the median of the input is returned.

.. index::
   single: Gastwirth estimator

Gastwirth Estimator
-------------------

Gastwirth's location estimator is a weighted sum of three order statistics,

.. only:: not texinfo

   .. math:: gastwirth = 0.3 \times Q_{\frac{1}{3}} + 0.4 \times Q_{\frac{1}{2}} + 0.3 \times Q_{\frac{2}{3}}

.. only:: texinfo

   ::

      gastwirth = 0.3 * Q_{1/3} + 0.4 * Q_{1/2} + 0.3 * Q_{2/3}

where :math:`Q_{\frac{1}{3}}` is the one-third quantile, :math:`Q_{\frac{1}{2}}` is the one-half
quantile (i.e. median), and :math:`Q_{\frac{2}{3}}` is the two-thirds quantile.

.. function:: double gsl_stats_gastwirth_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n)

   This function returns the Gastwirth location estimator of :data:`sorted_data`, a dataset
   of length :data:`n` with stride :data:`stride`.  The elements of the array
   must be in ascending numerical order.  There are no checks to see
   whether the data are sorted, so the function :func:`gsl_sort` should
   always be used first.

.. index::
   single: robust scale estimators
   single: scale estimation
   single: estimation, scale

Robust Scale Estimates
======================

A *robust scale estimate*, also known as a robust measure of scale, attempts to quantify
the statistical dispersion (variability, scatter, spread) in a set of data which may contain outliers.
In such datasets, the usual variance or standard deviation scale estimate can be rendered useless
by even a single outlier.

.. index::
   single: median absolute deviation

.. _sec_mad-statistic:

Median Absolute Deviation (MAD)
-------------------------------

The median absolute deviation (MAD) is defined as

.. only:: not texinfo

   .. math:: MAD = 1.4826 \times \textrm{median} \left\{ \left| x_i - \textrm{median} \left( x \right) \right| \right\}

.. only:: texinfo

   ::

      MAD = 1.4826 median { | x_i - median(x) | }

In words, first the median of all samples is computed. Then the median
is subtracted from all samples in the input to find the deviation of each sample
from the median. The median of all absolute deviations is then the MAD.
The factor :math:`1.4826` makes the MAD an unbiased estimator of the standard deviation for Gaussian data.
The median absolute deviation has an asymptotic efficiency of 37%.

.. function:: double gsl_stats_mad0 (const double data[], const size_t stride, const size_t n, double work[])
.. function:: double gsl_stats_mad (const double data[], const size_t stride, const size_t n, double work[])

   These functions return the median absolute deviation of :data:`data`, a dataset
   of length :data:`n` and stride :data:`stride`.
   The :code:`mad0` function calculates
   :math:`\textrm{median} \left\{ \left| x_i - \textrm{median} \left( x \right) \right| \right\}`
   (i.e. the :math:`MAD` statistic without the bias correction scale factor).
   These functions require additional workspace of size :code:`n` provided in :data:`work`.

.. index::
   single: Sn statistic

.. _sec_Sn-statistic:

:math:`S_n` Statistic
---------------------

The :math:`S_n` statistic developed by Croux and Rousseeuw is defined as

.. only:: not texinfo

   .. math:: S_n = 1.1926 \times c_n \times \textrm{median}_i \left\{ \textrm{median}_j \left( \left| x_i - x_j \right| \right) \right\}

.. only:: texinfo

   ::

      S_n = 1.1926 * c_n * median_i { median_j ( | x_i - x_j | ) }

For each sample :math:`x_i, 1 \le i \le n`, the median of the values :math:`\left| x_i - x_j \right|` is computed for all
:math:`x_j, 1 \le j \le n`. This yields :math:`n` values, whose median then gives the final :math:`S_n`.
The factor :math:`1.1926` makes :math:`S_n` an unbiased estimate of the standard deviation for Gaussian data.
The factor :math:`c_n` is a correction factor to correct bias in small sample sizes. :math:`S_n` has an asymptotic
efficiency of 58%.

.. function:: double gsl_stats_Sn0_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n, double work[])
.. function:: double gsl_stats_Sn_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n, double work[])

   These functions return the :math:`S_n` statistic of :data:`sorted_data`, a dataset
   of length :data:`n` with stride :data:`stride`.  The elements of the array
   must be in ascending numerical order.  There are no checks to see
   whether the data are sorted, so the function :func:`gsl_sort` should
   always be used first. The :code:`Sn0` function calculates
   :math:`\textrm{median}_i \left\{ \textrm{median}_j \left( \left| x_i - x_j \right| \right) \right\}`
   (i.e. the :math:`S_n` statistic without the bias correction scale factors).
   These functions require additional workspace of size
   :code:`n` provided in :data:`work`.

.. index::
   single: Qn statistic

.. _sec_Qn-statistic:

:math:`Q_n` Statistic
---------------------

The :math:`Q_n` statistic developed by Croux and Rousseeuw is defined as

.. only:: not texinfo

   .. math:: Q_n = 2.21914 \times d_n \times \left\{ \left| x_i - x_j \right|, i < j \right\}_{(k)}

.. only:: texinfo

   ::

      Q_n = 2.21914 * d_n * { | x_i - x_j |, i < j }_{(k)}

The factor :math:`2.21914` makes :math:`Q_n` an unbiased estimate of the standard deviation for Gaussian data.
The factor :math:`d_n` is a correction factor to correct bias in small sample sizes. The order statistic
is

.. only:: not texinfo

   .. math:: k = \left(
                   \begin{array}{c}
                     \left\lfloor \frac{n}{2} \right\rfloor + 1 \\
                     2
                   \end{array}
                 \right)

.. only:: texinfo

   ::

      k = ( floor(n/2) + 1 )
          (       2        )

:math:`Q_n` has an asymptotic efficiency of 82%.

.. function:: double gsl_stats_Qn0_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n, double work[], int work_int[])
              double gsl_stats_Qn_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n, double work[], int work_int[])

   These functions return the :math:`Q_n` statistic of :data:`sorted_data`, a dataset
   of length :data:`n` with stride :data:`stride`. The elements of the array
   must be in ascending numerical order.  There are no checks to see
   whether the data are sorted, so the function :func:`gsl_sort` should
   always be used first. The :code:`Qn0` function calculates
   :math:`\left\{ \left| x_i - x_j \right|, i < j \right\}_{(k)}`
   (i.e. :math:`Q_n` without the bias correction scale factors).
   These functions require additional workspace of size
   :code:`3n` provided in :data:`work` and integer workspace of size :code:`5n`
   provided in :data:`work_int`.

Examples
========

Here is a basic example of how to use the statistical functions:

.. include:: examples/stat.c
   :code:

The program should produce the following output,

.. include:: examples/stat.txt
   :code:

Here is an example using sorted data,

.. include:: examples/statsort.c
   :code:

This program should produce the following output,

.. include:: examples/statsort.txt
   :code:

References and Further Reading
==============================

The standard reference for almost any topic in statistics is the
multi-volume *Advanced Theory of Statistics* by Kendall and Stuart.

* Maurice Kendall, Alan Stuart, and J. Keith Ord.
  *The Advanced Theory of Statistics* (multiple volumes)
  reprinted as *Kendall's Advanced Theory of Statistics*.
  Wiley, ISBN 047023380X.

Many statistical concepts can be more easily understood by a Bayesian
approach.  The following book by Gelman, Carlin, Stern and Rubin gives a
comprehensive coverage of the subject.

* Andrew Gelman, John B. Carlin, Hal S. Stern, Donald B. Rubin.
  *Bayesian Data Analysis*.
  Chapman & Hall, ISBN 0412039915.

For physicists the Particle Data Group provides useful reviews of
Probability and Statistics in the "Mathematical Tools" section of its
Annual Review of Particle Physics. 

* *Review of Particle Properties*,
  R.M. Barnett et al., Physical Review D54, 1 (1996)

The Review of Particle Physics is available online at
the website http://pdg.lbl.gov/.

The following papers describe robust scale estimation,

* C. Croux and P. J. Rousseeuw, *Time-Efficient algorithms for two highly robust
  estimators of scale*, Comp. Stat., Physica, Heidelberg, 1992.
* P. J. Rousseeuw and C. Croux, *Explicit scale estimators with high breakdown point*,
  L1-Statistical Analysis and Related Methods, pp. 77-92, 1992.
