*****************
Digital Filtering
*****************

Introduction
============

The filters discussed in this chapter are based on the following moving data
window which is centered on :math:`i`-th sample:

.. only:: not texinfo

   .. math:: W_i^H = \left\{ x_{i-H}, \dots, x_i, \dots, x_{i+H} \right\}

.. only:: texinfo

   ::

      W_i^H = { x_{i-H}, ..., x_i, ..., x_{i+H} }

Here, :math:`H` is a non-negative integer called the *window half-length*, which
represents the number of samples before and after sample :math:`i`.
The total window length is :math:`K = 2 H + 1`.

Handling Endpoints
==================

When processing samples near the ends of the input signal, there will not
be enough samples to fill the window :math:`W_i^H` defined above.
Therefore the user must specify how to construct the windows near the end points.
This is done by passing an input argument of type :type:`gsl_filter_end_t`:

.. type:: gsl_filter_end_t

   This data type specifies how to construct windows near end points and can
   be selected from the following choices:

   .. macro:: GSL_FILTER_END_PADZERO

      With this option, a full window of length :math:`K` will be constructed
      by inserting zeros into the window near the signal end points. Effectively,
      the input signal is modified to

      .. only:: not texinfo

         .. math:: \tilde{x} = \{ \underbrace{0, \dots, 0}_{H \textrm{ zeros}}, x_1, x_2, \dots, x_{n-1}, x_n, \underbrace{0, \dots, 0}_{H \textrm{ zeros} } \}

      .. only:: texinfo

         ::

            x~ = { 0, ..., 0, x_1, x_2, ..., x_{n-1}, x_n, 0, ..., 0 }

      to ensure a well-defined window for all :math:`x_i`.

   .. macro:: GSL_FILTER_END_PADVALUE

      With this option, a full window of length :math:`K` will be constructed
      by padding the window with the first and last sample in the input signal.
      Effectively, the input signal is modified to

      .. only:: not texinfo

         .. math:: \tilde{x} = \{ \underbrace{x_1, \dots, x_1}_{H}, x_1, x_2, \dots, x_{n-1}, x_n, \underbrace{x_n, \dots, x_n}_{H} \}

      .. only:: texinfo

         ::

            x~ = { x_1, ..., x_1, x_1, x_2, ..., x_{n-1}, x_n, x_n, ..., x_n }

   .. macro:: GSL_FILTER_END_TRUNCATE

      With this option, no padding is performed, and the windows are simply truncated
      as the end points are approached.

Linear Digital Filters
======================

Gaussian Filter
---------------

The Gaussian filter convolves the input signal with a Gaussian kernel or window. This filter
is often used as a smoothing or noise reduction filter. The Gaussian kernel is
defined by

.. only:: not texinfo

   .. math:: G(k) = e^{-\frac{1}{2} \left( \alpha \frac{k}{(K-1)/2} \right)^2} = e^{-k^2/2\sigma^2}

.. only:: texinfo

   ::

      G(k) = e^{-1/2 ( \alpha k/((K-1)/2) )^2} = e^{-k^2/2\sigma^2}

for :math:`-(K-1)/2 \le k \le (K-1)/2`, and :math:`K` is the size of the kernel. The
parameter :math:`\alpha` specifies the number of standard deviations :math:`\sigma` desired
in the kernel. So for example setting :math:`\alpha = 3` would define a Gaussian window
of length :math:`K` which spans :math:`\pm 3 \sigma`. It is often more convenient to specify
the parameter :math:`\alpha` rather than the standard deviation :math:`\sigma` when constructing
the kernel, since a fixed value of :math:`\alpha` would correspond to the same shape of
Gaussian regardless of the size :math:`K`. The appropriate value of the standard deviation
depends on :math:`K` and is related to :math:`\alpha` as

.. only:: not texinfo

   .. math:: \sigma = \frac{K - 1}{2\alpha}

.. only:: texinfo

   ::

      \sigma = (K - 1)/(2 \alpha)

The routines below accept :math:`\alpha` as an input argument instead of :math:`\sigma`.

The Gaussian filter offers a convenient way of differentiating and smoothing an input signal
in a single pass. Using the derivative property of a convolution,

.. only:: not texinfo

   .. math:: \frac{d}{dt} \left( G * x \right) = \frac{dG}{dt} * x

.. only:: texinfo

   ::

      d/dt ( G * x ) = dG/dt * x

the input signal :math:`x(t)` can be smoothed and differentiated at the same time by
convolution with a derivative Gaussian kernel, which can be readily computed from the
analytic expression above. The same principle applies to higher order derivatives.

.. function:: gsl_filter_gaussian_workspace * gsl_filter_gaussian_alloc(const size_t K)

   This function initializes a workspace for Gaussian filtering using a kernel of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(K)`.

.. function:: void gsl_filter_gaussian_free(gsl_filter_gaussian_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_gaussian(const gsl_filter_end_t endtype, const double alpha, const size_t order, const gsl_vector * x, gsl_vector * y, gsl_filter_gaussian_workspace * w)

   This function applies a Gaussian filter parameterized by :data:`alpha` to the input vector :data:`x`,
   storing the output in :data:`y`. The derivative order is specified by :data:`order`, with
   :code:`0` corresponding to a Gaussian, :code:`1` corresponding to a first derivative
   Gaussian, and so on. The parameter :data:`endtype` specifies how the signal end points are handled.
   It is allowed for :data:`x` = :data:`y` for an in-place filter.

.. function:: int gsl_filter_gaussian_kernel(const double alpha, const size_t order, const int normalize, gsl_vector * kernel)

   This function constructs a Gaussian kernel parameterized by :data:`alpha` and
   stores the output in :data:`kernel`. The parameter :data:`order` specifies the
   derivative order, with :code:`0` corresponding to a Gaussian, :code:`1` corresponding
   to a first derivative Gaussian, and so on. If :data:`normalize` is set to :code:`1`, then
   the kernel will be normalized to sum to one on output. If :data:`normalize` is set to
   :code:`0`, no normalization is performed.

Nonlinear Digital Filters
=========================

The nonlinear digital filters described below are based on the window median, which is given
by

.. only:: not texinfo

   .. math:: m_i = \textrm{median} \left\{ W_i^H \right\} = \textrm{median} \left\{ x_{i-H}, \dots, x_i, \dots, x_{i+H} \right\}

.. only:: texinfo

   ::

      m_i = median { W_i^H } = median { x_{i-H}, ..., x_i, ..., x_{i+H} }

The median is considered robust to local outliers, unlike the mean.
Median filters can preserve sharp edges while at the same removing signal noise, and are used
in a wide range of applications.

Standard Median Filter
----------------------

The *standard median filter* (SMF) simply replaces the sample :math:`x_i` by the median
:math:`m_i` of the window :math:`W_i^H`: This filter has one tuning parameter given
by :math:`H`. The standard median filter is considered highly resistant to
local outliers and local noise in the data sequence :math:`\{x_i\}`.

.. function:: gsl_filter_median_workspace * gsl_filter_median_alloc(const size_t K)

   This function initializes a workspace for standard median filtering using a symmetric centered moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(7K)`.

.. function:: void gsl_filter_median_free(gsl_filter_median_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_median(const gsl_filter_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_filter_median_workspace * w)

   This function applies a standard median filter to the input :data:`x`, storing the output in :data:`y`.
   The parameter :data:`endtype` specifies how the signal end points are handled. It
   is allowed to have :data:`x` = :data:`y` for an in-place filter.

Recursive Median Filter
-----------------------

The *recursive median filter* (RMF) is a modification of the SMF to include previous filter outputs
in the window before computing the median. The filter's response is

.. only:: not texinfo

   .. math:: y_i = \textrm{median} \left( y_{i-H}, \dots, y_{i-1}, x_i, x_{i+1}, \dots, x_{i+H} \right)

.. only:: texinfo

   ::

      y_i = median ( y_{i-H}, ..., y_{i-1}, x_i, x_{i+1}, ..., x_{i+H} )

Sometimes, the SMF must be applied several times in a row to achieve adequate smoothing (i.e. a cascade filter).
The RMF, on the other hand, converges to a *root sequence* in one pass,
and can sometimes provide a smoother result than several passes of the SMF. A root sequence is an input which is
left unchanged by the filter.  So there is no need to apply a recursive median filter twice to an input vector.

.. function:: gsl_filter_rmedian_workspace * gsl_filter_rmedian_alloc(const size_t K)

   This function initializes a workspace for recursive median filtering using a symmetric centered moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(K)`.

.. function:: void gsl_filter_rmedian_free(gsl_filter_rmedian_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_rmedian(const gsl_filter_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_filter_rmedian_workspace * w)

   This function applies a recursive median filter to the input :data:`x`, storing the output in :data:`y`.
   The parameter :data:`endtype` specifies how the signal end points are handled. It
   is allowed to have :data:`x` = :data:`y` for an in-place filter.

Impulse Detection Filter
------------------------

Impulsive noise is characterized by short sequences of data points distinct from those in the
surrounding neighborhood. This section describes a powerful class of filters, also known as
*impulse rejection filters* and *decision-based filters*, designed to detect and remove such outliers from data.
The filter's response is given by

.. only:: not texinfo

   .. math:: y_i = \left\{
                     \begin{array}{cc}
                       x_i, & |x_i - m_i| \le t S_i \\
                       m_i, & |x_i - m_i| > t S_i
                     \end{array}
                   \right.

.. only:: texinfo

   ::

      y_i = { x_i, |x_i - m_i| <= t * S_i
            { m_i, |x_i - m_i| > t * S_i

where :math:`m_i` is the median value of the window :math:`W_i^H`, :math:`S_i` is a robust estimate
of the scatter or dispersion for the window :math:`W_i^H`, and :math:`t` is a tuning parameter specifying
the number of scale factors needed to determine that a point is an outlier. The main idea is that the median
:math:`m_i` will be unaffected by a small number of outliers in the window, and so a given
sample :math:`x_i` is tested to determine how far away it is from the median in terms of the local
scale estimate :math:`S_i`. Samples which are more than :math:`t` scale estimates away from the median
are labeled as outliers and replaced by the window median :math:`m_i`. Samples which are less than
:math:`t` scale estimates from the median are left unchanged by the filter.

Note that when :math:`t = 0`, the impulse detection filter is equivalent to the standard median filter. When
:math:`t \rightarrow \infty`, it becomes the identity filter. This means the impulse detection filter can
be viewed as a "less aggressive" version of the standard median filter, becoming less aggressive as :math:`t` is
increased. Note that this filter modifies only samples identified as outliers, while the standard median
filter changes all samples to the local median, regardless of whether they are outliers. This fact, plus
the additional flexibility offered by the additional tuning parameter :math:`t` can make the impulse detection filter
a better choice for some applications.

It is important to have a robust and accurate scale estimate :math:`S_i` in order to
detect impulse outliers even in the presence of noise. The window standard deviation is not
typically a good choice, as it can be significantly perturbed by the presence of even one outlier.
GSL offers the following choices (specified by a parameter of type :type:`gsl_filter_scale_t`) for
computing the scale estimate :math:`S_i`, all of which are robust to the presence of impulse outliers.

.. type:: gsl_filter_scale_t

   This type specifies how the scale estimate :math:`S_i` of the window :math:`W_i^H` is calculated.

   .. macro:: GSL_FILTER_SCALE_MAD

      This option specifies the median absolute deviation (MAD) scale estimate, defined by

      .. only:: not texinfo

         .. math:: S_i = 1.4826 \times \textrm{median} \left\{ | W_i^H - m_i | \right\}

      .. only:: texinfo

         ::

            S_i = 1.4826 median { | W_i^H - m_i | }

      This choice of scale estimate is also known as the *Hampel filter* in the statistical literature.
      See :ref:`here <sec_mad-statistic>` for more information.

   .. macro:: GSL_FILTER_SCALE_IQR

      This option specifies the interquartile range (IQR) scale estimate, defined as the difference between
      the 75th and 25th percentiles of the window :math:`W_i^H`,

      .. only:: not texinfo

         .. math:: S_i = 0.7413 \left( Q_{0.75} - Q_{0.25} \right)

      .. only:: texinfo

         ::

            S_i = 0.7413 ( Q_{0.75} - Q_{0.25} )

      where :math:`Q_p` is the p-quantile of the window :math:`W_i^H`. The idea is to throw away the largest
      and smallest 25% of the window samples (where the outliers would be), and estimate a scale from the middle 50%.
      The factor :math:`0.7413` provides an unbiased estimate of the standard deviation for Gaussian data.

   .. macro:: GSL_FILTER_SCALE_SN

      This option specifies the so-called :math:`S_n` statistic proposed by Croux and Rousseeuw.
      See :ref:`here <sec_Sn-statistic>` for more information.

   .. macro:: GSL_FILTER_SCALE_QN

      This option specifies the so-called :math:`Q_n` statistic proposed by Croux and Rousseeuw.
      See :ref:`here <sec_Qn-statistic>` for more information.

.. warning::

   While the scale estimates defined above are much less sensitive to outliers than the standard deviation,
   they can suffer from an effect called *implosion*. The standard deviation of a window :math:`W_i^H` will be zero
   if and only if all samples in the window are equal. However, it is possible for the MAD of a window
   to be zero even if all the samples in the window are not equal. For example, if :math:`K/2 + 1` or more
   of the :math:`K` samples in the window are equal to some value :math:`x^{*}`, then the window median will
   be equal to :math:`x^{*}`. Consequently, at least :math:`K/2 + 1` of the absolute deviations
   :math:`|x_j - x^{*}|` will be zero, and so the MAD will be zero. In such a case, the Hampel
   filter will act like the standard median filter regardless of the value of :math:`t`. Caution should also
   be exercised if dividing by :math:`S_i`.

.. function:: gsl_filter_impulse_workspace * gsl_filter_impulse_alloc(const size_t K)

   This function initializes a workspace for impulse detection filtering using a symmetric moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(6K)`.

.. function:: void gsl_filter_impulse_free(gsl_filter_impulse_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_impulse(const gsl_filter_end_t endtype, const gsl_filter_scale_t scale_type, const double t, const gsl_vector * x, gsl_vector * y, gsl_vector * xmedian, gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier, gsl_filter_impulse_workspace * w)

   These functions apply an impulse detection filter to the input vector :data:`x`, storing the filtered output
   in :data:`y`. The tuning parameter :math:`t` is provided in :data:`t`.
   The window medians :math:`m_i` are stored in :data:`xmedian` and the :math:`S_i` are stored in :data:`xsigma` on output.
   The number of outliers detected is stored in :data:`noutlier` on output, while
   the locations of flagged outliers are stored in the boolean array :data:`ioutlier`. The input
   :data:`ioutlier` may be :code:`NULL` if not desired. It  is allowed to have :data:`x` = :data:`y` for an
   in-place filter.

Examples
========

Gaussian Example 1
------------------

This example program illustrates the Gaussian filter applied to smoothing a time series of length
:math:`N = 500` with a kernel size of :math:`K = 51`. Three filters are applied with
parameters :math:`\alpha = 0.5, 3, 10`. The results are shown in :numref:`fig_filt-gaussian`.

.. _fig_filt-gaussian:

.. figure:: /images/gaussfilt.png
   :scale: 60%

   Top panel: Gaussian kernels (unnormalized) for :math:`\alpha = 0.5, 3, 10`.
   Bottom panel: Time series (gray) with Gaussian filter output for same :math:`\alpha`
   values.

We see that the filter corresponding to :math:`\alpha = 0.5` applies the most smoothing,
while :math:`\alpha = 10` corresponds to the least amount of smoothing.
The program is given below.

.. include:: examples/gaussfilt.c
   :code:

Gaussian Example 2
------------------

A common application of the Gaussian filter is to detect edges, or sudden jumps, in a noisy
input signal. It is used both for 1D edge detection in time series, as well as 2D edge
detection in images. Here we will examine a noisy time series of length :math:`N = 1000`
with a single edge. The input signal is defined as

.. only:: not texinfo

   .. math:: x(n) = e(n) +
               \left\{
                 \begin{array}{cc}
                   0, & n \le N/2 \\
                   0.5, & n > N/2
                 \end{array}
               \right.

.. only:: texinfo

   ::

      x(n) = e(n) + { 0,   n <= N/2
                    { 0.5, n >  N/2

where :math:`e(n)` is Gaussian random noise. The program smooths the input signal
with order :math:`0,1,` and :math:`2` Gaussian filters of length :math:`K = 61` with
:math:`\alpha = 3`. For comparison, the program also computes finite differences
of the input signal without smoothing. The results are shown in :numref:`fig_filt-gaussian2`.

.. _fig_filt-gaussian2:

.. figure:: /images/gaussfilt2.png
   :scale: 60%

   Top row: original input signal :math:`x(n)` (black) with Gaussian smoothed signal in red.
   Second row: First finite differences of input signal.
   Third row: Input signal smoothed with a first order Gaussian filter.
   Fourth row: Input signal smoothed with a second order Gaussian filter.

The finite difference approximation of the first derivative (second row) shows
the common problem with differentiating a noisy signal. The noise is amplified
and makes it extremely difficult to detect the sharp gradient at sample :math:`500`.
The third row shows the first order Gaussian smoothed signal with a clear peak
at the location of the edge. Alternatively, one could examine the second order
Gaussian smoothed signal (fourth row) and look for zero crossings to determine
the edge location.

The program is given below.

.. include:: examples/gaussfilt2.c
   :code:

Square Wave Signal Example
--------------------------

The following example program illustrates the median filters on a noisy
square wave signal. Median filters are well known for preserving sharp
edges in the input signal while reducing noise. The program constructs
a 5 Hz square wave signal with Gaussian noise added. Then the signal is
filtered with a standard median filter and recursive median filter using
a symmetric window of length :math:`K = 7`. The results are shown in
:numref:`fig_filt-edge`.

.. _fig_filt-edge:

.. figure:: /images/filt_edge.png
   :scale: 60%

   Original time series is in gray. The standard median filter output is in
   green and the recursive median filter output is in red.

Both filters preserve the sharp signal edges while reducing the noise. The
recursive median filter achieves a smoother result than the standard median
filter. The "blocky" nature of the output is characteristic of all median
filters. The program is given below.

.. include:: examples/filt_edge.c
   :code:

Impulse Detection Example
-------------------------

The following example program illustrates the impulse detection filter. First,
it constructs a sinusoid signal of length :math:`N = 1000` with Gaussian noise
added. Then, about 1% of the data are perturbed to represent large outliers. An
impulse detecting filter is applied with a window size :math:`K = 25` and
tuning parameter :math:`t = 4`, using the :math:`Q_n` statistic as the robust
measure of scale. The results are plotted in :numref:`fig_impulse`.

.. _fig_impulse:

.. figure:: /images/impulse.png
   :scale: 60%

   Original time series is in blue, filter output is in green, upper and
   lower intervals for detecting outliers are in red and yellow respectively.
   Detected outliers are marked with squares.

The program is given below.

.. include:: examples/impulse.c
   :code:

References and Further Reading
==============================

The following publications are relevant to the algorithms described
in this chapter,

* F. J. Harris, *On the use of windows for harmonic analysis with the discrete Fourier transform*,
  Proceedings of the IEEE, 66 (1), 1978.

* S-J. Ko, Y-H. Lee, and A. T. Fam. *Efficient implementation of one-dimensional recursive median filters*,
  IEEE transactions on circuits and systems 37.11 (1990): 1447-1450.

* R. K. Pearson and M. Gabbouj, *Nonlinear Digital Filtering with Python: An Introduction*.
  CRC Press, 2015.
