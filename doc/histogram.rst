.. index::
   single: histograms
   single: binning data

**********
Histograms
**********

This chapter describes functions for creating histograms.  Histograms
provide a convenient way of summarizing the distribution of a set of
data. A histogram consists of a set of *bins* which count the number
of events falling into a given range of a continuous variable :math:`x`.
In GSL the bins of a histogram contain floating-point numbers, so they
can be used to record both integer and non-integer distributions.  The
bins can use arbitrary sets of ranges (uniformly spaced bins are the
default).  Both one and two-dimensional histograms are supported.

Once a histogram has been created it can also be converted into a
probability distribution function.  The library provides efficient
routines for selecting random samples from probability distributions.
This can be useful for generating simulations based on real data.

The functions are declared in the header files :file:`gsl_histogram.h`
and :file:`gsl_histogram2d.h`.

The histogram struct
====================

A histogram is defined by the following struct,

.. type:: gsl_histogram

   ============================= ============================================================================
   size_t n                      This is the number of histogram bins
   double * range                The ranges of the bins are stored in an array of :code:`n+1` elements
                                 pointed to by range.
   double * bin                  The counts for each bin are stored in an array of :data:`n` elements
                                 pointed to by :data:`bin`.  The bins are floating-point numbers, so you can
                                 increment them by non-integer values if necessary.
   ============================= ============================================================================

   The range for :code:`bin[i]` is given by :code:`range[i]` to
   :code:`range[i+1]`.  For :math:`n` bins there are :code:`n+1` entries in the
   array :data:`range`.  Each bin is inclusive at the lower end and exclusive
   at the upper end.  Mathematically this means that the bins are defined by
   the following inequality,

   .. only:: not texinfo

      .. math:: \hbox{bin[i] corresponds to range[i]} \le x < \hbox{range[i+1]}

   .. only:: texinfo

      ::

         bin[i] corresponds to range[i] <= x < range[i+1]

   Here is a diagram of the correspondence between ranges and bins on the
   number-line for :math:`x`::

          [ bin[0] )[ bin[1] )[ bin[2] )[ bin[3] )[ bin[4] )
       ---|---------|---------|---------|---------|---------|---  x
        r[0]      r[1]      r[2]      r[3]      r[4]      r[5]

   In this picture the values of the :data:`range` array are denoted by
   :math:`r`.  On the left-hand side of each bin the square bracket
   :code:`[` denotes an inclusive lower bound 
   (:math:`r \le x`),
   and the round parentheses :code:`)` on the right-hand
   side denote an exclusive upper bound (:math:`x < r`).  Thus any samples
   which fall on the upper end of the histogram are excluded.  If you want
   to include this value for the last bin you will need to add an extra bin
   to your histogram.

   The :type:`gsl_histogram` struct and its associated functions are defined
   in the header file :file:`gsl_histogram.h`.

Histogram allocation
====================

The functions for allocating memory to a histogram follow the style of
:func:`malloc` and :func:`free`.  In addition they also perform their own
error checking.  If there is insufficient memory available to allocate a
histogram then the functions call the error handler (with an error
number of :macro:`GSL_ENOMEM`) in addition to returning a null pointer.
Thus if you use the library error handler to abort your program then it
isn't necessary to check every histogram :code:`alloc`.

.. function:: gsl_histogram * gsl_histogram_alloc (size_t n)

   This function allocates memory for a histogram with :data:`n` bins, and
   returns a pointer to a newly created :type:`gsl_histogram` struct.  If
   insufficient memory is available a null pointer is returned and the
   error handler is invoked with an error code of :macro:`GSL_ENOMEM`. The
   bins and ranges are not initialized, and should be prepared using one of
   the range-setting functions below in order to make the histogram ready
   for use.

.. @deftypefun {gsl_histogram *} gsl_histogram_calloc (size_t n)
.. This function allocates memory for a histogram with :data:`n` bins, and
.. returns a pointer to its newly initialized :type:`gsl_histogram` struct.
.. The bins are uniformly spaced with a total range of 
.. @c{$0 \le  x < n$}
.. @math{0 <=  x < n},
.. as shown in the table below.

.. @tex
.. \beforedisplay
.. $$
.. \matrix{
.. \hbox{bin[0]}&\hbox{corresponds to}& 0 \le x < 1\cr
.. \hbox{bin[1]}&\hbox{corresponds to}& 1 \le x < 2\cr
.. \dots&\dots&\dots\cr
.. \hbox{bin[n-1]}&\hbox{corresponds to}&n-1 \le x < n}
.. $$
.. \afterdisplay
.. @end tex
.. @ifinfo
.. @display
.. bin[0] corresponds to 0 <= x < 1
.. bin[1] corresponds to 1 <= x < 2
.. @dots{}
.. bin[n-1] corresponds to n-1 <= x < n
.. @end display
.. @end ifinfo
.. @noindent
.. The bins are initialized to zero so the histogram is ready for use.

.. If insufficient memory is available a null pointer is returned and the
.. error handler is invoked with an error code of :macro:`GSL_ENOMEM`.
.. @end deftypefun

.. @deftypefun {gsl_histogram *} gsl_histogram_calloc_uniform (size_t n, double xmin, double xmax)
.. This function allocates memory for a histogram with :data:`n` uniformly
.. spaced bins from :data:`xmin` to :data:`xmax`, and returns a pointer to the
.. newly initialized :type:`gsl_histogram` struct. 
.. If insufficient memory is available a null pointer is returned and the
.. error handler is invoked with an error code of :macro:`GSL_ENOMEM`.
.. @end deftypefun

.. @deftypefun {gsl_histogram *} gsl_histogram_calloc_range (size_t n, double * range)
.. This function allocates a histogram of size :data:`n` using the @math{n+1}
.. bin ranges specified by the array :data:`range`.
.. @end deftypefun

.. function:: int gsl_histogram_set_ranges (gsl_histogram * h, const double range[], size_t size)

   This function sets the ranges of the existing histogram :data:`h` using
   the array :data:`range` of size :data:`size`.  The values of the histogram
   bins are reset to zero.  The :data:`range` array should contain the
   desired bin limits.  The ranges can be arbitrary, subject to the
   restriction that they are monotonically increasing.

   The following example shows how to create a histogram with logarithmic
   bins with ranges [1,10), [10,100) and [100,1000)::

      gsl_histogram * h = gsl_histogram_alloc (3);

      /* bin[0] covers the range 1 <= x < 10 */
      /* bin[1] covers the range 10 <= x < 100 */
      /* bin[2] covers the range 100 <= x < 1000 */

      double range[4] = { 1.0, 10.0, 100.0, 1000.0 };

      gsl_histogram_set_ranges (h, range, 4);

   Note that the size of the :data:`range` array should be defined to be one
   element bigger than the number of bins.  The additional element is
   required for the upper value of the final bin.

.. function:: int gsl_histogram_set_ranges_uniform (gsl_histogram * h, double xmin, double xmax)

   This function sets the ranges of the existing histogram :data:`h` to cover
   the range :data:`xmin` to :data:`xmax` uniformly.  The values of the
   histogram bins are reset to zero.  The bin ranges are shown in the table
   below,

   .. only:: not texinfo

      .. math::

         \begin{array}{ccc}
           \hbox{bin[0]}&\hbox{corresponds to}& xmin \le  x < xmin + d \\
           \hbox{bin[1]} &\hbox{corresponds to}& xmin + d \le  x < xmin + 2 d \\
           \dots&\dots&\dots \\
           \hbox{bin[n-1]} & \hbox{corresponds to}& xmin + (n-1)d \le  x < xmax
         \end{array}

   .. only:: texinfo

      ::

         bin[0] corresponds to xmin <= x < xmin + d
         bin[1] corresponds to xmin + d <= x < xmin + 2 d
         ......
         bin[n-1] corresponds to xmin + (n-1)d <= x < xmax

   where :math:`d` is the bin spacing, :math:`d = (xmax-xmin)/n`.

.. function:: void gsl_histogram_free (gsl_histogram * h)

   This function frees the histogram :data:`h` and all of the memory
   associated with it.

Copying Histograms
==================

.. function:: int gsl_histogram_memcpy (gsl_histogram * dest, const gsl_histogram * src)

   This function copies the histogram :data:`src` into the pre-existing
   histogram :data:`dest`, making :data:`dest` into an exact copy of :data:`src`.
   The two histograms must be of the same size.

.. function:: gsl_histogram * gsl_histogram_clone (const gsl_histogram * src)

   This function returns a pointer to a newly created histogram which is an
   exact copy of the histogram :data:`src`.

Updating and accessing histogram elements
=========================================

There are two ways to access histogram bins, either by specifying an
:math:`x` coordinate or by using the bin-index directly.  The functions
for accessing the histogram through :math:`x` coordinates use a binary
search to identify the bin which covers the appropriate range.

.. function:: int gsl_histogram_increment (gsl_histogram * h, double x)

   This function updates the histogram :data:`h` by adding one (1.0) to the
   bin whose range contains the coordinate :data:`x`. 

   If :data:`x` lies in the valid range of the histogram then the function
   returns zero to indicate success.  If :data:`x` is less than the lower
   limit of the histogram then the function returns :macro:`GSL_EDOM`, and
   none of bins are modified.  Similarly, if the value of :data:`x` is greater
   than or equal to the upper limit of the histogram then the function
   returns :macro:`GSL_EDOM`, and none of the bins are modified.  The error
   handler is not called, however, since it is often necessary to compute
   histograms for a small range of a larger dataset, ignoring the values
   outside the range of interest.

.. function:: int gsl_histogram_accumulate (gsl_histogram * h, double x, double weight)

   This function is similar to :func:`gsl_histogram_increment` but increases
   the value of the appropriate bin in the histogram :data:`h` by the
   floating-point number :data:`weight`.

.. function:: double gsl_histogram_get (const gsl_histogram * h, size_t i)

   This function returns the contents of the :data:`i`-th bin of the histogram
   :data:`h`.  If :data:`i` lies outside the valid range of indices for the
   histogram then the error handler is called with an error code of
   :macro:`GSL_EDOM` and the function returns 0.

.. function:: int gsl_histogram_get_range (const gsl_histogram * h, size_t i, double * lower, double * upper)

   This function finds the upper and lower range limits of the :data:`i`-th
   bin of the histogram :data:`h`.  If the index :data:`i` is valid then the
   corresponding range limits are stored in :data:`lower` and :data:`upper`.
   The lower limit is inclusive (i.e. events with this coordinate are
   included in the bin) and the upper limit is exclusive (i.e. events with
   the coordinate of the upper limit are excluded and fall in the
   neighboring higher bin, if it exists).  The function returns 0 to
   indicate success.  If :data:`i` lies outside the valid range of indices for
   the histogram then the error handler is called and the function returns
   an error code of :macro:`GSL_EDOM`.

.. function:: double gsl_histogram_max (const gsl_histogram * h)
              double gsl_histogram_min (const gsl_histogram * h)
              size_t gsl_histogram_bins (const gsl_histogram * h)

   These functions return the maximum upper and minimum lower range limits
   and the number of bins of the histogram :data:`h`.  They provide a way of
   determining these values without accessing the :type:`gsl_histogram`
   struct directly.

.. function:: void gsl_histogram_reset (gsl_histogram * h)

   This function resets all the bins in the histogram :data:`h` to zero.

Searching histogram ranges
==========================

The following functions are used by the access and update routines to
locate the bin which corresponds to a given :math:`x` coordinate.

.. function:: int gsl_histogram_find (const gsl_histogram * h, double x, size_t * i)

   This function finds and sets the index :data:`i` to the bin number which
   covers the coordinate :data:`x` in the histogram :data:`h`.  The bin is
   located using a binary search. The search includes an optimization for
   histograms with uniform range, and will return the correct bin
   immediately in this case.  If :data:`x` is found in the range of the
   histogram then the function sets the index :data:`i` and returns
   :macro:`GSL_SUCCESS`.  If :data:`x` lies outside the valid range of the
   histogram then the function returns :macro:`GSL_EDOM` and the error
   handler is invoked.

.. index::
   single: histogram statistics
   single: statistics, from histogram
   single: maximum value, from histogram
   single: minimum value, from histogram

Histogram Statistics
====================

.. function:: double gsl_histogram_max_val (const gsl_histogram * h)

   This function returns the maximum value contained in the histogram bins.

.. function:: size_t gsl_histogram_max_bin (const gsl_histogram * h)

   This function returns the index of the bin containing the maximum
   value. In the case where several bins contain the same maximum value the
   smallest index is returned.

.. function:: double gsl_histogram_min_val (const gsl_histogram * h)

   This function returns the minimum value contained in the histogram bins.

.. function:: size_t gsl_histogram_min_bin (const gsl_histogram * h)

   This function returns the index of the bin containing the minimum
   value. In the case where several bins contain the same minimum value the
   smallest index is returned.

.. index::
   single: mean value, from histogram

.. function:: double gsl_histogram_mean (const gsl_histogram * h)

   This function returns the mean of the histogrammed variable, where the
   histogram is regarded as a probability distribution. Negative bin values
   are ignored for the purposes of this calculation.  The accuracy of the
   result is limited by the bin width.

.. index::
   single: standard deviation, from histogram
   single: variance, from histogram

.. function:: double gsl_histogram_sigma (const gsl_histogram * h)

   This function returns the standard deviation of the histogrammed
   variable, where the histogram is regarded as a probability
   distribution. Negative bin values are ignored for the purposes of this
   calculation. The accuracy of the result is limited by the bin width.

.. function:: double gsl_histogram_sum (const gsl_histogram * h)

   This function returns the sum of all bin values. Negative bin values
   are included in the sum.

Histogram Operations
====================

.. function:: int gsl_histogram_equal_bins_p (const gsl_histogram * h1, const gsl_histogram * h2)

   This function returns 1 if the all of the individual bin
   ranges of the two histograms are identical, and 0
   otherwise.

.. function:: int gsl_histogram_add (gsl_histogram * h1, const gsl_histogram * h2)

   This function adds the contents of the bins in histogram :data:`h2` to the
   corresponding bins of histogram :data:`h1`,  i.e. :math:`h'_1(i) = h_1(i) + h_2(i)`.
   The two histograms must have identical bin ranges.

.. function:: int gsl_histogram_sub (gsl_histogram * h1, const gsl_histogram * h2)

   This function subtracts the contents of the bins in histogram :data:`h2`
   from the corresponding bins of histogram :data:`h1`, i.e. :math:`h'_1(i) = h_1(i) - h_2(i)`.
   The two histograms must have identical bin ranges.

.. function:: int gsl_histogram_mul (gsl_histogram * h1, const gsl_histogram * h2)

   This function multiplies the contents of the bins of histogram :data:`h1`
   by the contents of the corresponding bins in histogram :data:`h2`,
   i.e. :math:`h'_1(i) = h_1(i) * h_2(i)`.  The two histograms must have
   identical bin ranges.

.. function:: int gsl_histogram_div (gsl_histogram * h1, const gsl_histogram * h2)

   This function divides the contents of the bins of histogram :data:`h1` by
   the contents of the corresponding bins in histogram :data:`h2`,
   i.e. :math:`h'_1(i) = h_1(i) / h_2(i)`.  The two histograms must have
   identical bin ranges.

.. function:: int gsl_histogram_scale (gsl_histogram * h, double scale)

   This function multiplies the contents of the bins of histogram :data:`h`
   by the constant :data:`scale`, i.e.
   
   .. only:: not texinfo
   
      .. math:: h'_1(i) = h_1(i) * \hbox{\it scale}

   .. only:: texinfo

      ::

         h'_1(i) = h_1(i) * scale

.. function:: int gsl_histogram_shift (gsl_histogram * h, double offset)

   This function shifts the contents of the bins of histogram :data:`h` by
   the constant :data:`offset`, i.e.
   
   .. only:: not texinfo
   
      .. math:: h'_1(i) = h_1(i) + \hbox{\it offset}

   .. only:: texinfo

      ::

         h'_1(i) = h_1(i) + offset

Reading and writing histograms
==============================

The library provides functions for reading and writing histograms to a file
as binary data or formatted text.

.. function:: int gsl_histogram_fwrite (FILE * stream, const gsl_histogram * h)

   This function writes the ranges and bins of the histogram :data:`h` to the
   stream :data:`stream` in binary format.  The return value is 0 for success
   and :macro:`GSL_EFAILED` if there was a problem writing to the file.  Since
   the data is written in the native binary format it may not be portable
   between different architectures.

.. function:: int gsl_histogram_fread (FILE * stream, gsl_histogram * h)

   This function reads into the histogram :data:`h` from the open stream
   :data:`stream` in binary format.  The histogram :data:`h` must be
   preallocated with the correct size since the function uses the number of
   bins in :data:`h` to determine how many bytes to read.  The return value is
   0 for success and :macro:`GSL_EFAILED` if there was a problem reading from
   the file.  The data is assumed to have been written in the native binary
   format on the same architecture.

.. function:: int gsl_histogram_fprintf (FILE * stream, const gsl_histogram * h, const char * range_format, const char * bin_format)

   This function writes the ranges and bins of the histogram :data:`h`
   line-by-line to the stream :data:`stream` using the format specifiers
   :data:`range_format` and :data:`bin_format`.  These should be one of the
   :code:`%g`, :code:`%e` or :code:`%f` formats for floating point
   numbers.  The function returns 0 for success and :macro:`GSL_EFAILED` if
   there was a problem writing to the file.  The histogram output is
   formatted in three columns, and the columns are separated by spaces,
   like this::

      range[0] range[1] bin[0]
      range[1] range[2] bin[1]
      range[2] range[3] bin[2]
      ....
      range[n-1] range[n] bin[n-1]

   The values of the ranges are formatted using :data:`range_format` and the
   value of the bins are formatted using :data:`bin_format`.  Each line
   contains the lower and upper limit of the range of the bins and the
   value of the bin itself.  Since the upper limit of one bin is the lower
   limit of the next there is duplication of these values between lines but
   this allows the histogram to be manipulated with line-oriented tools.

.. function:: int gsl_histogram_fscanf (FILE * stream, gsl_histogram * h)

   This function reads formatted data from the stream :data:`stream` into the
   histogram :data:`h`.  The data is assumed to be in the three-column format
   used by :func:`gsl_histogram_fprintf`.  The histogram :data:`h` must be
   preallocated with the correct length since the function uses the size of
   :data:`h` to determine how many numbers to read.  The function returns 0
   for success and :macro:`GSL_EFAILED` if there was a problem reading from
   the file.

.. index::
   single: resampling from histograms
   single: sampling from histograms
   single: probability distributions, from histograms

Resampling from histograms
==========================

A histogram made by counting events can be regarded as a measurement of
a probability distribution.  Allowing for statistical error, the height
of each bin represents the probability of an event where the value of
:math:`x` falls in the range of that bin.  The probability distribution
function has the one-dimensional form :math:`p(x)dx` where,

.. math:: p(x) = n_i / (N w_i)

In this equation :math:`n_i` is the number of events in the bin which
contains :math:`x`, :math:`w_i` is the width of the bin and :math:`N` is
the total number of events.  The distribution of events within each bin
is assumed to be uniform.

.. index::
   single: probability distribution, from histogram
   single: sampling from histograms
   single: random sampling from histograms
   single: histograms, random sampling from

The histogram probability distribution struct
=============================================

The probability distribution function for a histogram consists of a set
of *bins* which measure the probability of an event falling into a
given range of a continuous variable :math:`x`. A probability
distribution function is defined by the following struct, which actually
stores the cumulative probability distribution function.  This is the
natural quantity for generating samples via the inverse transform
method, because there is a one-to-one mapping between the cumulative
probability distribution and the range [0,1].  It can be shown that by
taking a uniform random number in this range and finding its
corresponding coordinate in the cumulative probability distribution we
obtain samples with the desired probability distribution.

.. type:: gsl_histogram_pdf

   ================================ =======================================================================
   :code:`size_t n`                 This is the number of bins used to approximate the probability
                                    distribution function. 
   :code:`double * range`           The ranges of the bins are stored in an array of :math:`n + 1`
                                    elements pointed to by :data:`range`.
   :code:`double * sum`             The cumulative probability for the bins is stored in an array of
                                    :data:`n` elements pointed to by :data:`sum`.
   ================================ =======================================================================

The following functions allow you to create a :type:`gsl_histogram_pdf`
struct which represents this probability distribution and generate
random samples from it.

.. function:: gsl_histogram_pdf * gsl_histogram_pdf_alloc (size_t n)

   This function allocates memory for a probability distribution with
   :data:`n` bins and returns a pointer to a newly initialized
   :type:`gsl_histogram_pdf` struct. If insufficient memory is available a
   null pointer is returned and the error handler is invoked with an error
   code of :macro:`GSL_ENOMEM`.

.. function:: int gsl_histogram_pdf_init (gsl_histogram_pdf * p, const gsl_histogram * h)

   This function initializes the probability distribution :data:`p` with
   the contents of the histogram :data:`h`. If any of the bins of :data:`h` are
   negative then the error handler is invoked with an error code of
   :macro:`GSL_EDOM` because a probability distribution cannot contain
   negative values.

.. function:: void gsl_histogram_pdf_free (gsl_histogram_pdf * p)

   This function frees the probability distribution function :data:`p` and
   all of the memory associated with it.

.. function:: double gsl_histogram_pdf_sample (const gsl_histogram_pdf * p, double r)

   This function uses :data:`r`, a uniform random number between zero and
   one, to compute a single random sample from the probability distribution
   :data:`p`.  The algorithm used to compute the sample :math:`s` is given by
   the following formula,

   .. only:: not texinfo

      .. math:: s = \hbox{range}[i] + \delta * (\hbox{range}[i+1] - \hbox{range}[i])

   .. only:: texinfo

      ::

         s = range[i] + delta * (range[i+1] - range[i])

   where :math:`i` is the index which satisfies 
   :math:`sum[i] \le  r < sum[i+1]`
   and :math:`delta` is 
   :math:`(r - sum[i])/(sum[i+1] - sum[i])`.

Example programs for histograms
===============================

The following program shows how to make a simple histogram of a column
of numerical data supplied on :code:`stdin`.  The program takes three
arguments, specifying the upper and lower bounds of the histogram and
the number of bins.  It then reads numbers from :code:`stdin`, one line at
a time, and adds them to the histogram.  When there is no more data to
read it prints out the accumulated histogram using
:func:`gsl_histogram_fprintf`.

.. include:: examples/histogram.c
   :code:

Here is an example of the program in use.  We generate 10000 random
samples from a Cauchy distribution with a width of 30 and histogram
them over the range -100 to 100, using 200 bins::

  $ gsl-randist 0 10000 cauchy 30 
     | gsl-histogram -- -100 100 200 > histogram.dat

:numref:`fig_histogram` shows the familiar shape of the
Cauchy distribution and the fluctuations caused by the finite sample
size.

.. _fig_histogram:

.. figure:: /images/histogram.png
   :scale: 60%

   Histogram output from example program

.. index::
   single: two dimensional histograms
   single: 2D histograms

Two dimensional histograms
==========================

A two dimensional histogram consists of a set of *bins* which count
the number of events falling in a given area of the :math:`(x,y)`
plane.  The simplest way to use a two dimensional histogram is to record
two-dimensional position information, :math:`n(x,y)`.  Another possibility
is to form a *joint distribution* by recording related
variables.  For example a detector might record both the position of an
event (:math:`x`) and the amount of energy it deposited :math:`E`.  These
could be histogrammed as the joint distribution :math:`n(x,E)`.

The 2D histogram struct
=======================

Two dimensional histograms are defined by the following struct,

.. type:: gsl_histogram2d

   =========================== ============================================================================
   :code:`size_t nx, ny`       This is the number of histogram bins in the x and y directions.
   :code:`double * xrange`     The ranges of the bins in the x-direction are stored in an array of
                               :code:`nx + 1` elements pointed to by :data:`xrange`.
   :code:`double * yrange`     The ranges of the bins in the y-direction are stored in an array of
                               :code:`ny + 1` elements pointed to by :data:`yrange`.
   :code:`double * bin`        The counts for each bin are stored in an array pointed to by :data:`bin`.
                               The bins are floating-point numbers, so you can increment them by
                               non-integer values if necessary.  The array :data:`bin` stores the two
                               dimensional array of bins in a single block of memory according to the
                               mapping :code:`bin(i,j)` = :code:`bin[i * ny + j]`.
   =========================== ============================================================================

The range for :code:`bin(i,j)` is given by :code:`xrange[i]` to
:code:`xrange[i+1]` in the x-direction and :code:`yrange[j]` to
:code:`yrange[j+1]` in the y-direction.  Each bin is inclusive at the lower
end and exclusive at the upper end.  Mathematically this means that the
bins are defined by the following inequality,

.. only:: not texinfo

   .. math::

      \begin{array}{cc}
        \hbox{bin(i,j) corresponds to} & \hbox{\it xrange}[i] \le x < \hbox{\it xrange}[i+1] \\
        \hbox{and} & \hbox{\it yrange}[j] \le y < \hbox{\it yrange}[j+1]
      \end{array}

.. only:: texinfo

   ::

      bin(i,j) corresponds to xrange[i] <= x < xrange[i+1]
                          and yrange[j] <= y < yrange[j+1]

Note that any samples which fall on the upper sides of the histogram are
excluded.  If you want to include these values for the side bins you will
need to add an extra row or column to your histogram.

The :type:`gsl_histogram2d` struct and its associated functions are
defined in the header file :file:`gsl_histogram2d.h`.

2D Histogram allocation
=======================

The functions for allocating memory to a 2D histogram follow the style
of :func:`malloc` and :func:`free`.  In addition they also perform their
own error checking.  If there is insufficient memory available to
allocate a histogram then the functions call the error handler (with
an error number of :macro:`GSL_ENOMEM`) in addition to returning a null
pointer.  Thus if you use the library error handler to abort your program
then it isn't necessary to check every 2D histogram :code:`alloc`.

.. function:: gsl_histogram2d * gsl_histogram2d_alloc (size_t nx, size_t ny)

   This function allocates memory for a two-dimensional histogram with
   :data:`nx` bins in the x direction and :data:`ny` bins in the y direction.
   The function returns a pointer to a newly created :type:`gsl_histogram2d`
   struct. If insufficient memory is available a null pointer is returned
   and the error handler is invoked with an error code of
   :macro:`GSL_ENOMEM`. The bins and ranges must be initialized with one of
   the functions below before the histogram is ready for use.

.. @deftypefun {gsl_histogram2d *} gsl_histogram2d_calloc (size_t nx, size_t ny)
.. This function allocates memory for a two-dimensional histogram with
.. :data:`nx` bins in the x direction and :data:`ny` bins in the y
.. direction.  The function returns a pointer to a newly initialized
.. :type:`gsl_histogram2d` struct.  The bins are uniformly spaced with a
.. total range of 
.. @c{$0 \le  x < nx$}
.. @math{0 <= x < nx} in the x-direction and 
.. @c{$0 \le  y < ny$} 
.. @math{0 <=  y < ny} in the y-direction, as shown in the table below.
.. 
.. The bins are initialized to zero so the histogram is ready for use.
.. 
.. If insufficient memory is available a null pointer is returned and the
.. error handler is invoked with an error code of :macro:`GSL_ENOMEM`.
.. @end deftypefun
.. 
.. @deftypefun {gsl_histogram2d *} gsl_histogram2d_calloc_uniform (size_t nx, size_t ny, double xmin, double xmax, double ymin, double ymax)
.. This function allocates a histogram of size :data:`nx`-by-:data:`ny` which
.. uniformly covers the ranges :data:`xmin` to :data:`xmax` and :data:`ymin` to
.. :data:`ymax` in the :math:`x` and :math:`y` directions respectively.
.. @end deftypefun
.. 
.. @deftypefun {gsl_histogram2d *} gsl_histogram2d_calloc_range (size_t nx, size_t ny, double * xrange, double * yrange)
.. This function allocates a histogram of size :data:`nx`-by-:data:`ny` using
.. the @math{nx+1} and @math{ny+1} bin ranges specified by the arrays
.. :data:`xrange` and :data:`yrange`.
.. @end deftypefun

.. function:: int gsl_histogram2d_set_ranges (gsl_histogram2d * h,  const double xrange[], size_t xsize, const double yrange[], size_t ysize)

   This function sets the ranges of the existing histogram :data:`h` using
   the arrays :data:`xrange` and :data:`yrange` of size :data:`xsize` and
   :data:`ysize` respectively.  The values of the histogram bins are reset to
   zero.

.. function:: int gsl_histogram2d_set_ranges_uniform (gsl_histogram2d * h, double xmin, double xmax, double ymin, double ymax)

   This function sets the ranges of the existing histogram :data:`h` to cover
   the ranges :data:`xmin` to :data:`xmax` and :data:`ymin` to :data:`ymax`
   uniformly.  The values of the histogram bins are reset to zero.

.. function:: void gsl_histogram2d_free (gsl_histogram2d * h)

   This function frees the 2D histogram :data:`h` and all of the memory
   associated with it.

Copying 2D Histograms
=====================

.. function:: int gsl_histogram2d_memcpy (gsl_histogram2d * dest, const gsl_histogram2d * src)

   This function copies the histogram :data:`src` into the pre-existing
   histogram :data:`dest`, making :data:`dest` into an exact copy of :data:`src`.
   The two histograms must be of the same size.

.. function:: gsl_histogram2d * gsl_histogram2d_clone (const gsl_histogram2d * src)

   This function returns a pointer to a newly created histogram which is an
   exact copy of the histogram :data:`src`.

Updating and accessing 2D histogram elements
============================================

You can access the bins of a two-dimensional histogram either by
specifying a pair of :math:`(x,y)` coordinates or by using the bin
indices :math:`(i,j)` directly.  The functions for accessing the histogram
through :math:`(x,y)` coordinates use binary searches in the x and y
directions to identify the bin which covers the appropriate range.

.. function:: int gsl_histogram2d_increment (gsl_histogram2d * h, double x, double y)

   This function updates the histogram :data:`h` by adding one (1.0) to the
   bin whose x and y ranges contain the coordinates (:data:`x`, :data:`y`).

   If the point :math:`(x,y)` lies inside the valid ranges of the
   histogram then the function returns zero to indicate success.  If
   :math:`(x,y)` lies outside the limits of the histogram then the
   function returns :macro:`GSL_EDOM`, and none of the bins are modified.  The
   error handler is not called, since it is often necessary to compute
   histograms for a small range of a larger dataset, ignoring any
   coordinates outside the range of interest.

.. function:: int gsl_histogram2d_accumulate (gsl_histogram2d * h, double x, double y, double weight)

   This function is similar to :func:`gsl_histogram2d_increment` but increases
   the value of the appropriate bin in the histogram :data:`h` by the
   floating-point number :data:`weight`.

.. function:: double gsl_histogram2d_get (const gsl_histogram2d * h, size_t i, size_t j)

   This function returns the contents of the (:data:`i`, :data:`j`)-th bin of the
   histogram :data:`h`.  If (:data:`i`, :data:`j`) lies outside the valid range of
   indices for the histogram then the error handler is called with an error
   code of :macro:`GSL_EDOM` and the function returns 0.

.. function:: int gsl_histogram2d_get_xrange (const gsl_histogram2d * h, size_t i, double * xlower, double * xupper)
              int gsl_histogram2d_get_yrange (const gsl_histogram2d * h, size_t j, double * ylower, double * yupper)

   These functions find the upper and lower range limits of the :data:`i`-th
   and :data:`j`-th bins in the x and y directions of the histogram :data:`h`.
   The range limits are stored in :data:`xlower` and :data:`xupper` or
   :data:`ylower` and :data:`yupper`.  The lower limits are inclusive
   (i.e. events with these coordinates are included in the bin) and the
   upper limits are exclusive (i.e. events with the value of the upper
   limit are not included and fall in the neighboring higher bin, if it
   exists).  The functions return 0 to indicate success.  If :data:`i` or
   :data:`j` lies outside the valid range of indices for the histogram then
   the error handler is called with an error code of :macro:`GSL_EDOM`.

.. function:: double gsl_histogram2d_xmax (const gsl_histogram2d * h)
              double gsl_histogram2d_xmin (const gsl_histogram2d * h)
              size_t gsl_histogram2d_nx (const gsl_histogram2d * h)
              double gsl_histogram2d_ymax (const gsl_histogram2d * h)
              double gsl_histogram2d_ymin (const gsl_histogram2d * h)
              size_t gsl_histogram2d_ny (const gsl_histogram2d * h)

   These functions return the maximum upper and minimum lower range limits
   and the number of bins for the x and y directions of the histogram
   :data:`h`.  They provide a way of determining these values without
   accessing the :type:`gsl_histogram2d` struct directly.

.. function:: void gsl_histogram2d_reset (gsl_histogram2d * h)

   This function resets all the bins of the histogram :data:`h` to zero.

Searching 2D histogram ranges
=============================

The following functions are used by the access and update routines to
locate the bin which corresponds to a given :math:`(x,y)` coordinate.

.. function:: int gsl_histogram2d_find (const gsl_histogram2d * h, double x, double y, size_t * i, size_t * j)

   This function finds and sets the indices :data:`i` and :data:`j` to
   the bin which covers the coordinates (:data:`x`, :data:`y`). The bin is
   located using a binary search.  The search includes an optimization for
   histograms with uniform ranges, and will return the correct bin immediately
   in this case. If :math:`(x,y)` is found then the function sets the
   indices (:data:`i`, :data:`j`) and returns :macro:`GSL_SUCCESS`.  If
   :math:`(x,y)` lies outside the valid range of the histogram then the
   function returns :macro:`GSL_EDOM` and the error handler is invoked.

2D Histogram Statistics
=======================

.. function:: double gsl_histogram2d_max_val (const gsl_histogram2d * h)

   This function returns the maximum value contained in the histogram bins.

.. function:: void gsl_histogram2d_max_bin (const gsl_histogram2d * h, size_t * i, size_t * j)

   This function finds the indices of the bin containing the maximum value
   in the histogram :data:`h` and stores the result in (:data:`i`, :data:`j`). In
   the case where several bins contain the same maximum value the first bin
   found is returned.

.. function:: double gsl_histogram2d_min_val (const gsl_histogram2d * h)

   This function returns the minimum value contained in the histogram bins.

.. function:: void gsl_histogram2d_min_bin (const gsl_histogram2d * h, size_t * i, size_t * j)

   This function finds the indices of the bin containing the minimum value
   in the histogram :data:`h` and stores the result in (:data:`i`, :data:`j`). In
   the case where several bins contain the same maximum value the first bin
   found is returned.

.. function:: double gsl_histogram2d_xmean (const gsl_histogram2d * h)

   This function returns the mean of the histogrammed x variable, where the
   histogram is regarded as a probability distribution. Negative bin values
   are ignored for the purposes of this calculation.

.. function:: double gsl_histogram2d_ymean (const gsl_histogram2d * h)

   This function returns the mean of the histogrammed y variable, where the
   histogram is regarded as a probability distribution. Negative bin values
   are ignored for the purposes of this calculation.

.. function:: double gsl_histogram2d_xsigma (const gsl_histogram2d * h)

   This function returns the standard deviation of the histogrammed
   x variable, where the histogram is regarded as a probability
   distribution. Negative bin values are ignored for the purposes of this
   calculation.

.. function:: double gsl_histogram2d_ysigma (const gsl_histogram2d * h)

   This function returns the standard deviation of the histogrammed
   y variable, where the histogram is regarded as a probability
   distribution. Negative bin values are ignored for the purposes of this
   calculation.

.. function:: double gsl_histogram2d_cov (const gsl_histogram2d * h)

   This function returns the covariance of the histogrammed x and y
   variables, where the histogram is regarded as a probability
   distribution. Negative bin values are ignored for the purposes of this
   calculation.

.. function:: double gsl_histogram2d_sum (const gsl_histogram2d * h)

   This function returns the sum of all bin values. Negative bin values
   are included in the sum.

2D Histogram Operations
=======================

.. function:: int gsl_histogram2d_equal_bins_p (const gsl_histogram2d * h1, const gsl_histogram2d * h2)

   This function returns 1 if all the individual bin ranges of the two
   histograms are identical, and 0 otherwise.

.. function:: int gsl_histogram2d_add (gsl_histogram2d * h1, const gsl_histogram2d * h2)

   This function adds the contents of the bins in histogram :data:`h2` to the
   corresponding bins of histogram :data:`h1`,
   i.e. :math:`h'_1(i,j) = h_1(i,j) + h_2(i,j)`.
   The two histograms must have identical bin ranges.

.. function:: int gsl_histogram2d_sub (gsl_histogram2d * h1, const gsl_histogram2d * h2)

   This function subtracts the contents of the bins in histogram :data:`h2` from the
   corresponding bins of histogram :data:`h1`,
   i.e. :math:`h'_1(i,j) = h_1(i,j) - h_2(i,j)`.
   The two histograms must have identical bin ranges.

.. function:: int gsl_histogram2d_mul (gsl_histogram2d * h1, const gsl_histogram2d * h2)

   This function multiplies the contents of the bins of histogram :data:`h1`
   by the contents of the corresponding bins in histogram :data:`h2`,
   i.e. :math:`h'_1(i,j) = h_1(i,j) * h_2(i,j)`.
   The two histograms must have identical bin ranges.

.. function:: int gsl_histogram2d_div (gsl_histogram2d * h1, const gsl_histogram2d * h2)

   This function divides the contents of the bins of histogram :data:`h1`
   by the contents of the corresponding bins in histogram :data:`h2`,
   i.e. :math:`h'_1(i,j) = h_1(i,j) / h_2(i,j)`.
   The two histograms must have identical bin ranges.

.. function:: int gsl_histogram2d_scale (gsl_histogram2d * h, double scale)

   This function multiplies the contents of the bins of histogram :data:`h`
   by the constant :data:`scale`, i.e.
   
   .. only:: not texinfo
   
      .. math:: h'_1(i,j) = h_1(i,j) * \hbox{\it scale}

   .. only:: texinfo

      ::

         h'_1(i,j) = h_1(i,j) scale

.. function:: int gsl_histogram2d_shift (gsl_histogram2d * h, double offset)

   This function shifts the contents of the bins of histogram :data:`h`
   by the constant :data:`offset`, i.e.
   
   .. only:: not texinfo
   
      .. math:: h'_1(i,j) = h_1(i,j) + \hbox{\it offset}

   .. only:: texinfo

      ::

         h'_1(i,j) = h_1(i,j) + offset

Reading and writing 2D histograms
=================================

The library provides functions for reading and writing two dimensional
histograms to a file as binary data or formatted text.

.. function:: int gsl_histogram2d_fwrite (FILE * stream, const gsl_histogram2d * h)

   This function writes the ranges and bins of the histogram :data:`h` to the
   stream :data:`stream` in binary format.  The return value is 0 for success
   and :macro:`GSL_EFAILED` if there was a problem writing to the file.  Since
   the data is written in the native binary format it may not be portable
   between different architectures.

.. function:: int gsl_histogram2d_fread (FILE * stream, gsl_histogram2d * h)

   This function reads into the histogram :data:`h` from the stream
   :data:`stream` in binary format.  The histogram :data:`h` must be
   preallocated with the correct size since the function uses the number of
   x and y bins in :data:`h` to determine how many bytes to read.  The return
   value is 0 for success and :macro:`GSL_EFAILED` if there was a problem
   reading from the file.  The data is assumed to have been written in the
   native binary format on the same architecture.

.. function:: int gsl_histogram2d_fprintf (FILE * stream, const gsl_histogram2d * h, const char * range_format, const char * bin_format)

   This function writes the ranges and bins of the histogram :data:`h`
   line-by-line to the stream :data:`stream` using the format specifiers
   :data:`range_format` and :data:`bin_format`.  These should be one of the
   :code:`%g`, :code:`%e` or :code:`%f` formats for floating point
   numbers.  The function returns 0 for success and :macro:`GSL_EFAILED` if
   there was a problem writing to the file.  The histogram output is
   formatted in five columns, and the columns are separated by spaces,
   like this::

      xrange[0] xrange[1] yrange[0] yrange[1] bin(0,0)
      xrange[0] xrange[1] yrange[1] yrange[2] bin(0,1)
      xrange[0] xrange[1] yrange[2] yrange[3] bin(0,2)
      ....
      xrange[0] xrange[1] yrange[ny-1] yrange[ny] bin(0,ny-1)

      xrange[1] xrange[2] yrange[0] yrange[1] bin(1,0)
      xrange[1] xrange[2] yrange[1] yrange[2] bin(1,1)
      xrange[1] xrange[2] yrange[1] yrange[2] bin(1,2)
      ....
      xrange[1] xrange[2] yrange[ny-1] yrange[ny] bin(1,ny-1)

      ....

      xrange[nx-1] xrange[nx] yrange[0] yrange[1] bin(nx-1,0)
      xrange[nx-1] xrange[nx] yrange[1] yrange[2] bin(nx-1,1)
      xrange[nx-1] xrange[nx] yrange[1] yrange[2] bin(nx-1,2)
      ....
      xrange[nx-1] xrange[nx] yrange[ny-1] yrange[ny] bin(nx-1,ny-1)

   Each line contains the lower and upper limits of the bin and the
   contents of the bin.  Since the upper limits of the each bin are the
   lower limits of the neighboring bins there is duplication of these
   values but this allows the histogram to be manipulated with
   line-oriented tools.

.. function:: int gsl_histogram2d_fscanf (FILE * stream, gsl_histogram2d * h)

   This function reads formatted data from the stream :data:`stream` into the
   histogram :data:`h`.  The data is assumed to be in the five-column format
   used by :func:`gsl_histogram2d_fprintf`.  The histogram :data:`h` must be
   preallocated with the correct lengths since the function uses the sizes
   of :data:`h` to determine how many numbers to read.  The function returns 0
   for success and :macro:`GSL_EFAILED` if there was a problem reading from
   the file.

Resampling from 2D histograms
=============================

As in the one-dimensional case, a two-dimensional histogram made by
counting events can be regarded as a measurement of a probability
distribution.  Allowing for statistical error, the height of each bin
represents the probability of an event where (:math:`x`, :math:`y`) falls in
the range of that bin.  For a two-dimensional histogram the probability
distribution takes the form :math:`p(x,y) dx dy` where,

.. math:: p(x,y) = n_{ij} / (N A_{ij})

In this equation :math:`n_{ij}`
is the number of events in the bin which
contains :math:`(x,y)`, :math:`A_{ij}`
is the area of the bin and :math:`N` is
the total number of events.  The distribution of events within each bin
is assumed to be uniform.

.. type:: gsl_histogram2d_pdf

   ============================= ===========================================================================
   :code:`size_t nx, ny`         This is the number of histogram bins used to approximate the probability
                                 distribution function in the x and y directions.
   :code:`double * xrange`       The ranges of the bins in the x-direction are stored in an array of
                                 :code:`nx + 1` elements pointed to by :data:`xrange`.
   :code:`double * yrange`       The ranges of the bins in the y-direction are stored in an array of
                                 :code:`ny + 1` pointed to by :data:`yrange`.
   :code:`double * sum`          The cumulative probability for the bins is stored in an array of
                                 :data:`nx` * :data:`ny` elements pointed to by :data:`sum`.
   ============================= ===========================================================================

The following functions allow you to create a :type:`gsl_histogram2d_pdf`
struct which represents a two dimensional probability distribution and
generate random samples from it.

.. function:: gsl_histogram2d_pdf * gsl_histogram2d_pdf_alloc (size_t nx, size_t ny)

   This function allocates memory for a two-dimensional probability
   distribution of size :data:`nx`-by-:data:`ny` and returns a pointer to a
   newly initialized :type:`gsl_histogram2d_pdf` struct. If insufficient
   memory is available a null pointer is returned and the error handler is
   invoked with an error code of :macro:`GSL_ENOMEM`.

.. function:: int gsl_histogram2d_pdf_init (gsl_histogram2d_pdf * p, const gsl_histogram2d * h)

   This function initializes the two-dimensional probability distribution
   calculated :data:`p` from the histogram :data:`h`.  If any of the bins of
   :data:`h` are negative then the error handler is invoked with an error
   code of :macro:`GSL_EDOM` because a probability distribution cannot
   contain negative values.

.. function:: void gsl_histogram2d_pdf_free (gsl_histogram2d_pdf * p)

   This function frees the two-dimensional probability distribution
   function :data:`p` and all of the memory associated with it.

.. function:: int gsl_histogram2d_pdf_sample (const gsl_histogram2d_pdf * p, double r1, double r2, double * x, double * y)

   This function uses two uniform random numbers between zero and one,
   :data:`r1` and :data:`r2`, to compute a single random sample from the
   two-dimensional probability distribution :data:`p`.

Example programs for 2D histograms
==================================

This program demonstrates two features of two-dimensional histograms.
First a 10-by-10 two-dimensional histogram is created with x and y running
from 0 to 1.  Then a few sample points are added to the histogram, at
(0.3,0.3) with a height of 1, at (0.8,0.1) with a height of 5 and at
(0.7,0.9) with a height of 0.5.  This histogram with three events is
used to generate a random sample of 1000 simulated events, which are
printed out.

.. include:: examples/histogram2d.c
   :code:

The following plot shows the distribution of the simulated events.  Using
a higher resolution grid we can see the original underlying histogram
and also the statistical fluctuations caused by the events being
uniformly distributed over the area of the original bins.

.. figure:: /images/histogram2d.png
   :scale: 60%

   Distribution of simulated events from example program
