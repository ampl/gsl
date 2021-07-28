.. index::
   single: FFT
   single: Fast Fourier Transforms, see FFT
   single: Fourier Transforms, see FFT
   single: Discrete Fourier Transforms, see FFT
   single: DFTs, see FFT

******************************
Fast Fourier Transforms (FFTs)
******************************

.. include:: include.rst

This chapter describes functions for performing Fast Fourier Transforms
(FFTs).  The library includes radix-2 routines (for lengths which are a
power of two) and mixed-radix routines (which work for any length).  For
efficiency there are separate versions of the routines for real data and
for complex data.  The mixed-radix routines are a reimplementation of the
|fftpack| library of Paul Swarztrauber.  Fortran code for |fftpack| is
available on Netlib (|fftpack| also includes some routines for sine and
cosine transforms but these are currently not available in GSL).  For
details and derivations of the underlying algorithms consult the
document "GSL FFT Algorithms" (see :ref:`References and Further Reading <fft-references>`)

.. index:: FFT mathematical definition

Mathematical Definitions
========================

Fast Fourier Transforms are efficient algorithms for
calculating the discrete Fourier transform (DFT),

.. math:: x_j = \sum_{k=0}^{n-1} z_k \exp(-2 \pi i j k / n) 

The DFT usually arises as an approximation to the continuous Fourier
transform when functions are sampled at discrete intervals in space or
time.  The naive evaluation of the discrete Fourier transform is a
matrix-vector multiplication :math:`W\vec{z}`.
A general matrix-vector multiplication takes
:math:`O(n^2)` operations for :math:`n` data-points.  Fast Fourier
transform algorithms use a divide-and-conquer strategy to factorize the
matrix :math:`W` into smaller sub-matrices, corresponding to the integer
factors of the length :math:`n`.  If :math:`n` can be factorized into a
product of integers :math:`f_1 f_2 \ldots f_m`
then the DFT can be computed in :math:`O(n \sum f_i)`
operations.  For a radix-2 FFT this gives an operation count of
:math:`O(n \log_2 n)`.

All the FFT functions offer three types of transform: forwards, inverse
and backwards, based on the same mathematical definitions.  The
definition of the *forward Fourier transform*,
:math:`x = \hbox{FFT}(z)`, is,

.. math:: x_j = \sum_{k=0}^{n-1} z_k \exp(-2 \pi i j k / n) 

and the definition of the *inverse Fourier transform*,
:math:`x = \hbox{IFFT}(z)`, is,

.. math:: z_j = {1 \over n} \sum_{k=0}^{n-1} x_k \exp(2 \pi i j k / n).

The factor of :math:`1/n` makes this a true inverse.  For example, a call
to :func:`gsl_fft_complex_forward` followed by a call to
:func:`gsl_fft_complex_inverse` should return the original data (within
numerical errors).

In general there are two possible choices for the sign of the
exponential in the transform/ inverse-transform pair. GSL follows the
same convention as |fftpack|, using a negative exponential for the forward
transform.  The advantage of this convention is that the inverse
transform recreates the original function with simple Fourier
synthesis.  Numerical Recipes uses the opposite convention, a positive
exponential in the forward transform.

The *backwards FFT* is simply our terminology for an unscaled
version of the inverse FFT,

.. math:: z^{backwards}_j = \sum_{k=0}^{n-1} x_k \exp(2 \pi i j k / n)

When the overall scale of the result is unimportant it is often
convenient to use the backwards FFT instead of the inverse to save
unnecessary divisions.

.. index::
   single: FFT, complex data

Overview of complex data FFTs
=============================

The inputs and outputs for the complex FFT routines are *packed arrays*
of floating point numbers.  In a packed array the real and
imaginary parts of each complex number are placed in alternate
neighboring elements.  For example, the following definition of a packed
array of length 6::

  double x[3*2];
  gsl_complex_packed_array data = x;

can be used to hold an array of three complex numbers, :code:`z[3]`, in
the following way::

  data[0] = Re(z[0])
  data[1] = Im(z[0])
  data[2] = Re(z[1])
  data[3] = Im(z[1])
  data[4] = Re(z[2])
  data[5] = Im(z[2])

The array indices for the data have the same ordering as those
in the definition of the DFT---i.e. there are no index transformations
or permutations of the data.

A *stride* parameter allows the user to perform transforms on the
elements :code:`z[stride*i]` instead of :code:`z[i]`.  A stride greater
than 1 can be used to take an in-place FFT of the column of a matrix. A
stride of 1 accesses the array without any additional spacing between
elements.  

To perform an FFT on a vector argument, such as :code:`gsl_vector_complex * v`,
use the following definitions (or their equivalents) when calling
the functions described in this chapter::

  gsl_complex_packed_array data = v->data;
  size_t stride = v->stride;
  size_t n = v->size;

For physical applications it is important to remember that the index
appearing in the DFT does not correspond directly to a physical
frequency.  If the time-step of the DFT is :math:`\Delta` then the
frequency-domain includes both positive and negative frequencies,
ranging from :math:`-1/(2\Delta)` through 0 to :math:`+1/(2\Delta)`.  The
positive frequencies are stored from the beginning of the array up to
the middle, and the negative frequencies are stored backwards from the
end of the array.

Here is a table which shows the layout of the array :data:`data`, and the
correspondence between the time-domain data :math:`z`, and the
frequency-domain data :math:`x`::

  index    z               x = FFT(z)

  0        z(t = 0)        x(f = 0)
  1        z(t = 1)        x(f = 1/(n Delta))
  2        z(t = 2)        x(f = 2/(n Delta))
  .        ........        ..................
  n/2      z(t = n/2)      x(f = +1/(2 Delta),
                                 -1/(2 Delta))
  .        ........        ..................
  n-3      z(t = n-3)      x(f = -3/(n Delta))
  n-2      z(t = n-2)      x(f = -2/(n Delta))
  n-1      z(t = n-1)      x(f = -1/(n Delta))

When :math:`n` is even the location :math:`n/2` contains the most positive
and negative frequencies (:math:`+1/(2 \Delta)`, :math:`-1/(2 \Delta)`)
which are equivalent.  If :math:`n` is odd then general structure of the
table above still applies, but :math:`n/2` does not appear.

.. index::
   single: FFT of complex data, radix-2 algorithm
   single: Radix-2 FFT, complex data

Radix-2 FFT routines for complex data
=====================================

The radix-2 algorithms described in this section are simple and compact,
although not necessarily the most efficient.  They use the Cooley-Tukey
algorithm to compute in-place complex FFTs for lengths which are a power
of 2---no additional storage is required.  The corresponding
self-sorting mixed-radix routines offer better performance at the
expense of requiring additional working space.

All the functions described in this section are declared in the header file :file:`gsl_fft_complex.h`.

.. function:: int gsl_fft_complex_radix2_forward (gsl_complex_packed_array data, size_t stride, size_t n)
              int gsl_fft_complex_radix2_transform (gsl_complex_packed_array data, size_t stride, size_t n, gsl_fft_direction sign)
              int gsl_fft_complex_radix2_backward (gsl_complex_packed_array data, size_t stride, size_t n)
              int gsl_fft_complex_radix2_inverse (gsl_complex_packed_array data, size_t stride, size_t n)

   These functions compute forward, backward and inverse FFTs of length
   :data:`n` with stride :data:`stride`, on the packed complex array :data:`data`
   using an in-place radix-2 decimation-in-time algorithm.  The length of
   the transform :data:`n` is restricted to powers of two.  For the
   :code:`transform` version of the function the :data:`sign` argument can be
   either :code:`forward` (:math:`-1`) or :code:`backward` (:math:`+1`).

   The functions return a value of :macro:`GSL_SUCCESS` if no errors were
   detected, or :macro:`GSL_EDOM` if the length of the data :data:`n` is not a
   power of two.

.. function:: int gsl_fft_complex_radix2_dif_forward (gsl_complex_packed_array data, size_t stride, size_t n)
              int gsl_fft_complex_radix2_dif_transform (gsl_complex_packed_array data, size_t stride, size_t n, gsl_fft_direction sign)
              int gsl_fft_complex_radix2_dif_backward (gsl_complex_packed_array data, size_t stride, size_t n)
              int gsl_fft_complex_radix2_dif_inverse (gsl_complex_packed_array data, size_t stride, size_t n)

   These are decimation-in-frequency versions of the radix-2 FFT functions.

Here is an example program which computes the FFT of a short pulse in a
sample of length 128.  To make the resulting Fourier transform real the
pulse is defined for equal positive and negative times (:math:`-10 \dots 10`),
where the negative times wrap around the end of the array.

.. include:: examples/fft.c
   :code:

Note that we have assumed that the program is using the default error
handler (which calls :func:`abort` for any errors).  If you are not using
a safe error handler you would need to check the return status of
:func:`gsl_fft_complex_radix2_forward`.

The transformed data is rescaled by :math:`1/\sqrt n` so that it fits on
the same plot as the input.  Only the real part is shown, by the choice
of the input data the imaginary part is zero.  Allowing for the
wrap-around of negative times at :math:`t=128`, and working in units of
:math:`k/n`, the DFT approximates the continuum Fourier transform, giving
a modulated sine function.

.. math:: \int_{-a}^{+a} e^{-2 \pi i k x} dx = {\sin(2\pi k a) \over\pi k}

The output of the example program is plotted in :numref:`fig_fft-complex-radix2`.

.. _fig_fft-complex-radix2:

.. figure:: /images/fft-complex-radix2.png
   :scale: 60%

   A pulse and its discrete Fourier transform, output from
   the example program.

.. index::
   single: FFT of complex data, mixed-radix algorithm
   single: Mixed-radix FFT, complex data

Mixed-radix FFT routines for complex data
=========================================

This section describes mixed-radix FFT algorithms for complex data.  The
mixed-radix functions work for FFTs of any length.  They are a
reimplementation of Paul Swarztrauber's Fortran |fftpack| library.
The theory is explained in the review article "Self-sorting
Mixed-radix FFTs" by Clive Temperton.  The routines here use the same
indexing scheme and basic algorithms as |fftpack|.

The mixed-radix algorithm is based on sub-transform modules---highly
optimized small length FFTs which are combined to create larger FFTs.
There are efficient modules for factors of 2, 3, 4, 5, 6 and 7.  The
modules for the composite factors of 4 and 6 are faster than combining
the modules for :math:`2*2` and :math:`2*3`.

For factors which are not implemented as modules there is a fall-back to
a general length-:math:`n` module which uses Singleton's method for
efficiently computing a DFT. This module is :math:`O(n^2)`, and slower
than a dedicated module would be but works for any length :math:`n`.  Of
course, lengths which use the general length-:math:`n` module will still
be factorized as much as possible.  For example, a length of 143 will be
factorized into :math:`11*13`.  Large prime factors are the worst case
scenario, e.g. as found in :math:`n=2*3*99991`, and should be avoided
because their :math:`O(n^2)` scaling will dominate the run-time (consult
the document "GSL FFT Algorithms" included in the GSL distribution
if you encounter this problem).

The mixed-radix initialization function :func:`gsl_fft_complex_wavetable_alloc`
returns the list of factors chosen by the library for a given length
:math:`n`.  It can be used to check how well the length has been
factorized, and estimate the run-time.  To a first approximation the
run-time scales as :math:`n \sum f_i`, where the :math:`f_i` are the
factors of :math:`n`.  For programs under user control you may wish to
issue a warning that the transform will be slow when the length is
poorly factorized.  If you frequently encounter data lengths which
cannot be factorized using the existing small-prime modules consult
"GSL FFT Algorithms" for details on adding support for other
factors.

.. First, the space for the trigonometric lookup tables and scratch area is
.. allocated by a call to one of the :code:`alloc` functions.  We
.. call the combination of factorization, scratch space and trigonometric
.. lookup arrays a *wavetable*.  It contains the sine and cosine
.. waveforms for the all the frequencies that will be used in the FFT.

.. The wavetable is initialized by a call to the corresponding :code:`init`
.. function.  It factorizes the data length, using the implemented
.. subtransforms as preferred factors wherever possible.  The trigonometric
.. lookup table for the chosen factorization is also computed.

.. An FFT is computed by a call to one of the :code:`forward`,
.. :code:`backward` or :code:`inverse` functions, with the data, length and
.. wavetable as arguments.

All the functions described in this section are declared in the header
file :file:`gsl_fft_complex.h`.

.. function:: gsl_fft_complex_wavetable * gsl_fft_complex_wavetable_alloc (size_t n)

   This function prepares a trigonometric lookup table for a complex FFT of
   length :data:`n`. The function returns a pointer to the newly allocated
   :type:`gsl_fft_complex_wavetable` if no errors were detected, and a null
   pointer in the case of error.  The length :data:`n` is factorized into a
   product of subtransforms, and the factors and their trigonometric
   coefficients are stored in the wavetable. The trigonometric coefficients
   are computed using direct calls to :code:`sin` and :code:`cos`, for
   accuracy.  Recursion relations could be used to compute the lookup table
   faster, but if an application performs many FFTs of the same length then
   this computation is a one-off overhead which does not affect the final
   throughput.

   The wavetable structure can be used repeatedly for any transform of the
   same length.  The table is not modified by calls to any of the other FFT
   functions.  The same wavetable can be used for both forward and backward
   (or inverse) transforms of a given length.

.. function:: void gsl_fft_complex_wavetable_free (gsl_fft_complex_wavetable * wavetable)

   This function frees the memory associated with the wavetable
   :data:`wavetable`.  The wavetable can be freed if no further FFTs of the
   same length will be needed.

These functions operate on a :type:`gsl_fft_complex_wavetable` structure
which contains internal parameters for the FFT.  It is not necessary to
set any of the components directly but it can sometimes be useful to
examine them.  For example, the chosen factorization of the FFT length
is given and can be used to provide an estimate of the run-time or
numerical error. The wavetable structure is declared in the header file
:file:`gsl_fft_complex.h`.

.. type:: gsl_fft_complex_wavetable

   This is a structure that holds the factorization and trigonometric
   lookup tables for the mixed radix fft algorithm.  It has the following
   components:

   ================================= ==============================================================================================
   :code:`size_t n`                  This is the number of complex data points
   :code:`size_t nf`                 This is the number of factors that the length :code:`n` was decomposed into.
   :code:`size_t factor[64]`         This is the array of factors.  Only the first :code:`nf` elements are used. 
   :code:`gsl_complex * trig`        This is a pointer to a preallocated trigonometric lookup table of :code:`n` complex elements.
   :code:`gsl_complex * twiddle[64]` This is an array of pointers into :code:`trig`, giving the twiddle factors for each pass.
   ================================= ==============================================================================================

.. (FIXME: factor[64] is a fixed length array and therefore probably in
.. violation of the GNU Coding Standards).

.. type:: gsl_fft_complex_workspace

   The mixed radix algorithms require additional working space to hold
   the intermediate steps of the transform.

.. function:: gsl_fft_complex_workspace * gsl_fft_complex_workspace_alloc (size_t n)

   This function allocates a workspace for a complex transform of length
   :data:`n`.

.. function:: void gsl_fft_complex_workspace_free (gsl_fft_complex_workspace * workspace)

   This function frees the memory associated with the workspace
   :data:`workspace`. The workspace can be freed if no further FFTs of the
   same length will be needed.

The following functions compute the transform,

.. function:: int gsl_fft_complex_forward (gsl_complex_packed_array data, size_t stride, size_t n, const gsl_fft_complex_wavetable * wavetable, gsl_fft_complex_workspace * work)
              int gsl_fft_complex_transform (gsl_complex_packed_array data, size_t stride, size_t n, const gsl_fft_complex_wavetable * wavetable, gsl_fft_complex_workspace * work, gsl_fft_direction sign)
              int gsl_fft_complex_backward (gsl_complex_packed_array data, size_t stride, size_t n, const gsl_fft_complex_wavetable * wavetable, gsl_fft_complex_workspace * work)
              int gsl_fft_complex_inverse (gsl_complex_packed_array data, size_t stride, size_t n, const gsl_fft_complex_wavetable * wavetable, gsl_fft_complex_workspace * work)

   These functions compute forward, backward and inverse FFTs of length
   :data:`n` with stride :data:`stride`, on the packed complex array
   :data:`data`, using a mixed radix decimation-in-frequency algorithm.
   There is no restriction on the length :data:`n`.  Efficient modules are
   provided for subtransforms of length 2, 3, 4, 5, 6 and 7.  Any remaining
   factors are computed with a slow, :math:`O(n^2)`, general-:math:`n`
   module. The caller must supply a :data:`wavetable` containing the
   trigonometric lookup tables and a workspace :data:`work`.  For the
   :code:`transform` version of the function the :data:`sign` argument can be
   either :code:`forward` (:math:`-1`) or :code:`backward` (:math:`+1`).

   The functions return a value of :code:`0` if no errors were detected. The
   following :data:`gsl_errno` conditions are defined for these functions:

   =================================== =========================================================================================================
   :macro:`GSL_EDOM`                   The length of the data :data:`n` is not a positive integer (i.e. :data:`n` is zero).
   :macro:`GSL_EINVAL`                 The length of the data :data:`n` and the length used to compute the given :data:`wavetable` do not match.
   =================================== =========================================================================================================

Here is an example program which computes the FFT of a short pulse in a
sample of length 630 (:math:`=2*3*3*5*7`) using the mixed-radix
algorithm.

.. include:: examples/fftmr.c
   :code:

Note that we have assumed that the program is using the default
:code:`gsl` error handler (which calls :func:`abort` for any errors).  If
you are not using a safe error handler you would need to check the
return status of all the :code:`gsl` routines.

.. index:: FFT of real data

Overview of real data FFTs
==========================

The functions for real data are similar to those for complex data.
However, there is an important difference between forward and inverse
transforms.  The Fourier transform of a real sequence is not real.  It is
a complex sequence with a special symmetry:

.. math:: z_k = z_{n-k}^*

A sequence with this symmetry is called *conjugate-complex* or
*half-complex*.  This different structure requires different
storage layouts for the forward transform (from real to half-complex)
and inverse transform (from half-complex back to real).  As a
consequence the routines are divided into two sets: functions in
:code:`gsl_fft_real` which operate on real sequences and functions in
:code:`gsl_fft_halfcomplex` which operate on half-complex sequences.

Functions in :code:`gsl_fft_real` compute the frequency coefficients of a
real sequence.  The half-complex coefficients :math:`c` of a real sequence
:math:`x` are given by Fourier analysis,

.. math:: c_k = \sum_{j=0}^{n-1} x_j \exp(-2 \pi i j k /n)

Functions in :code:`gsl_fft_halfcomplex` compute inverse or backwards
transforms.  They reconstruct real sequences by Fourier synthesis from
their half-complex frequency coefficients, :math:`c`,

.. math:: x_j = {1 \over n} \sum_{k=0}^{n-1} c_k \exp(2 \pi i j k /n)

The symmetry of the half-complex sequence implies that only half of the
complex numbers in the output need to be stored.  The remaining half can
be reconstructed using the half-complex symmetry condition. This works
for all lengths, even and odd---when the length is even the middle value
where :math:`k=n/2` is also real.  Thus only :data:`n` real numbers are
required to store the half-complex sequence, and the transform of a real
sequence can be stored in the same size array as the original data.

The precise storage arrangements depend on the algorithm, and are
different for radix-2 and mixed-radix routines.  The radix-2 function
operates in-place, which constrains the locations where each element can
be stored.  The restriction forces real and imaginary parts to be stored
far apart.  The mixed-radix algorithm does not have this restriction, and
it stores the real and imaginary parts of a given term in neighboring
locations (which is desirable for better locality of memory accesses).

.. index::
   single: FFT of real data, radix-2 algorithm
   single: Radix-2 FFT for real data

Radix-2 FFT routines for real data
==================================

This section describes radix-2 FFT algorithms for real data.  They use
the Cooley-Tukey algorithm to compute in-place FFTs for lengths which
are a power of 2. 

The radix-2 FFT functions for real data are declared in the header files
:file:`gsl_fft_real.h` 

.. function:: int gsl_fft_real_radix2_transform (double data[], size_t stride, size_t n)

   This function computes an in-place radix-2 FFT of length :data:`n` and
   stride :data:`stride` on the real array :data:`data`.  The output is a
   half-complex sequence, which is stored in-place.  The arrangement of the
   half-complex terms uses the following scheme: for :math:`k < n/2` the
   real part of the :math:`k`-th term is stored in location :math:`k`, and
   the corresponding imaginary part is stored in location :math:`n-k`.  Terms
   with :math:`k > n/2` can be reconstructed using the symmetry 
   :math:`z_k = z^*_{n-k}`.
   The terms for :math:`k=0` and :math:`k=n/2` are both purely
   real, and count as a special case.  Their real parts are stored in
   locations :math:`0` and :math:`n/2` respectively, while their imaginary
   parts which are zero are not stored.

   The following table shows the correspondence between the output
   :data:`data` and the equivalent results obtained by considering the input
   data as a complex sequence with zero imaginary part (assuming :data:`stride` = 1})::

      complex[0].real    =    data[0] 
      complex[0].imag    =    0 
      complex[1].real    =    data[1] 
      complex[1].imag    =    data[n-1]
      ...............         ................
      complex[k].real    =    data[k]
      complex[k].imag    =    data[n-k] 
      ...............         ................
      complex[n/2].real  =    data[n/2]
      complex[n/2].imag  =    0
      ...............         ................
      complex[k'].real   =    data[k]        k' = n - k
      complex[k'].imag   =   -data[n-k] 
      ...............         ................
      complex[n-1].real  =    data[1]
      complex[n-1].imag  =   -data[n-1]

   Note that the output data can be converted into the full complex
   sequence using the function :func:`gsl_fft_halfcomplex_radix2_unpack`
   described below.

The radix-2 FFT functions for halfcomplex data are declared in the
header file :file:`gsl_fft_halfcomplex.h`.

.. function:: int gsl_fft_halfcomplex_radix2_inverse (double data[], size_t stride, size_t n)
              int gsl_fft_halfcomplex_radix2_backward (double data[], size_t stride, size_t n)

   These functions compute the inverse or backwards in-place radix-2 FFT of
   length :data:`n` and stride :data:`stride` on the half-complex sequence
   :data:`data` stored according the output scheme used by
   :func:`gsl_fft_real_radix2`.  The result is a real array stored in natural
   order.

.. function:: int gsl_fft_halfcomplex_radix2_unpack (const double halfcomplex_coefficient[], gsl_complex_packed_array complex_coefficient, size_t stride, size_t n)

   This function converts :data:`halfcomplex_coefficient`, an array of
   half-complex coefficients as returned by :func:`gsl_fft_real_radix2_transform`, into an ordinary complex array, :data:`complex_coefficient`.  It fills in the
   complex array using the symmetry :math:`z_k = z_{n-k}^*`
   to reconstruct the redundant elements.  The algorithm for the conversion
   is::

      complex_coefficient[0].real = halfcomplex_coefficient[0];
      complex_coefficient[0].imag = 0.0;

      for (i = 1; i < n - i; i++)
        {
          double hc_real = halfcomplex_coefficient[i*stride];
          double hc_imag = halfcomplex_coefficient[(n-i)*stride];
          complex_coefficient[i*stride].real = hc_real;
          complex_coefficient[i*stride].imag = hc_imag;
          complex_coefficient[(n - i)*stride].real = hc_real;
          complex_coefficient[(n - i)*stride].imag = -hc_imag;
        }

      if (i == n - i)
        {
          complex_coefficient[i*stride].real = halfcomplex_coefficient[(n - 1)*stride];
          complex_coefficient[i*stride].imag = 0.0;
        }

.. index::
   single: FFT of real data, mixed-radix algorithm
   single: Mixed-radix FFT, real data

Mixed-radix FFT routines for real data
======================================

This section describes mixed-radix FFT algorithms for real data.  The
mixed-radix functions work for FFTs of any length.  They are a
reimplementation of the real-FFT routines in the Fortran |fftpack| library
by Paul Swarztrauber.  The theory behind the algorithm is explained in
the article "Fast Mixed-Radix Real Fourier Transforms" by Clive
Temperton.  The routines here use the same indexing scheme and basic
algorithms as |fftpack|.

The functions use the |fftpack| storage convention for half-complex
sequences.  In this convention the half-complex transform of a real
sequence is stored with frequencies in increasing order, starting at
zero, with the real and imaginary parts of each frequency in neighboring
locations.  When a value is known to be real the imaginary part is not
stored.  The imaginary part of the zero-frequency component is never
stored.  It is known to be zero (since the zero frequency component is
simply the sum of the input data (all real)).  For a sequence of even
length the imaginary part of the frequency :math:`n/2` is not stored
either, since the symmetry :math:`z_k = z_{n-k}^*`
implies that this is purely real too.

The storage scheme is best shown by some examples.  The table below
shows the output for an odd-length sequence, :math:`n=5`.  The two columns
give the correspondence between the 5 values in the half-complex
sequence returned by :func:`gsl_fft_real_transform`, :code:`halfcomplex[]` and the
values :code:`complex[]` that would be returned if the same real input
sequence were passed to :func:`gsl_fft_complex_backward` as a complex
sequence (with imaginary parts set to :code:`0`)::

  complex[0].real  =  halfcomplex[0] 
  complex[0].imag  =  0
  complex[1].real  =  halfcomplex[1] 
  complex[1].imag  =  halfcomplex[2]
  complex[2].real  =  halfcomplex[3]
  complex[2].imag  =  halfcomplex[4]
  complex[3].real  =  halfcomplex[3]
  complex[3].imag  = -halfcomplex[4]
  complex[4].real  =  halfcomplex[1]
  complex[4].imag  = -halfcomplex[2]

The upper elements of the :code:`complex` array, :code:`complex[3]` and
:code:`complex[4]` are filled in using the symmetry condition.  The
imaginary part of the zero-frequency term :code:`complex[0].imag` is
known to be zero by the symmetry.

The next table shows the output for an even-length sequence, :math:`n=6`.
In the even case there are two values which are purely real::

  complex[0].real  =  halfcomplex[0]
  complex[0].imag  =  0
  complex[1].real  =  halfcomplex[1] 
  complex[1].imag  =  halfcomplex[2] 
  complex[2].real  =  halfcomplex[3] 
  complex[2].imag  =  halfcomplex[4] 
  complex[3].real  =  halfcomplex[5] 
  complex[3].imag  =  0 
  complex[4].real  =  halfcomplex[3] 
  complex[4].imag  = -halfcomplex[4]
  complex[5].real  =  halfcomplex[1] 
  complex[5].imag  = -halfcomplex[2] 

The upper elements of the :code:`complex` array, :code:`complex[4]` and
:code:`complex[5]` are filled in using the symmetry condition.  Both
:code:`complex[0].imag` and :code:`complex[3].imag` are known to be zero.

All these functions are declared in the header files
:file:`gsl_fft_real.h` and :file:`gsl_fft_halfcomplex.h`.

.. type:: gsl_fft_real_wavetable
          gsl_fft_halfcomplex_wavetable

   These data structures contain lookup tables for an FFT of a fixed size.

.. function:: gsl_fft_real_wavetable * gsl_fft_real_wavetable_alloc (size_t n)
              gsl_fft_halfcomplex_wavetable * gsl_fft_halfcomplex_wavetable_alloc (size_t n)

   These functions prepare trigonometric lookup tables for an FFT of size
   :math:`n` real elements.  The functions return a pointer to the newly
   allocated struct if no errors were detected, and a null pointer in the
   case of error.  The length :data:`n` is factorized into a product of
   subtransforms, and the factors and their trigonometric coefficients are
   stored in the wavetable. The trigonometric coefficients are computed
   using direct calls to :code:`sin` and :code:`cos`, for accuracy.
   Recursion relations could be used to compute the lookup table faster,
   but if an application performs many FFTs of the same length then
   computing the wavetable is a one-off overhead which does not affect the
   final throughput.

   The wavetable structure can be used repeatedly for any transform of the
   same length.  The table is not modified by calls to any of the other FFT
   functions.  The appropriate type of wavetable must be used for forward
   real or inverse half-complex transforms.

.. function:: void gsl_fft_real_wavetable_free (gsl_fft_real_wavetable * wavetable)
              void gsl_fft_halfcomplex_wavetable_free (gsl_fft_halfcomplex_wavetable * wavetable)

   These functions free the memory associated with the wavetable
   :data:`wavetable`. The wavetable can be freed if no further FFTs of the
   same length will be needed.

The mixed radix algorithms require additional working space to hold
the intermediate steps of the transform,

.. type:: gsl_fft_real_workspace

   This workspace contains parameters needed to compute a real FFT.

.. function:: gsl_fft_real_workspace * gsl_fft_real_workspace_alloc (size_t n)

   This function allocates a workspace for a real transform of length
   :data:`n`.  The same workspace can be used for both forward real and inverse
   halfcomplex transforms.

.. function:: void gsl_fft_real_workspace_free (gsl_fft_real_workspace * workspace)

   This function frees the memory associated with the workspace
   :data:`workspace`. The workspace can be freed if no further FFTs of the
   same length will be needed.

The following functions compute the transforms of real and half-complex
data,

.. function:: int gsl_fft_real_transform (double data[], size_t stride, size_t n, const gsl_fft_real_wavetable * wavetable, gsl_fft_real_workspace * work)
              int gsl_fft_halfcomplex_transform (double data[], size_t stride, size_t n, const gsl_fft_halfcomplex_wavetable * wavetable, gsl_fft_real_workspace * work)

   These functions compute the FFT of :data:`data`, a real or half-complex
   array of length :data:`n`, using a mixed radix decimation-in-frequency
   algorithm.  For :func:`gsl_fft_real_transform` :data:`data` is an array of
   time-ordered real data.  For :func:`gsl_fft_halfcomplex_transform`
   :data:`data` contains Fourier coefficients in the half-complex ordering
   described above.  There is no restriction on the length :data:`n`.
   Efficient modules are provided for subtransforms of length 2, 3, 4 and
   5.  Any remaining factors are computed with a slow, :math:`O(n^2)`,
   general-n module.  The caller must supply a :data:`wavetable` containing
   trigonometric lookup tables and a workspace :data:`work`. 

.. function:: int gsl_fft_real_unpack (const double real_coefficient[], gsl_complex_packed_array complex_coefficient, size_t stride, size_t n)

   This function converts a single real array, :data:`real_coefficient` into
   an equivalent complex array, :data:`complex_coefficient`, (with imaginary
   part set to zero), suitable for :code:`gsl_fft_complex` routines.  The
   algorithm for the conversion is simply::

      for (i = 0; i < n; i++)
        {
          complex_coefficient[i*stride].real = real_coefficient[i*stride];
          complex_coefficient[i*stride].imag = 0.0;
        }

.. function:: int gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[], gsl_complex_packed_array complex_coefficient, size_t stride, size_t n)

   This function converts :data:`halfcomplex_coefficient`, an array of
   half-complex coefficients as returned by :func:`gsl_fft_real_transform`, into an
   ordinary complex array, :data:`complex_coefficient`.  It fills in the
   complex array using the symmetry :math:`z_k = z_{n-k}^*`
   to reconstruct the redundant elements.  The algorithm for the conversion
   is::

      complex_coefficient[0].real = halfcomplex_coefficient[0];
      complex_coefficient[0].imag = 0.0;

      for (i = 1; i < n - i; i++)
        {
          double hc_real = halfcomplex_coefficient[(2 * i - 1)*stride];
          double hc_imag = halfcomplex_coefficient[(2 * i)*stride];
          complex_coefficient[i*stride].real = hc_real;
          complex_coefficient[i*stride].imag = hc_imag;
          complex_coefficient[(n - i)*stride].real = hc_real;
          complex_coefficient[(n - i)*stride].imag = -hc_imag;
        }

      if (i == n - i)
        {
          complex_coefficient[i*stride].real = halfcomplex_coefficient[(n - 1)*stride];
          complex_coefficient[i*stride].imag = 0.0;
        }

Here is an example program using :func:`gsl_fft_real_transform` and
:func:`gsl_fft_halfcomplex_inverse`.  It generates a real signal in the
shape of a square pulse.  The pulse is Fourier transformed to frequency
space, and all but the lowest ten frequency components are removed from
the array of Fourier coefficients returned by
:func:`gsl_fft_real_transform`.

The remaining Fourier coefficients are transformed back to the
time-domain, to give a filtered version of the square pulse.  Since
Fourier coefficients are stored using the half-complex symmetry both
positive and negative frequencies are removed and the final filtered
signal is also real.

.. include:: examples/fftreal.c
   :code:

The program output is shown in :numref:`fig_fft-real-mixedradix`.

.. _fig_fft-real-mixedradix:

.. figure:: /images/fft-real-mixedradix.png
   :scale: 100%

   Low-pass filtered version of a real pulse, output from the example program.

.. _fft-references:

References and Further Reading
==============================

A good starting point for learning more about the FFT is the following review
article,

* P. Duhamel and M. Vetterli.
  Fast Fourier transforms: A tutorial review and a state of the art.
  Signal Processing, 19:259--299, 1990.

To find out about the algorithms used in the GSL routines you may want
to consult the document "GSL FFT Algorithms" (it is included
in GSL, as :file:`doc/fftalgorithms.tex`).  This has general information
on FFTs and explicit derivations of the implementation for each
routine.  There are also references to the relevant literature.  For
convenience some of the more important references are reproduced below.

There are several introductory books on the FFT with example programs,
such as "The Fast Fourier Transform" by Brigham and "DFT/FFT
and Convolution Algorithms" by Burrus and Parks,

* E. Oran Brigham. "The Fast Fourier Transform".  Prentice Hall, 1974.

* C. S. Burrus and T. W. Parks.  "DFT/FFT and Convolution Algorithms",
  Wiley, 1984.

Both these introductory books cover the radix-2 FFT in some detail.
The mixed-radix algorithm at the heart of the |fftpack| routines is
reviewed in Clive Temperton's paper,

* Clive Temperton, Self-sorting mixed-radix fast Fourier transforms,
  Journal of Computational Physics, 52(1):1--23, 1983.

The derivation of FFTs for real-valued data is explained in the
following two articles,

* Henrik V. Sorenson, Douglas L. Jones, Michael T. Heideman, and C. Sidney
  Burrus.  Real-valued fast Fourier transform algorithms.
  "IEEE Transactions on Acoustics, Speech, and Signal Processing",
  ASSP-35(6):849--863, 1987.

* Clive Temperton.  Fast mixed-radix real Fourier transforms.
  "Journal of Computational Physics", 52:340--350, 1983.

In 1979 the IEEE published a compendium of carefully-reviewed Fortran
FFT programs in "Programs for Digital Signal Processing".  It is a
useful reference for implementations of many different FFT
algorithms,

* Digital Signal Processing Committee and IEEE Acoustics, Speech, and Signal
  Processing Committee, editors.
  Programs for Digital Signal Processing. IEEE Press, 1979.

For large-scale FFT work we recommend the use of the dedicated FFTW library
by Frigo and Johnson.  The FFTW library is self-optimizing---it
automatically tunes itself for each hardware platform in order to
achieve maximum performance.  It is available under the GNU GPL.

* FFTW Website, http://www.fftw.org/

The source code for |fftpack| is available from http://www.netlib.org/fftpack/
