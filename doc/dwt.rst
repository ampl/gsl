.. index::
   single: DWT, see wavelet transforms
   single: wavelet transforms
   single: transforms, wavelet

******************
Wavelet Transforms
******************

This chapter describes functions for performing Discrete Wavelet
Transforms (DWTs).  The library includes wavelets for real data in both
one and two dimensions.  The wavelet functions are declared in the header
files :file:`gsl_wavelet.h` and :file:`gsl_wavelet2d.h`.

.. index::
   single: DWT, mathematical definition

Definitions
===========

The continuous wavelet transform and its inverse are defined by
the relations,

.. math:: w(s, \tau) = \int_{-\infty}^\infty f(t) * \psi^*_{s,\tau}(t) dt

and,

.. math:: f(t) = \int_0^\infty ds \int_{-\infty}^\infty w(s, \tau) * \psi_{s,\tau}(t) d\tau

where the basis functions :math:`\psi_{s,\tau}`
are obtained by scaling
and translation from a single function, referred to as the *mother wavelet*.

The discrete version of the wavelet transform acts on equally-spaced
samples, with fixed scaling and translation steps (:math:`s`,
:math:`\tau`).  The frequency and time axes are sampled *dyadically*
on scales of :math:`2^j` through a level parameter :math:`j`.

..  The wavelet :math:`\psi`
..  can be expressed in terms of a scaling function :math:`\varphi`,
..
..  @tex
..  \beforedisplay
..  $$
..  \psi(2^{j-1},t) = \sum_{k=0}^{2^j-1} g_j(k) * \bar{\varphi}(2^j t-k)
..  $$
..  \afterdisplay
..  @end tex
..  @ifinfo
..  @example
..  \psi(2^@{j-1@},t) = \sum_@{k=0@}^@{2^j-1@} g_j(k) * \bar@{\varphi@}(2^j t-k)
..  @end example
..  @end ifinfo
..  @noindent
..  and
..
..  @tex
..  \beforedisplay
..  $$
..  \varphi(2^{j-1},t) = \sum_{k=0}^{2^j-1} h_j(k) * \bar{\varphi}(2^j t-k)
..  $$
..  \afterdisplay
..  @end tex
..  @ifinfo
..  @example
..  \varphi(2^@{j-1@},t) = \sum_@{k=0@}^@{2^j-1@} h_j(k) * \bar@{\varphi@}(2^j t-k)
..  @end example
..  @end ifinfo
..  @noindent
..  The functions :math:`\psi` and :math:`\varphi` are related through the
..  coefficients
..  @c{$g_{n} = (-1)^n h_{L-1-n}$}
..  @math{g_@{n@} = (-1)^n h_@{L-1-n@}}
..  for @c{$n=0 \dots L-1$}
..  @math{n=0 ... L-1},
..  where :math:`L` is the total number of coefficients.  The two sets of
..  coefficients :math:`h_j` and :math:`g_i` define the scaling function and
.. the wavelet.  

The resulting family of functions :math:`\{\psi_{j,n}\}`
constitutes an orthonormal basis for square-integrable signals.  
The discrete wavelet transform is an :math:`O(N)` algorithm, and is also
referred to as the *fast wavelet transform*.

.. index:: DWT initialization

Initialization
==============

.. type:: gsl_wavelet

   This structure contains the filter coefficients
   defining the wavelet and any associated offset parameters.

.. function:: gsl_wavelet * gsl_wavelet_alloc (const gsl_wavelet_type * T, size_t k)

   This function allocates and initializes a wavelet object of type
   :data:`T`.  The parameter :data:`k` selects the specific member of the
   wavelet family.  A null pointer is returned if insufficient memory is
   available or if a unsupported member is selected.

The following wavelet types are implemented:

.. type:: gsl_wavelet_type

   .. index::
      single: Daubechies wavelets
      single: maximal phase, Daubechies wavelets

   .. var:: gsl_wavelet_type * gsl_wavelet_daubechies
            gsl_wavelet_type * gsl_wavelet_daubechies_centered

      This is the Daubechies wavelet family of maximum phase with :math:`k/2`
      vanishing moments.  The implemented wavelets are 
      :math:`k=4, 6, \dots, 20`, with :data:`k` even.

   .. index:: Haar wavelets

   .. var:: gsl_wavelet_type * gsl_wavelet_haar
            gsl_wavelet_type * gsl_wavelet_haar_centered

      This is the Haar wavelet.  The only valid choice of :math:`k` for the
      Haar wavelet is :math:`k=2`.

   .. index::
      single: biorthogonal wavelets
      single: B-spline wavelets

   .. var:: gsl_wavelet_type * gsl_wavelet_bspline
            gsl_wavelet_type * gsl_wavelet_bspline_centered

      This is the biorthogonal B-spline wavelet family of order :math:`(i,j)`.  
      The implemented values of :math:`k = 100*i + j` are 103, 105, 202, 204,
      206, 208, 301, 303, 305 307, 309.

The centered forms of the wavelets align the coefficients of the various
sub-bands on edges.  Thus the resulting visualization of the
coefficients of the wavelet transform in the phase plane is easier to
understand.

.. function:: const char * gsl_wavelet_name (const gsl_wavelet * w)

   This function returns a pointer to the name of the wavelet family for
   :data:`w`.

..  @deftypefun {void} gsl_wavelet_print (const gsl_wavelet * w)
..  This function prints the filter coefficients (@code{**h1}, @code{**g1}, @code{**h2}, @code{**g2}) of the wavelet object :data:`w`.
..  @end deftypefun

.. function:: void gsl_wavelet_free (gsl_wavelet * w)

   This function frees the wavelet object :data:`w`.

.. type:: gsl_wavelet_workspace

   This structure contains scratch space of the
   same size as the input data and is used to hold intermediate results
   during the transform.

.. function:: gsl_wavelet_workspace * gsl_wavelet_workspace_alloc (size_t n)

   This function allocates a workspace for the discrete wavelet transform.
   To perform a one-dimensional transform on :data:`n` elements, a workspace
   of size :data:`n` must be provided.  For two-dimensional transforms of
   :data:`n`-by-:data:`n` matrices it is sufficient to allocate a workspace of
   size :data:`n`, since the transform operates on individual rows and
   columns. A null pointer is returned if insufficient memory is available.

.. function:: void gsl_wavelet_workspace_free (gsl_wavelet_workspace * work)

   This function frees the allocated workspace :data:`work`.

Transform Functions
===================

This sections describes the actual functions performing the discrete
wavelet transform.  Note that the transforms use periodic boundary
conditions.  If the signal is not periodic in the sample length then
spurious coefficients will appear at the beginning and end of each level
of the transform.

.. index::
   single: DWT, one dimensional

Wavelet transforms in one dimension
-----------------------------------

.. function:: int gsl_wavelet_transform (const gsl_wavelet * w, double * data, size_t stride, size_t n, gsl_wavelet_direction dir, gsl_wavelet_workspace * work)
              int gsl_wavelet_transform_forward (const gsl_wavelet * w, double * data, size_t stride, size_t n, gsl_wavelet_workspace * work)
              int gsl_wavelet_transform_inverse (const gsl_wavelet * w, double * data, size_t stride, size_t n, gsl_wavelet_workspace * work)

   These functions compute in-place forward and inverse discrete wavelet
   transforms of length :data:`n` with stride :data:`stride` on the array
   :data:`data`. The length of the transform :data:`n` is restricted to powers
   of two.  For the :code:`transform` version of the function the argument
   :data:`dir` can be either :code:`forward` (:math:`+1`) or :code:`backward`
   (:math:`-1`).  A workspace :data:`work` of length :data:`n` must be provided.

   For the forward transform, the elements of the original array are 
   replaced by the discrete wavelet
   transform :math:`f_i \rightarrow w_{j,k}`
   in a packed triangular storage layout, 
   where :data:`j` is the index of the level 
   :math:`j = 0 \dots J-1`
   and
   :data:`k` is the index of the coefficient within each level,
   :math:`k = 0 \dots 2^j - 1`.
   The total number of levels is :math:`J = \log_2(n)`.  The output data
   has the following form,

   .. math:: (s_{-1,0}, d_{0,0}, d_{1,0}, d_{1,1}, d_{2,0},\cdots, d_{j,k},\cdots, d_{J-1,2^{J-1} - 1}) 

   where the first element is the smoothing coefficient :math:`s_{-1,0}`,
   followed by the detail coefficients :math:`d_{j,k}`
   for each level
   :math:`j`.  The backward transform inverts these coefficients to obtain 
   the original data.

   These functions return a status of :macro:`GSL_SUCCESS` upon successful
   completion.  :macro:`GSL_EINVAL` is returned if :data:`n` is not an integer
   power of 2 or if insufficient workspace is provided.

.. index::
   single: DWT, two dimensional

Wavelet transforms in two dimension
-----------------------------------

The library provides functions to perform two-dimensional discrete
wavelet transforms on square matrices.  The matrix dimensions must be an
integer power of two.  There are two possible orderings of the rows and
columns in the two-dimensional wavelet transform, referred to as the
"standard" and "non-standard" forms.

The "standard" transform performs a complete discrete wavelet
transform on the rows of the matrix, followed by a separate complete
discrete wavelet transform on the columns of the resulting
row-transformed matrix.  This procedure uses the same ordering as a
two-dimensional Fourier transform.

The "non-standard" transform is performed in interleaved passes on the
rows and columns of the matrix for each level of the transform.  The
first level of the transform is applied to the matrix rows, and then to
the matrix columns.  This procedure is then repeated across the rows and
columns of the data for the subsequent levels of the transform, until
the full discrete wavelet transform is complete.  The non-standard form
of the discrete wavelet transform is typically used in image analysis.

The functions described in this section are declared in the header file
:file:`gsl_wavelet2d.h`.

.. function:: int gsl_wavelet2d_transform (const gsl_wavelet * w, double * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_direction dir, gsl_wavelet_workspace * work)
              int gsl_wavelet2d_transform_forward (const gsl_wavelet * w, double * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work)
              int gsl_wavelet2d_transform_inverse (const gsl_wavelet * w, double * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work)

   These functions compute two-dimensional in-place forward and inverse
   discrete wavelet transforms in standard form on the
   array :data:`data` stored in row-major form with dimensions :data:`size1`
   and :data:`size2` and physical row length :data:`tda`.  The dimensions must
   be equal (square matrix) and are restricted to powers of two.  For the
   :code:`transform` version of the function the argument :data:`dir` can be
   either :code:`forward` (:math:`+1`) or :code:`backward` (:math:`-1`).  A
   workspace :data:`work` of the appropriate size must be provided.  On exit,
   the appropriate elements of the array :data:`data` are replaced by their
   two-dimensional wavelet transform.

   The functions return a status of :macro:`GSL_SUCCESS` upon successful
   completion.  :macro:`GSL_EINVAL` is returned if :data:`size1` and
   :data:`size2` are not equal and integer powers of 2, or if insufficient
   workspace is provided.

.. function:: int gsl_wavelet2d_transform_matrix (const gsl_wavelet * w, gsl_matrix * m, gsl_wavelet_direction dir, gsl_wavelet_workspace * work)
              int gsl_wavelet2d_transform_matrix_forward (const gsl_wavelet * w, gsl_matrix * m, gsl_wavelet_workspace * work)
              int gsl_wavelet2d_transform_matrix_inverse (const gsl_wavelet * w, gsl_matrix * m, gsl_wavelet_workspace * work)

   These functions compute the two-dimensional in-place wavelet transform
   on a matrix :data:`m`.

.. function:: int gsl_wavelet2d_nstransform (const gsl_wavelet * w, double * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_direction dir, gsl_wavelet_workspace * work)
              int gsl_wavelet2d_nstransform_forward (const gsl_wavelet * w, double * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work)
              int gsl_wavelet2d_nstransform_inverse (const gsl_wavelet * w, double * data, size_t tda, size_t size1, size_t size2, gsl_wavelet_workspace * work)

   These functions compute the two-dimensional wavelet transform in
   non-standard form.

.. function:: int gsl_wavelet2d_nstransform_matrix (const gsl_wavelet * w, gsl_matrix * m, gsl_wavelet_direction dir, gsl_wavelet_workspace * work)
              int gsl_wavelet2d_nstransform_matrix_forward (const gsl_wavelet * w, gsl_matrix * m, gsl_wavelet_workspace * work)
              int gsl_wavelet2d_nstransform_matrix_inverse (const gsl_wavelet * w, gsl_matrix * m, gsl_wavelet_workspace * work)

   These functions compute the non-standard form of the two-dimensional
   in-place wavelet transform on a matrix :data:`m`.

Examples
========

The following program demonstrates the use of the one-dimensional
wavelet transform functions.  It computes an approximation to an input
signal (of length 256) using the 20 largest components of the wavelet
transform, while setting the others to zero.

.. include:: examples/dwt.c
   :code:

The output can be used with the GNU plotutils :code:`graph` program::

  $ ./a.out ecg.dat > dwt.txt
  $ graph -T ps -x 0 256 32 -h 0.3 -a dwt.txt > dwt.ps

:numref:`fig_dwt` shows an original and compressed version of a sample ECG
recording from the MIT-BIH Arrhythmia Database, part of the PhysioNet
archive of public-domain of medical datasets.

.. _fig_dwt:

.. figure:: /images/dwt.png
   :scale: 60%

   Original (upper) and wavelet-compressed (lower) ECG signals, using the
   20 largest components of the Daubechies(4) discrete wavelet transform.

References and Further Reading
==============================

The mathematical background to wavelet transforms is covered in the
original lectures by Daubechies,

* Ingrid Daubechies.
  Ten Lectures on Wavelets.
  *CBMS-NSF Regional Conference Series in Applied Mathematics* (1992), 
  SIAM, ISBN 0898712742.

An easy to read introduction to the subject with an emphasis on the
application of the wavelet transform in various branches of science is,

* Paul S. Addison. *The Illustrated Wavelet Transform Handbook*.
  Institute of Physics Publishing (2002), ISBN 0750306920.

For extensive coverage of signal analysis by wavelets, wavelet packets
and local cosine bases see,

* S. G. Mallat. *A wavelet tour of signal processing* (Second
  edition). Academic Press (1999), ISBN 012466606X.

The concept of multiresolution analysis underlying the wavelet transform
is described in,

* S. G. Mallat.
  Multiresolution Approximations and Wavelet Orthonormal Bases of L^2(R).
  *Transactions of the American Mathematical Society*, 315(1), 1989, 69--87.

* S. G. Mallat.
  A Theory for Multiresolution Signal Decomposition---The Wavelet Representation.
  *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 11, 1989,
  674--693. 

The coefficients for the individual wavelet families implemented by the
library can be found in the following papers,

* I. Daubechies.
  Orthonormal Bases of Compactly Supported Wavelets.
  *Communications on Pure and Applied Mathematics*, 41 (1988) 909--996.

* A. Cohen, I. Daubechies, and J.-C. Feauveau.
  Biorthogonal Bases of Compactly Supported Wavelets.
  *Communications on Pure and Applied Mathematics*, 45 (1992)
  485--560.

The PhysioNet archive of physiological datasets can be found online at
http://www.physionet.org/ and is described in the following
paper,

* Goldberger et al.  
  PhysioBank, PhysioToolkit, and PhysioNet: Components
  of a New Research Resource for Complex Physiologic
  Signals. 
  *Circulation* 101(23):e215-e220 2000.
