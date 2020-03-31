.. index:: permutations

************
Permutations
************

.. include:: include.rst

This chapter describes functions for creating and manipulating
permutations. A permutation :math:`p` is represented by an array of
:math:`n` integers in the range 0 to :math:`n-1`, where each value
:math:`p_i` occurs once and only once.  The application of a permutation
:math:`p` to a vector :math:`v` yields a new vector :math:`v'` where
:math:`v'_i = v_{p_i}`.
For example, the array :math:`(0,1,3,2)` represents a permutation
which exchanges the last two elements of a four element vector.
The corresponding identity permutation is :math:`(0,1,2,3)`.   

Note that the permutations produced by the linear algebra routines
correspond to the exchange of matrix columns, and so should be considered
as applying to row-vectors in the form :math:`v' = v P` rather than
column-vectors, when permuting the elements of a vector.

The functions described in this chapter are defined in the header file
:file:`gsl_permutation.h`.

The Permutation struct
======================

.. type:: gsl_permutation

   A permutation is defined by a structure containing two components, the size
   of the permutation and a pointer to the permutation array.  The elements
   of the permutation array are all of type :code:`size_t`.  The
   :type:`gsl_permutation` structure looks like this::

      typedef struct
      {
        size_t size;
        size_t * data;
      } gsl_permutation;

Permutation allocation
======================

.. function:: gsl_permutation * gsl_permutation_alloc (size_t n)

   This function allocates memory for a new permutation of size :data:`n`.
   The permutation is not initialized and its elements are undefined.  Use
   the function :func:`gsl_permutation_calloc` if you want to create a
   permutation which is initialized to the identity. A null pointer is
   returned if insufficient memory is available to create the permutation.

.. function:: gsl_permutation * gsl_permutation_calloc (size_t n)

   This function allocates memory for a new permutation of size :data:`n` and
   initializes it to the identity. A null pointer is returned if
   insufficient memory is available to create the permutation.

.. index:: identity permutation

.. function:: void gsl_permutation_init (gsl_permutation * p)

   This function initializes the permutation :data:`p` to the identity, i.e.
   :math:`(0, 1, 2, \dots, n - 1)`.

.. function:: void gsl_permutation_free (gsl_permutation * p)

   This function frees all the memory used by the permutation :data:`p`.

.. function:: int gsl_permutation_memcpy (gsl_permutation * dest, const gsl_permutation * src)

   This function copies the elements of the permutation :data:`src` into the
   permutation :data:`dest`.  The two permutations must have the same size.

Accessing permutation elements
==============================

The following functions can be used to access and manipulate
permutations.

.. function:: size_t gsl_permutation_get (const gsl_permutation * p, const size_t i)

   This function returns the value of the :data:`i`-th element of the
   permutation :data:`p`.  If :data:`i` lies outside the allowed range of 0 to
   :math:`n - 1` then the error handler is invoked and 0 is returned. |inlinefn|

.. index::
   single: exchanging permutation elements
   single: swapping permutation elements

.. function:: int gsl_permutation_swap (gsl_permutation * p, const size_t i, const size_t j)

   This function exchanges the :data:`i`-th and :data:`j`-th elements of the
   permutation :data:`p`.

Permutation properties
======================

.. function:: size_t gsl_permutation_size (const gsl_permutation * p)

   This function returns the size of the permutation :data:`p`.

.. function:: size_t * gsl_permutation_data (const gsl_permutation * p)

   This function returns a pointer to the array of elements in the
   permutation :data:`p`.

.. index::
   single: checking permutation for validity
   single: testing permutation for validity

.. function:: int gsl_permutation_valid (const gsl_permutation * p)

   This function checks that the permutation :data:`p` is valid.  The :code:`n`
   elements should contain each of the numbers 0 to :code:`n - 1` once and only
   once.

Permutation functions
=====================

.. index:: reversing a permutation

.. function:: void gsl_permutation_reverse (gsl_permutation * p)

   This function reverses the elements of the permutation :data:`p`.

.. index:: inverting a permutation

.. function:: int gsl_permutation_inverse (gsl_permutation * inv, const gsl_permutation * p)

   This function computes the inverse of the permutation :data:`p`, storing
   the result in :data:`inv`.

.. index:: iterating through permutations

.. function:: int gsl_permutation_next (gsl_permutation * p)

   This function advances the permutation :data:`p` to the next permutation
   in lexicographic order and returns :macro:`GSL_SUCCESS`.  If no further
   permutations are available it returns :macro:`GSL_FAILURE` and leaves
   :data:`p` unmodified.  Starting with the identity permutation and
   repeatedly applying this function will iterate through all possible
   permutations of a given order.

.. function:: int gsl_permutation_prev (gsl_permutation * p)

   This function steps backwards from the permutation :data:`p` to the
   previous permutation in lexicographic order, returning
   :macro:`GSL_SUCCESS`.  If no previous permutation is available it returns
   :macro:`GSL_FAILURE` and leaves :data:`p` unmodified.

Applying Permutations
=====================

The following functions are defined in the header files :file:`gsl_permute.h`
and :file:`gsl_permute_vector.h`.

.. function:: int gsl_permute (const size_t * p, double * data, size_t stride, size_t n)

   This function applies the permutation :data:`p` to the array :data:`data` of
   size :data:`n` with stride :data:`stride`.

.. function:: int gsl_permute_inverse (const size_t * p, double * data, size_t stride, size_t n)

   This function applies the inverse of the permutation :data:`p` to the
   array :data:`data` of size :data:`n` with stride :data:`stride`.

.. function:: int gsl_permute_vector (const gsl_permutation * p, gsl_vector * v)

   This function applies the permutation :data:`p` to the elements of the
   vector :data:`v`, considered as a row-vector acted on by a permutation
   matrix from the right, :math:`v' = v P`.  The :math:`j`-th column of the
   permutation matrix :math:`P` is given by the :math:`p_j`-th column of the
   identity matrix. The permutation :data:`p` and the vector :data:`v` must
   have the same length.

.. function:: int gsl_permute_vector_inverse (const gsl_permutation * p, gsl_vector * v)

   This function applies the inverse of the permutation :data:`p` to the
   elements of the vector :data:`v`, considered as a row-vector acted on by
   an inverse permutation matrix from the right, :math:`v' = v P^T`.  Note
   that for permutation matrices the inverse is the same as the transpose.
   The :math:`j`-th column of the permutation matrix :math:`P` is given by
   the :math:`p_j`-th column of the identity matrix. The permutation :data:`p`
   and the vector :data:`v` must have the same length.

.. function:: int gsl_permute_matrix (const gsl_permutation * p, gsl_matrix * A)

   This function applies the permutation :data:`p` to the matrix :data:`A` from
   the right, :math:`A' = A P`.  The :math:`j`-th column of the
   permutation matrix :math:`P` is given by the :math:`p_j`-th column of the
   identity matrix. This effectively permutes the columns of :data:`A` according
   to the permutation :data:`p`, and so the number of columns of :data:`A` must
   equal the size of the permutation :data:`p`.

.. function:: int gsl_permutation_mul (gsl_permutation * p, const gsl_permutation * pa, const gsl_permutation * pb)

   This function combines the two permutations :data:`pa` and :data:`pb` into a
   single permutation :data:`p`, where :math:`p = pa * pb`
   The permutation :data:`p` is equivalent to applying :data:`pb` first and
   then :data:`pa`.

Reading and writing permutations
================================

The library provides functions for reading and writing permutations to a
file as binary data or formatted text.

.. function:: int gsl_permutation_fwrite (FILE * stream, const gsl_permutation * p)

   This function writes the elements of the permutation :data:`p` to the
   stream :data:`stream` in binary format.  The function returns
   :macro:`GSL_EFAILED` if there was a problem writing to the file.  Since the
   data is written in the native binary format it may not be portable
   between different architectures.

.. function:: int gsl_permutation_fread (FILE * stream, gsl_permutation * p)

   This function reads into the permutation :data:`p` from the open stream
   :data:`stream` in binary format.  The permutation :data:`p` must be
   preallocated with the correct length since the function uses the size of
   :data:`p` to determine how many bytes to read.  The function returns
   :macro:`GSL_EFAILED` if there was a problem reading from the file.  The
   data is assumed to have been written in the native binary format on the
   same architecture.

.. function:: int gsl_permutation_fprintf (FILE * stream, const gsl_permutation * p, const char * format)

   This function writes the elements of the permutation :data:`p`
   line-by-line to the stream :data:`stream` using the format specifier
   :data:`format`, which should be suitable for a type of :data:`size_t`. 
   In ISO C99 the type modifier :code:`z` represents :code:`size_t`, so
   :code:`"%zu\n"` is a suitable format [#f1]_.
   The function returns :macro:`GSL_EFAILED` if there was a problem writing
   to the file.

.. function:: int gsl_permutation_fscanf (FILE * stream, gsl_permutation * p)

   This function reads formatted data from the stream :data:`stream` into the
   permutation :data:`p`.  The permutation :data:`p` must be preallocated with
   the correct length since the function uses the size of :data:`p` to
   determine how many numbers to read.  The function returns
   :macro:`GSL_EFAILED` if there was a problem reading from the file.

Permutations in cyclic form
===========================

A permutation can be represented in both *linear* and *cyclic*
notations.  The functions described in this section convert between the
two forms.  The linear notation is an index mapping, and has already
been described above.  The cyclic notation expresses a permutation as a
series of circular rearrangements of groups of elements, or
*cycles*.

For example, under the cycle (1 2 3), 1 is replaced by 2, 2 is replaced
by 3 and 3 is replaced by 1 in a circular fashion. Cycles of different
sets of elements can be combined independently, for example (1 2 3) (4
5) combines the cycle (1 2 3) with the cycle (4 5), which is an exchange
of elements 4 and 5.  A cycle of length one represents an element which
is unchanged by the permutation and is referred to as a *singleton*.

It can be shown that every permutation can be decomposed into
combinations of cycles.  The decomposition is not unique, but can always
be rearranged into a standard *canonical form* by a reordering of
elements.  The library uses the canonical form defined in Knuth's
*Art of Computer Programming* (Vol 1, 3rd Ed, 1997) Section 1.3.3,
p.178.

The procedure for obtaining the canonical form given by Knuth is,

#. Write all singleton cycles explicitly
#. Within each cycle, put the smallest number first
#. Order the cycles in decreasing order of the first number in the cycle.

For example, the linear representation (2 4 3 0 1) is represented as (1
4) (0 2 3) in canonical form. The permutation corresponds to an
exchange of elements 1 and 4, and rotation of elements 0, 2 and 3.

The important property of the canonical form is that it can be
reconstructed from the contents of each cycle without the brackets. In
addition, by removing the brackets it can be considered as a linear
representation of a different permutation. In the example given above
the permutation (2 4 3 0 1) would become (1 4 0 2 3).  This mapping has
many applications in the theory of permutations.

.. function:: int gsl_permutation_linear_to_canonical (gsl_permutation * q, const gsl_permutation * p)

   This function computes the canonical form of the permutation :data:`p` and
   stores it in the output argument :data:`q`.

.. function:: int gsl_permutation_canonical_to_linear (gsl_permutation * p, const gsl_permutation * q)

   This function converts a permutation :data:`q` in canonical form back into
   linear form storing it in the output argument :data:`p`.

.. function:: size_t gsl_permutation_inversions (const gsl_permutation * p)

   This function counts the number of inversions in the permutation
   :data:`p`.  An inversion is any pair of elements that are not in order.
   For example, the permutation 2031 has three inversions, corresponding to
   the pairs (2,0) (2,1) and (3,1).  The identity permutation has no
   inversions.

.. function:: size_t gsl_permutation_linear_cycles (const gsl_permutation * p)

   This function counts the number of cycles in the permutation :data:`p`, given in linear form.

.. function:: size_t gsl_permutation_canonical_cycles (const gsl_permutation * q)

   This function counts the number of cycles in the permutation :data:`q`, given in canonical form.

Examples
========

The example program below creates a random permutation (by shuffling the
elements of the identity) and finds its inverse.

.. include:: examples/permshuffle.c
   :code:

Here is the output from the program::

   $ ./a.out 
   initial permutation: 0 1 2 3 4 5 6 7 8 9
    random permutation: 1 3 5 2 7 6 0 4 9 8
   inverse permutation: 6 0 3 1 7 2 5 4 9 8

The random permutation :code:`p[i]` and its inverse :code:`q[i]` are
related through the identity :code:`p[q[i]] = i`, which can be verified
from the output.

The next example program steps forwards through all possible third order
permutations, starting from the identity,

.. include:: examples/permseq.c
   :code:

Here is the output from the program::

   $ ./a.out 
    0 1 2
    0 2 1
    1 0 2
    1 2 0
    2 0 1
    2 1 0

The permutations are generated in lexicographic order.  To reverse the
sequence, begin with the final permutation (which is the reverse of the
identity) and replace :func:`gsl_permutation_next` with
:func:`gsl_permutation_prev`.

References and Further Reading
==============================

The subject of permutations is covered extensively in the following,

* Donald E. Knuth, The Art of Computer Programming: Sorting and
  Searching (Vol 3, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850.

For the definition of the *canonical form* see,

* Donald E. Knuth, The Art of Computer Programming: Fundamental
  Algorithms (Vol 1, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850.
  Section 1.3.3, An Unusual Correspondence, p.178--179.

.. rubric:: Footnotes

.. [#f1] In versions of the GNU C library prior to the ISO C99 standard, 
         the type modifier :code:`Z` was used instead.
