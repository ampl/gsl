.. index:: multisets

*********
Multisets
*********

.. include:: include.rst

This chapter describes functions for creating and manipulating multisets. A
multiset :math:`c` is represented by an array of :math:`k` integers in the range
0 to :math:`n - 1`, where each value :math:`c_i` may occur more than once.  The
multiset :math:`c` corresponds to indices of :math:`k` elements chosen from an
:math:`n` element vector with replacement.  In mathematical terms, :math:`n` is
the cardinality of the multiset while :math:`k` is the maximum multiplicity of
any value.  Multisets are useful, for example, when iterating over the indices
of a :math:`k`-th order symmetric tensor in :math:`n`-space.

The functions described in this chapter are defined in the header file
:file:`gsl_multiset.h`.

The Multiset struct
===================

.. type:: gsl_multiset

   A multiset is defined by a structure containing three components, the
   values of :math:`n` and :math:`k`, and a pointer to the multiset array.
   The elements of the multiset array are all of type :code:`size_t`, and
   are stored in increasing order.  The :type:`gsl_multiset` structure
   looks like this::

      typedef struct
      {
        size_t n;
        size_t k;
        size_t *data;
      } gsl_multiset;

Multiset allocation
===================

.. function:: gsl_multiset * gsl_multiset_alloc (size_t n, size_t k)

   This function allocates memory for a new multiset with parameters :data:`n`,
   :data:`k`.  The multiset is not initialized and its elements are undefined.  Use
   the function :func:`gsl_multiset_calloc` if you want to create a multiset which
   is initialized to the lexicographically first multiset element. A null pointer
   is returned if insufficient memory is available to create the multiset.

.. function:: gsl_multiset * gsl_multiset_calloc (size_t n, size_t k)

   This function allocates memory for a new multiset with parameters :data:`n`,
   :data:`k` and initializes it to the lexicographically first multiset element. A
   null pointer is returned if insufficient memory is available to create the
   multiset.

.. function:: void gsl_multiset_init_first (gsl_multiset * c)

   This function initializes the multiset :data:`c` to the lexicographically first
   multiset element, i.e. :math:`0` repeated :math:`k` times.

.. function:: void gsl_multiset_init_last (gsl_multiset * c)

   This function initializes the multiset :data:`c` to the lexicographically last
   multiset element, i.e. :math:`n-1` repeated :math:`k` times.

.. function:: void gsl_multiset_free (gsl_multiset * c)

   This function frees all the memory used by the multiset :data:`c`.

.. function:: int gsl_multiset_memcpy (gsl_multiset * dest, const gsl_multiset * src)

   This function copies the elements of the multiset :data:`src` into the
   multiset :data:`dest`.  The two multisets must have the same size.

Accessing multiset elements
===========================

The following function can be used to access the elements of a multiset.

.. function:: size_t gsl_multiset_get (const gsl_multiset * c, const size_t i)

   This function returns the value of the :data:`i`-th element of the
   multiset :data:`c`.  If :data:`i` lies outside the allowed range of 0 to
   :math:`k - 1` then the error handler is invoked and 0 is returned. |inlinefn|

Multiset properties
===================

.. function:: size_t gsl_multiset_n (const gsl_multiset * c)

   This function returns the range (:math:`n`) of the multiset :data:`c`.

.. function:: size_t gsl_multiset_k (const gsl_multiset * c)

   This function returns the number of elements (:math:`k`) in the multiset :data:`c`.

.. function:: size_t * gsl_multiset_data (const gsl_multiset * c)

   This function returns a pointer to the array of elements in the
   multiset :data:`c`.

.. index::
   single: checking multiset for validity
   single: testing multiset for validity

.. function:: int gsl_multiset_valid (gsl_multiset * c)

   This function checks that the multiset :data:`c` is valid.  The :data:`k`
   elements should lie in the range 0 to :math:`n - 1`, with each
   value occurring in nondecreasing order.

Multiset functions
==================

.. index:: iterating through multisets

.. function:: int gsl_multiset_next (gsl_multiset * c)

   This function advances the multiset :data:`c` to the next multiset element in
   lexicographic order and returns :macro:`GSL_SUCCESS`.  If no further multisets
   elements are available it returns :macro:`GSL_FAILURE` and leaves :data:`c`
   unmodified.  Starting with the first multiset and repeatedly applying this
   function will iterate through all possible multisets of a given order.

.. function:: int gsl_multiset_prev (gsl_multiset * c)

   This function steps backwards from the multiset :data:`c` to the previous
   multiset element in lexicographic order, returning :macro:`GSL_SUCCESS`.  If no
   previous multiset is available it returns :macro:`GSL_FAILURE` and leaves :data:`c`
   unmodified.

Reading and writing multisets
=============================

The library provides functions for reading and writing multisets to a
file as binary data or formatted text.

.. function:: int gsl_multiset_fwrite (FILE * stream, const gsl_multiset * c)

   This function writes the elements of the multiset :data:`c` to the
   stream :data:`stream` in binary format.  The function returns
   :macro:`GSL_EFAILED` if there was a problem writing to the file.  Since the
   data is written in the native binary format it may not be portable
   between different architectures.

.. function:: int gsl_multiset_fread (FILE * stream, gsl_multiset * c)

   This function reads elements from the open stream :data:`stream` into the
   multiset :data:`c` in binary format.  The multiset :data:`c` must be
   preallocated with correct values of :math:`n` and :math:`k` since the
   function uses the size of :data:`c` to determine how many bytes to read.
   The function returns :macro:`GSL_EFAILED` if there was a problem reading
   from the file.  The data is assumed to have been written in the native
   binary format on the same architecture.

.. function:: int gsl_multiset_fprintf (FILE * stream, const gsl_multiset * c, const char * format)

   This function writes the elements of the multiset :data:`c`
   line-by-line to the stream :data:`stream` using the format specifier
   :data:`format`, which should be suitable for a type of :code:`size_t`.
   In ISO C99 the type modifier :code:`z` represents :code:`size_t`, so
   :code:`"%zu\n"` is a suitable format [#f1]_.
   The function returns :macro:`GSL_EFAILED` if there was a problem writing to the file.

.. function:: int gsl_multiset_fscanf (FILE * stream, gsl_multiset * c)

   This function reads formatted data from the stream :data:`stream` into the
   multiset :data:`c`.  The multiset :data:`c` must be preallocated with
   correct values of :math:`n` and :math:`k` since the function uses the size of :data:`c` to
   determine how many numbers to read.  The function returns
   :macro:`GSL_EFAILED` if there was a problem reading from the file.

Examples
========

The example program below prints all multisets elements containing the values
:math:`{0,1,2,3}` ordered by size.  Multiset elements of the same size are
ordered lexicographically.

.. include:: examples/multiset.c
   :code:

Here is the output from the program,

.. include:: examples/multiset.txt
   :code:

All 70 multisets are generated and sorted lexicographically.

.. rubric:: Footnotes

.. [#f1] In versions of the GNU C library prior to the ISO C99 standard,
         the type modifier :code:`Z` was used instead.
