.. index:: combinations

************
Combinations
************

.. include:: include.rst

This chapter describes functions for creating and manipulating
combinations. A combination :math:`c` is represented by an array of
:math:`k` integers in the range 0 to :math:`n - 1`, where each value
:math:`c_i` occurs at most once.  The combination :math:`c` corresponds to
indices of :math:`k` elements chosen from an :math:`n` element vector.
Combinations are useful for iterating over all :math:`k`-element subsets
of a set.

The functions described in this chapter are defined in the header file
:file:`gsl_combination.h`.

The Combination struct
======================

.. type:: gsl_combination

   A combination is defined by a structure containing three components, the
   values of :math:`n` and :math:`k`, and a pointer to the combination array.
   The elements of the combination array are all of type :code:`size_t`, and
   are stored in increasing order.  The :type:`gsl_combination` structure
   looks like this::

      typedef struct
      {
        size_t n;
        size_t k;
        size_t *data;
      } gsl_combination;

Combination allocation
======================

.. function:: gsl_combination * gsl_combination_alloc (size_t n, size_t k)

   This function allocates memory for a new combination with parameters
   :data:`n`, :data:`k`.  The combination is not initialized and its elements
   are undefined.  Use the function :func:`gsl_combination_calloc` if you
   want to create a combination which is initialized to the
   lexicographically first combination. A null pointer is returned if
   insufficient memory is available to create the combination.

.. function:: gsl_combination * gsl_combination_calloc (size_t n, size_t k)

   This function allocates memory for a new combination with parameters
   :data:`n`, :data:`k` and initializes it to the lexicographically first
   combination. A null pointer is returned if insufficient memory is
   available to create the combination.

.. function:: void gsl_combination_init_first (gsl_combination * c)

   This function initializes the combination :data:`c` to the
   lexicographically first combination, i.e.  :math:`(0, 1, 2, \dots, k - 1)`.

.. function:: void gsl_combination_init_last (gsl_combination * c)

   This function initializes the combination :data:`c` to the
   lexicographically last combination, i.e.  :math:`(n - k, n - k + 1, \dots, n - 1)`.

.. function:: void gsl_combination_free (gsl_combination * c)

   This function frees all the memory used by the combination :data:`c`.

.. function:: int gsl_combination_memcpy (gsl_combination * dest, const gsl_combination * src)

   This function copies the elements of the combination :data:`src` into the
   combination :data:`dest`.  The two combinations must have the same size.

Accessing combination elements
==============================

The following function can be used to access the elements of a combination.

.. function:: size_t gsl_combination_get (const gsl_combination * c, const size_t i)

   This function returns the value of the :data:`i`-th element of the
   combination :data:`c`.  If :data:`i` lies outside the allowed range of 0 to
   :math:`k - 1` then the error handler is invoked and 0 is returned. |inlinefn|

Combination properties
======================

.. function:: size_t gsl_combination_n (const gsl_combination * c)

   This function returns the range (:math:`n`) of the combination c.

.. function:: size_t gsl_combination_k (const gsl_combination * c)

   This function returns the number of elements (:math:`k`) in the combination :data:`c`.

.. function:: size_t * gsl_combination_data (const gsl_combination * c)

   This function returns a pointer to the array of elements in the
   combination :data:`c`.

.. index::
   single: checking combination for validity
   single: testing combination for validity

.. function:: int gsl_combination_valid (gsl_combination * c)

   This function checks that the combination :data:`c` is valid.  The :data:`k`
   elements should lie in the range 0 to :math:`n - 1`, with each
   value occurring once at most and in increasing order.

Combination functions
=====================

.. index:: iterating through combinations

.. function:: int gsl_combination_next (gsl_combination * c)

   This function advances the combination :data:`c` to the next combination
   in lexicographic order and returns :macro:`GSL_SUCCESS`.  If no further
   combinations are available it returns :macro:`GSL_FAILURE` and leaves
   :data:`c` unmodified.  Starting with the first combination and
   repeatedly applying this function will iterate through all possible
   combinations of a given order.

.. function:: int gsl_combination_prev (gsl_combination * c)

   This function steps backwards from the combination :data:`c` to the
   previous combination in lexicographic order, returning
   :macro:`GSL_SUCCESS`.  If no previous combination is available it returns
   :macro:`GSL_FAILURE` and leaves :data:`c` unmodified.

Reading and writing combinations
================================

The library provides functions for reading and writing combinations to a
file as binary data or formatted text.

.. function:: int gsl_combination_fwrite (FILE * stream, const gsl_combination * c)

   This function writes the elements of the combination :data:`c` to the
   stream :data:`stream` in binary format.  The function returns
   :macro:`GSL_EFAILED` if there was a problem writing to the file.  Since the
   data is written in the native binary format it may not be portable
   between different architectures.

.. function:: int gsl_combination_fread (FILE * stream, gsl_combination * c)

   This function reads elements from the open stream :data:`stream` into the
   combination :data:`c` in binary format.  The combination :data:`c` must be
   preallocated with correct values of :math:`n` and :math:`k` since the
   function uses the size of :data:`c` to determine how many bytes to read.
   The function returns :macro:`GSL_EFAILED` if there was a problem reading
   from the file.  The data is assumed to have been written in the native
   binary format on the same architecture.

.. function:: int gsl_combination_fprintf (FILE * stream, const gsl_combination * c, const char * format)

   This function writes the elements of the combination :data:`c`
   line-by-line to the stream :data:`stream` using the format specifier
   :data:`format`, which should be suitable for a type of :code:`size_t`.
   In ISO C99 the type modifier :code:`z` represents :code:`size_t`, so
   :code:`"%zu\n"` is a suitable format [#f1]_.
   The function returns
   :macro:`GSL_EFAILED` if there was a problem writing to the file.

.. function:: int gsl_combination_fscanf (FILE * stream, gsl_combination * c)

   This function reads formatted data from the stream :data:`stream` into the
   combination :data:`c`.  The combination :data:`c` must be preallocated with
   correct values of :math:`n` and :math:`k` since the function uses the size of :data:`c` to
   determine how many numbers to read.  The function returns
   :macro:`GSL_EFAILED` if there was a problem reading from the file.

Examples
========

The example program below prints all subsets of the set
:math:`{0,1,2,3}` ordered by size.  Subsets of the same size are
ordered lexicographically.

.. include:: examples/combination.c
   :code:

Here is the output from the program,

.. include:: examples/combination.txt
   :code:

All 16 subsets are generated, and the subsets of each size are sorted
lexicographically.

References and Further Reading
==============================

Further information on combinations can be found in,

* Donald L. Kreher, Douglas R. Stinson, Combinatorial Algorithms:
  Generation, Enumeration and Search, 1998, CRC Press LLC, ISBN
  084933988X

.. rubric:: Footnotes

.. [#f1] In versions of the GNU C library prior to the ISO C99 standard, 
         the type modifier :code:`Z` was used instead.
