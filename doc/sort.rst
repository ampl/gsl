.. index:: sorting, heapsort

*******
Sorting
*******

This chapter describes functions for sorting data, both directly and
indirectly (using an index).  All the functions use the *heapsort*
algorithm.  Heapsort is an :math:`O(N \log N)` algorithm which operates
in-place and does not require any additional storage.  It also provides
consistent performance, the running time for its worst-case (ordered
data) being not significantly longer than the average and best cases.
Note that the heapsort algorithm does not preserve the relative ordering
of equal elements---it is an *unstable* sort.  However the resulting
order of equal elements will be consistent across different platforms
when using these functions.

Sorting objects
===============

The following function provides a simple alternative to the standard
library function :func:`qsort`.  It is intended for systems lacking
:func:`qsort`, not as a replacement for it.  The function :func:`qsort`
should be used whenever possible, as it will be faster and can provide
stable ordering of equal elements.  Documentation for :func:`qsort` is
available in the GNU C Library Reference Manual.

The functions described in this section are defined in the header file
:file:`gsl_heapsort.h`.

.. index::
   single: comparison functions, definition

.. function:: void gsl_heapsort (void * array, size_t count, size_t size, gsl_comparison_fn_t compare)

   This function sorts the :data:`count` elements of the array :data:`array`,
   each of size :data:`size`, into ascending order using the comparison
   function :data:`compare`.  The type of the comparison function is defined by

   .. type:: gsl_comparison_fn_t

      ::

         int (*gsl_comparison_fn_t) (const void * a, const void * b)

   A comparison function should return a negative integer if the first
   argument is less than the second argument, :code:`0` if the two arguments
   are equal and a positive integer if the first argument is greater than
   the second argument.

   For example, the following function can be used to sort doubles into
   ascending numerical order.

   ::

      int
      compare_doubles (const double * a, const double * b)
      {
        if (*a > *b)
          return 1;
        else if (*a < *b)
          return -1;
        else
          return 0;
      }

   The appropriate function call to perform the sort is::

      gsl_heapsort (array, count, sizeof(double), compare_doubles);

   Note that unlike :func:`qsort` the heapsort algorithm cannot be made into
   a stable sort by pointer arithmetic.  The trick of comparing pointers for
   equal elements in the comparison function does not work for the heapsort
   algorithm.  The heapsort algorithm performs an internal rearrangement of
   the data which destroys its initial ordering.

.. index:: indirect sorting

.. function:: int gsl_heapsort_index (size_t * p, const void * array, size_t count, size_t size, gsl_comparison_fn_t compare)

   This function indirectly sorts the :data:`count` elements of the array
   :data:`array`, each of size :data:`size`, into ascending order using the
   comparison function :data:`compare`.  The resulting permutation is stored
   in :data:`p`, an array of length :data:`n`.  The elements of :data:`p` give the
   index of the array element which would have been stored in that position
   if the array had been sorted in place.  The first element of :data:`p`
   gives the index of the least element in :data:`array`, and the last
   element of :data:`p` gives the index of the greatest element in
   :data:`array`.  The array itself is not changed.

Sorting vectors
===============

The following functions will sort the elements of an array or vector,
either directly or indirectly.  They are defined for all real and integer
types using the normal suffix rules.  For example, the :code:`float`
versions of the array functions are :func:`gsl_sort_float` and
:func:`gsl_sort_float_index`.  The corresponding vector functions are
:func:`gsl_sort_vector_float` and :func:`gsl_sort_vector_float_index`.  The
prototypes are available in the header files :file:`gsl_sort_float.h`
:file:`gsl_sort_vector_float.h`.  The complete set of prototypes can be
included using the header files :file:`gsl_sort.h` and
:file:`gsl_sort_vector.h`.

There are no functions for sorting complex arrays or vectors, since the
ordering of complex numbers is not uniquely defined.  To sort a complex
vector by magnitude compute a real vector containing the magnitudes
of the complex elements, and sort this vector indirectly.  The resulting
index gives the appropriate ordering of the original complex vector.

.. index::
   single: sorting vector elements
   single: vector, sorting elements of

.. function:: void gsl_sort (double * data, const size_t stride, size_t n)

   This function sorts the :data:`n` elements of the array :data:`data` with
   stride :data:`stride` into ascending numerical order.

.. function:: void gsl_sort2 (double * data1, const size_t stride1, double * data2, const size_t stride2, size_t n)

   This function sorts the :data:`n` elements of the array :data:`data1` with
   stride :data:`stride1` into ascending numerical order, while making the
   same rearrangement of the array :data:`data2` with stride :data:`stride2`,
   also of size :data:`n`.

.. function:: void gsl_sort_vector (gsl_vector * v)

   This function sorts the elements of the vector :data:`v` into ascending
   numerical order.

.. function:: void gsl_sort_vector2 (gsl_vector * v1, gsl_vector * v2)

   This function sorts the elements of the vector :data:`v1` into ascending
   numerical order, while making the same rearrangement of the vector :data:`v2`.

.. index::
   single: indirect sorting, of vector elements

.. function:: void gsl_sort_index (size_t * p, const double * data, size_t stride, size_t n)

   This function indirectly sorts the :data:`n` elements of the array
   :data:`data` with stride :data:`stride` into ascending order, storing the
   resulting permutation in :data:`p`.  The array :data:`p` must be allocated with
   a sufficient length to store the :data:`n` elements of the permutation.
   The elements of :data:`p` give the index of the array element which would
   have been stored in that position if the array had been sorted in place.
   The array :data:`data` is not changed.

.. function:: int gsl_sort_vector_index (gsl_permutation * p, const gsl_vector * v)

   This function indirectly sorts the elements of the vector :data:`v` into
   ascending order, storing the resulting permutation in :data:`p`.  The
   elements of :data:`p` give the index of the vector element which would
   have been stored in that position if the vector had been sorted in
   place.  The first element of :data:`p` gives the index of the least element
   in :data:`v`, and the last element of :data:`p` gives the index of the
   greatest element in :data:`v`.  The vector :data:`v` is not changed.

Selecting the k smallest or largest elements
============================================

The functions described in this section select the :math:`k` smallest
or largest elements of a data set of size :math:`N`.  The routines use an
:math:`O(kN)` direct insertion algorithm which is suited to subsets that
are small compared with the total size of the dataset. For example, the
routines are useful for selecting the 10 largest values from one million
data points, but not for selecting the largest 100,000 values.  If the
subset is a significant part of the total dataset it may be faster
to sort all the elements of the dataset directly with an :math:`O(N \log N)`
algorithm and obtain the smallest or largest values that way.

.. function:: int gsl_sort_smallest (double * dest, size_t k, const double * src, size_t stride, size_t n)

   This function copies the :data:`k` smallest elements of the array
   :data:`src`, of size :data:`n` and stride :data:`stride`, in ascending
   numerical order into the array :data:`dest`.  The size :data:`k` of the subset must be
   less than or equal to :data:`n`.  The data :data:`src` is not modified by
   this operation.

.. function:: int gsl_sort_largest (double * dest, size_t k, const double * src, size_t stride, size_t n)

   This function copies the :data:`k` largest elements of the array
   :data:`src`, of size :data:`n` and stride :data:`stride`, in descending
   numerical order into the array :data:`dest`. :data:`k` must be
   less than or equal to :data:`n`. The data :data:`src` is not modified by
   this operation.

.. function:: int gsl_sort_vector_smallest (double * dest, size_t k, const gsl_vector * v)
              int gsl_sort_vector_largest (double * dest, size_t k, const gsl_vector * v)

   These functions copy the :data:`k` smallest or largest elements of the
   vector :data:`v` into the array :data:`dest`. :data:`k`
   must be less than or equal to the length of the vector :data:`v`.

The following functions find the indices of the :math:`k` smallest or
largest elements of a dataset.

.. function:: int gsl_sort_smallest_index (size_t * p, size_t k, const double * src, size_t stride, size_t n)

   This function stores the indices of the :data:`k` smallest elements of
   the array :data:`src`, of size :data:`n` and stride :data:`stride`, in the
   array :data:`p`.  The indices are chosen so that the corresponding data is
   in ascending numerical order.  :data:`k` must be
   less than or equal to :data:`n`. The data :data:`src` is not modified by
   this operation.

.. function:: int gsl_sort_largest_index (size_t * p, size_t k, const double * src, size_t stride, size_t n)

   This function stores the indices of the :data:`k` largest elements of
   the array :data:`src`, of size :data:`n` and stride :data:`stride`, in the
   array :data:`p`.  The indices are chosen so that the corresponding data is
   in descending numerical order.  :data:`k` must be
   less than or equal to :data:`n`. The data :data:`src` is not modified by
   this operation.

.. function:: int gsl_sort_vector_smallest_index (size_t * p, size_t k, const gsl_vector * v)
              int gsl_sort_vector_largest_index (size_t * p, size_t k, const gsl_vector * v)

   These functions store the indices of the :data:`k` smallest or largest
   elements of the vector :data:`v` in the array :data:`p`. :data:`k` must be less than or equal to the length of the vector
   :data:`v`.

Computing the rank
==================

The *rank* of an element is its order in the sorted data.  The rank
is the inverse of the index permutation, :math:`p`.  It can be computed
using the following algorithm::

  for (i = 0; i < p->size; i++) 
    {
      size_t pi = p->data[i];
      rank->data[pi] = i;
    }

This can be computed directly from the function
:code:`gsl_permutation_inverse(rank,p)`.

The following function will print the rank of each element of the vector
:math:`v`::

  void
  print_rank (gsl_vector * v)
  {
    size_t i;
    size_t n = v->size;
    gsl_permutation * perm = gsl_permutation_alloc(n);
    gsl_permutation * rank = gsl_permutation_alloc(n);

    gsl_sort_vector_index (perm, v);
    gsl_permutation_inverse (rank, perm);

    for (i = 0; i < n; i++)
      {
        double vi = gsl_vector_get(v, i);
        printf ("element = %d, value = %g, rank = %d\n",
                 i, vi, rank->data[i]);
      }

    gsl_permutation_free (perm);
    gsl_permutation_free (rank);
  }

Examples
========

The following example shows how to use the permutation :math:`p` to print
the elements of the vector :math:`v` in ascending order::

  gsl_sort_vector_index (p, v);

  for (i = 0; i < v->size; i++)
    {
      double vpi = gsl_vector_get (v, p->data[i]);
      printf ("order = %d, value = %g\n", i, vpi);
    }

The next example uses the function :func:`gsl_sort_smallest` to select
the 5 smallest numbers from 100000 uniform random variates stored in an
array,

.. include:: examples/sortsmall.c
   :code:

The output lists the 5 smallest values, in ascending order,

.. include:: examples/sortsmall.txt
   :code:

References and Further Reading
==============================

The subject of sorting is covered extensively in the following,

* Donald E. Knuth, The Art of Computer Programming: Sorting and
  Searching (Vol 3, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850.

The Heapsort algorithm is described in the following book,

* Robert Sedgewick, Algorithms in C, Addison-Wesley, 
  ISBN 0201514257.
