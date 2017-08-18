.. index::
   single: sparse matrices
   single: matrices, sparse

***************
Sparse Matrices
***************

This chapter describes functions for the construction and
manipulation of sparse matrices, matrices which are populated
primarily with zeros and contain only a few non-zero elements.
Sparse matrices often appear in the solution of partial
differential equations. It is beneficial to use specialized
data structures and algorithms for storing and working with
sparse matrices, since dense matrix algorithms and structures
can be very slow and use huge amounts of memory when applied
to sparse matrices.

The header file :file:`gsl_spmatrix.h` contains the prototypes for the
sparse matrix functions and related declarations.

.. index::
   single: sparse matrices, overview

Overview
========

These routines provide support for constructing and manipulating
sparse matrices in GSL, using an API similar to :type:`gsl_matrix`.
The basic structure is called :type:`gsl_spmatrix`. There are
three supported storage formats for sparse matrices: the triplet,
compressed column storage (CCS) and compressed row storage (CRS)
formats. The triplet format stores triplets :math:`(i,j,x)` for each
non-zero element of the matrix. This notation means that the
:math:`(i,j)` element of the matrix :math:`A`
is :math:`A_{ij} = x`. Compressed column storage stores each column of
non-zero values in the sparse matrix in a continuous memory block, keeping
pointers to the beginning of each column in that memory block, and storing
the row indices of each non-zero element. Compressed row storage stores
each row of non-zero values in a continuous memory block, keeping pointers
to the beginning of each row in the block and storing the column indices
of each non-zero element. The triplet format is ideal
for adding elements to the sparse matrix structure while it is being
constructed, while the compressed storage formats are better suited for
matrix-matrix multiplication or linear solvers.

.. type:: gsl_spmatrix

   This structure is defined as::

      typedef struct
      {
        size_t size1;
        size_t size2;
        size_t *i;
        double *data;
        size_t *p;
        size_t nzmax;
        size_t nz;
        gsl_spmatrix_tree *tree_data;
        void *work;
        size_t sptype;
      } gsl_spmatrix;

   This defines a :data:`size1`-by-:data:`size2` sparse matrix. The number of non-zero
   elements currently in the matrix is given by :data:`nz`. For the triplet
   representation, :data:`i`, :data:`p`, and :data:`data` are arrays of size :data:`nz`
   which contain the row indices, column indices, and element value, respectively.
   So if :math:`data[k] = A(i,j)`, then :math:`i = i[k]` and :math:`j = p[k]`.

   For compressed column storage, :data:`i` and :data:`data` are arrays of size
   :data:`nz` containing the row indices and element values, identical to the triplet
   case. :data:`p` is an array of size :data:`size2` + 1 where :code:`p[j]` points
   to the index in :data:`data` of the start of column :data:`j`. Thus, if
   :math:`data[k] = A(i,j)`, then :math:`i = i[k]` and :math:`p[j] <= k < p[j+1]`.

   For compressed row storage, :data:`i` and :data:`data` are arrays of size
   :data:`nz` containing the column indices and element values, identical to the triplet
   case. :data:`p` is an array of size :data:`size1` + 1 where :code:`p[i]` points
   to the index in :data:`data` of the start of row :data:`i`. Thus, if
   :math:`data[k] = A(i,j)`, then :math:`j = i[k]` and :math:`p[i] <= k < p[i+1]`.

   The parameter :data:`tree_data` is a binary tree structure used in the triplet
   representation, specifically a balanced AVL tree. This speeds up element
   searches and duplicate detection during the matrix assembly process.
   The parameter :data:`work` is additional workspace needed for various operations like
   converting from triplet to compressed storage. :data:`sptype` indicates
   the type of storage format being used (triplet, CCS or CRS).

   The compressed storage format defined above makes it very simple
   to interface with sophisticated external linear solver libraries
   which accept compressed storage input. The user can simply
   pass the arrays :data:`i`, :data:`p`, and :data:`data` as the
   inputs to external libraries.

.. index::
   single: sparse matrices, allocation

Allocation
==========

The functions for allocating memory for a sparse matrix follow the style of
:func:`malloc` and :func:`free`. They also perform their own error checking. If
there is insufficient memory available to allocate a matrix then the functions
call the GSL error handler with an error code of :macro:`GSL_ENOMEM` in addition
to returning a null pointer.

.. function:: gsl_spmatrix * gsl_spmatrix_alloc (const size_t n1, const size_t n2)

   This function allocates a sparse matrix of size :data:`n1`-by-:data:`n2` and
   initializes it to all zeros. If the size of the matrix is not known at allocation
   time, both :data:`n1` and :data:`n2` may be set to 1, and they will automatically
   grow as elements are added to the matrix. This function sets the
   matrix to the triplet representation, which is the easiest for adding
   and accessing matrix elements. This function tries to make a reasonable guess
   for the number of non-zero elements (:data:`nzmax`) which will be added to the matrix by
   assuming a sparse density of :math:`10\%`. The function
   :func:`gsl_spmatrix_alloc_nzmax` can be used if this number is known more
   accurately. The workspace is of size :math:`O(nzmax)`.

.. function:: gsl_spmatrix * gsl_spmatrix_alloc_nzmax (const size_t n1, const size_t n2, const size_t nzmax, const size_t sptype)

   This function allocates a sparse matrix of size :data:`n1`-by-:data:`n2` and
   initializes it to all zeros. If the size of the matrix is not known at allocation
   time, both :data:`n1` and :data:`n2` may be set to 1, and they will automatically
   grow as elements are added to the matrix. The parameter :data:`nzmax` specifies
   the maximum number of non-zero elements which will be added to the matrix.
   It does not need to be precisely known in advance, since storage space will 
   automatically grow using :func:`gsl_spmatrix_realloc` if :data:`nzmax` is not
   large enough. Accurate knowledge of this parameter reduces the number of
   reallocation calls required. The parameter :data:`sptype` specifies the
   storage format of the sparse matrix. Possible values are

   .. macro:: GSL_SPMATRIX_TRIPLET

      This flag specifies triplet storage.

   .. macro:: GSL_SPMATRIX_CCS

      This flag specifies compressed column storage.

   .. macro:: GSL_SPMATRIX_CRS

      This flag specifies compressed row storage.

   The allocated :type:`gsl_spmatrix` structure is of size :math:`O(nzmax)`.

.. function:: int gsl_spmatrix_realloc (const size_t nzmax, gsl_spmatrix * m)

   This function reallocates the storage space for :data:`m` to accomodate
   :data:`nzmax` non-zero elements. It is typically called internally by
   :func:`gsl_spmatrix_set` if the user wants to add more elements to the
   sparse matrix than the previously specified :data:`nzmax`.

.. function:: void gsl_spmatrix_free (gsl_spmatrix * m)

   This function frees the memory associated with the sparse matrix :data:`m`.

.. index::
   single: sparse matrices, accessing elements

Accessing Matrix Elements
=========================

.. function:: double gsl_spmatrix_get (const gsl_spmatrix * m, const size_t i, const size_t j)

   This function returns element (:data:`i`, :data:`j`) of the matrix :data:`m`.
   The matrix may be in triplet or compressed format.

.. function:: int gsl_spmatrix_set (gsl_spmatrix * m, const size_t i, const size_t j, const double x)

   This function sets element (:data:`i`, :data:`j`) of the matrix :data:`m` to
   the value :data:`x`. The matrix must be in triplet representation.

.. function:: double * gsl_spmatrix_ptr (gsl_spmatrix * m, const size_t i, const size_t j)

   This function returns a pointer to the (:data:`i`, :data:`j`) element of the matrix :data:`m`.
   If the (:data:`i`, :data:`j`) element is not explicitly stored in the matrix,
   a null pointer is returned.

.. index::
   single: sparse matrices, initializing elements

Initializing Matrix Elements
============================

Since the sparse matrix format only stores the non-zero elements, it is automatically
initialized to zero upon allocation. The function :func:`gsl_spmatrix_set_zero` may
be used to re-initialize a matrix to zero after elements have been added to it.

.. function:: int gsl_spmatrix_set_zero (gsl_spmatrix * m)

   This function sets (or resets) all the elements of the matrix :data:`m` to zero.

.. index::
   single: sparse matrices, reading
   single: sparse matrices, writing

Reading and Writing Matrices
============================

.. function:: int gsl_spmatrix_fwrite (FILE * stream, const gsl_spmatrix * m)

   This function writes the elements of the matrix :data:`m` to the stream
   :data:`stream` in binary format.  The return value is 0 for success and
   :macro:`GSL_EFAILED` if there was a problem writing to the file.  Since the
   data is written in the native binary format it may not be portable
   between different architectures.

.. function:: int gsl_spmatrix_fread (FILE * stream, gsl_spmatrix * m)

   This function reads into the matrix :data:`m` from the open stream
   :data:`stream` in binary format.  The matrix :data:`m` must be preallocated
   with the correct storage format, dimensions and have a sufficiently large :data:`nzmax`
   in order to read in all matrix elements, otherwise :macro:`GSL_EBADLEN`
   is returned. The return value is 0 for success and
   :macro:`GSL_EFAILED` if there was a problem reading from the file.  The
   data is assumed to have been written in the native binary format on the
   same architecture.

.. function:: int gsl_spmatrix_fprintf (FILE * stream, const gsl_spmatrix * m, const char * format)

   This function writes the elements of the matrix :data:`m` line-by-line to
   the stream :data:`stream` using the format specifier :data:`format`, which
   should be one of the :code:`%g`, :code:`%e` or :code:`%f` formats for
   floating point numbers.  The function returns 0 for success and
   :macro:`GSL_EFAILED` if there was a problem writing to the file. The
   input matrix :data:`m` may be in any storage format, and the output file
   will be written in MatrixMarket format.

.. function:: gsl_spmatrix * gsl_spmatrix_fscanf (FILE * stream)

   This function reads sparse matrix data in the MatrixMarket format
   from the stream :data:`stream` and stores it in a newly allocated matrix
   which is returned in triplet format.  The function returns 0 for success and
   :macro:`GSL_EFAILED` if there was a problem reading from the file. The
   user should free the returned matrix when it is no longer needed.

.. index::
   single: sparse matrices, copying

Copying Matrices
================

.. function:: int gsl_spmatrix_memcpy (gsl_spmatrix * dest, const gsl_spmatrix * src)

   This function copies the elements of the sparse matrix :data:`src` into
   :data:`dest`. The two matrices must have the same dimensions and be in the
   same storage format.

.. index::
   single: sparse matrices, exchanging rows and columns

Exchanging Rows and Columns
===========================

.. function:: int gsl_spmatrix_transpose_memcpy (gsl_spmatrix * dest, const gsl_spmatrix * src)

   This function copies the transpose of the sparse matrix :data:`src` into
   :data:`dest`. The dimensions of :data:`dest` must match the transpose of the
   matrix :data:`src`. Also, both matrices must use the same sparse storage
   format.

.. function:: int gsl_spmatrix_transpose (gsl_spmatrix * m)

   This function replaces the matrix :data:`m` by its transpose,
   preserving the storage format of the input matrix. Currently,
   only triplet matrix inputs are supported.

.. function:: int gsl_spmatrix_transpose2 (gsl_spmatrix * m)

   This function replaces the matrix :data:`m` by its transpose, but
   changes the storage format for compressed matrix inputs. Since
   compressed column storage is the transpose of compressed row storage,
   this function simply converts a CCS matrix to CRS and vice versa.
   This is the most efficient way to transpose a compressed storage
   matrix, but the user should note that the storage format of their
   compressed matrix will change on output. For triplet matrices,
   the output matrix is also in triplet storage.

.. index::
   single: sparse matrices, operations

Matrix Operations
=================

.. function:: int gsl_spmatrix_add (gsl_spmatrix * c, const gsl_spmatrix * a, const gsl_spmatrix * b)

   This function computes the sum :math:`c = a + b`. The three matrices must
   have the same dimensions and be stored in a compressed format.

.. function:: int gsl_spmatrix_scale (gsl_spmatrix * m, const double x)

   This function scales all elements of the matrix :data:`m` by the constant
   factor :data:`x`. The result :math:`m(i,j) \leftarrow x m(i,j)` is stored in :data:`m`.

.. index::
   single: sparse matrices, properties

Matrix Properties
=================

.. function:: size_t gsl_spmatrix_nnz (const gsl_spmatrix * m)

   This function returns the number of non-zero elements in :data:`m`.

.. function:: int gsl_spmatrix_equal (const gsl_spmatrix * a, const gsl_spmatrix * b)

   This function returns 1 if the matrices :data:`a` and :data:`b` are equal (by comparison of
   element values) and 0 otherwise. The matrices :data:`a` and :data:`b` must be in the same
   sparse storage format for comparison.

.. index::
   single: sparse matrices, min/max elements

Finding Maximum and Minimum Elements
====================================

.. function:: int gsl_spmatrix_minmax (const gsl_spmatrix * m, double * min_out, double * max_out)

   This function returns the minimum and maximum elements of the matrix
   :data:`m`, storing them in :data:`min_out` and :data:`max_out`, and searching
   only the non-zero values.

.. index::
   single: sparse matrices, compression

Compressed Format
=================

GSL supports compressed column storage (CCS) and compressed row storage (CRS)
formats.

.. function:: gsl_spmatrix * gsl_spmatrix_ccs (const gsl_spmatrix * T)

   This function creates a sparse matrix in compressed column format
   from the input sparse matrix :data:`T` which must be in triplet format.
   A pointer to a newly allocated matrix is returned. The calling function
   should free the newly allocated matrix when it is no longer needed.

.. function:: gsl_spmatrix * gsl_spmatrix_crs (const gsl_spmatrix * T)

   This function creates a sparse matrix in compressed row format
   from the input sparse matrix :data:`T` which must be in triplet format.
   A pointer to a newly allocated matrix is returned. The calling function
   should free the newly allocated matrix when it is no longer needed.

.. index::
   single: sparse matrices, conversion

Conversion Between Sparse and Dense Matrices
============================================

The :type:`gsl_spmatrix` structure can be converted into the dense :type:`gsl_matrix`
format and vice versa with the following routines.

.. function:: int gsl_spmatrix_d2sp (gsl_spmatrix * S, const gsl_matrix * A)

   This function converts the dense matrix :data:`A` into sparse triplet format
   and stores the result in :data:`S`.

.. function:: int gsl_spmatrix_sp2d (gsl_matrix * A, const gsl_spmatrix * S)

   This function converts the sparse matrix :data:`S` into a dense matrix and
   stores the result in :data:`A`. :data:`S` must be in triplet format.

.. index::
   single: sparse matrices, examples

Examples
========

The following example program builds a 5-by-4 sparse matrix
and prints it in triplet, compressed column, and compressed
row format. The matrix which is constructed is

.. only:: not texinfo

   .. math::

      \left(
        \begin{array}{cccc}
          0 & 0 & 3.1 & 4.6 \\
          1 & 0 & 7.2 & 0 \\
          0 & 0 & 0 & 0 \\
          2.1 & 2.9 & 0 & 8.5 \\
          4.1 & 0 & 0 & 0
        \end{array}
      \right)

.. only:: texinfo

   ::

     [ 0    0  3.1  4.6 ]
     [ 1    0  7.2   0  ]
     [ 0    0   0    0  ]
     [ 2.1 2.9  0   8.5 ]
     [ 4.1  0   0    0  ]

The output of the program is::

  printing all matrix elements:
  A(0,0) = 0
  A(0,1) = 0
  A(0,2) = 3.1
  A(0,3) = 4.6
  A(1,0) = 1
  .
  .
  .
  A(4,0) = 4.1
  A(4,1) = 0
  A(4,2) = 0
  A(4,3) = 0
  matrix in triplet format (i,j,Aij):
  (0, 2, 3.1)
  (0, 3, 4.6)
  (1, 0, 1.0)
  (1, 2, 7.2)
  (3, 0, 2.1)
  (3, 1, 2.9)
  (3, 3, 8.5)
  (4, 0, 4.1)
  matrix in compressed column format:
  i = [ 1, 3, 4, 3, 0, 1, 0, 3, ]
  p = [ 0, 3, 4, 6, 8, ]
  d = [ 1, 2.1, 4.1, 2.9, 3.1, 7.2, 4.6, 8.5, ]
  matrix in compressed row format:
  i = [ 2, 3, 0, 2, 0, 1, 3, 0, ]
  p = [ 0, 2, 4, 4, 7, 8, ]
  d = [ 3.1, 4.6, 1, 7.2, 2.1, 2.9, 8.5, 4.1, ]

We see in the compressed column output, the data array stores
each column contiguously, the array :math:`i` stores
the row index of the corresponding data element, and the
array :math:`p` stores the index of the start of each column in the
data array. Similarly, for the compressed row output, the
data array stores each row contiguously, the array :math:`i`
stores the column index of the corresponding data element, and
the :math:`p` array stores the index of the start of each row
in the data array.

.. include:: examples/spmatrix.c
   :code:

.. index::
   single: sparse matrices, references

References and Further Reading
==============================

The algorithms used by these functions are described in the
following sources,

* Davis, T. A., Direct Methods for Sparse Linear Systems, SIAM, 2006.

* CSparse software library, https://www.cise.ufl.edu/research/sparse/CSparse
