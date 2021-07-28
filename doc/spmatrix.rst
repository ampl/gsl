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
can be prohibitively slow and use huge amounts of memory when applied
to sparse matrices.

The header file :file:`gsl_spmatrix.h` contains the prototypes for the
sparse matrix functions and related declarations.

.. index::
   single: sparse matrices, data types

Data types
==========

All the functions are available for each of the standard data-types.
The versions for :code:`double` have the prefix :code:`gsl_spmatrix`,
Similarly the versions for single-precision :code:`float` arrays have the prefix
:code:`gsl_spmatrix_float`.  The full list of available types is given
below,

================================ ====================
Prefix                           Type
================================ ====================
gsl_spmatrix                     double         
gsl_spmatrix_float               float         
gsl_spmatrix_long_double         long double   
gsl_spmatrix_int                 int           
gsl_spmatrix_uint                unsigned int  
gsl_spmatrix_long                long          
gsl_spmatrix_ulong               unsigned long 
gsl_spmatrix_short               short         
gsl_spmatrix_ushort              unsigned short
gsl_spmatrix_char                char          
gsl_spmatrix_uchar               unsigned char 
gsl_spmatrix_complex             complex double        
gsl_spmatrix_complex_float       complex float         
gsl_spmatrix_complex_long_double complex long double   
================================ ====================

.. index::
   single: sparse matrices, storage formats

Sparse Matrix Storage Formats
=============================

GSL currently supports three storage formats for sparse matrices:
the coordinate (COO) representation, compressed sparse column (CSC)
and compressed sparse row (CSR) formats. These are discussed in more
detail below. In order to illustrate the different storage formats,
the following sections will reference this :math:`M`-by-:math:`N`
sparse matrix, with :math:`M=4` and :math:`N=5`:

.. math:: \begin{pmatrix}
            9 & 0 & 0 & 0 & -3 \\
            4 & 7 & 0 & 0 & 0 \\
            0 & 8 & -1 & 8 & 0 \\
            4 & 0 & 5 & 6 & 0
          \end{pmatrix}

The number of non-zero elements in the matrix, also abbreviated as
:code:`nnz` is equal to :math:`10` in this case.

.. index::
   single: sparse matrices, coordinate format
   single: sparse matrices, triplet format

.. _sec_spmatrix-coo:

Coordinate Storage (COO)
------------------------

The coordinate storage format, also known as *triplet format*,
stores triplets :math:`(i,j,x)` for each non-zero element of the matrix.
This notation means that the :math:`(i,j)` element of the matrix :math:`A`
is :math:`A_{ij} = x`. The matrix is stored using three arrays of the same
length, representing the row indices, column indices, and matrix data. For
the reference matrix above, one possible storage format is:

==== == == == == == == == == == ==
data  9  7  4  8 -3 -1  8  5  6  4 
row   0  1  1  2  0  2  2  3  3  3
col   0  1  0  1  4  2  3  2  3  0
==== == == == == == == == == == ==

Note that this representation is not unique - the coordinate triplets may
appear in any ordering and would still represent the same sparse matrix.
The length of the three arrays is equal to the number of non-zero elements
in the matrix, :code:`nnz`, which in this case is :math:`10`. The coordinate format is
extremely convenient for sparse matrix *assembly*,
the process of adding new elements, or changing existing elements, in a
sparse matrix. However, it is generally not suitable for the efficient
implementation of matrix-matrix products, or matrix factorization algorithms.
For these applications it is better to use one of the compressed formats
discussed below.

In order to faciliate efficient sparse matrix assembly, GSL stores
the coordinate data in a balanced binary search tree, specifically an AVL
tree, in addition to the three arrays
discussed above. This allows GSL to efficiently determine whether an
entry :math:`(i,j)` already exists in the matrix, and to replace an existing
matrix entry with a new value, without needing to search unsorted arrays.

.. index::
   single: sparse matrices, compressed sparse column
   single: sparse matrices, compressed column storage

.. _sec_spmatrix-csc:

Compressed Sparse Column (CSC)
------------------------------

Compressed sparse column storage stores each column of
non-zero values in the sparse matrix in a continuous memory block, keeping
pointers to the beginning of each column in that memory block, and storing
the row indices of each non-zero element. For the reference matrix above,
these arrays look like

======= == == == == == == == == == ==
data     9  4  4  7  8 -1  5  8  6 -3
row      0  1  3  1  2  2  3  2  3  0
col_ptr  0  3  5  7  9 10
======= == == == == == == == == == ==

The :code:`data` and :code:`row` arrays are of length :code:`nnz` and
are the same as the COO storage format. The :code:`col_ptr` array
has length :math:`N+1`, and :code:`col_ptr[j]` gives the index in
:code:`data` of the start of column :code:`j`. Therefore, the
:math:`j`-th column of the matrix is stored in
:code:`data[col_ptr[j]]`, :code:`data[col_ptr[j] + 1]`, ...,
:code:`data[col_ptr[j+1] - 1]`.
The last element of :code:`col_ptr` is :code:`nnz`.

.. index::
   single: sparse matrices, compressed sparse row
   single: sparse matrices, compressed row storage

.. _sec_spmatrix-csr:

Compressed Sparse Row (CSR)
---------------------------

Compressed row storage stores each row of non-zero values in a
continuous memory block, keeping pointers to the beginning of each
row in the block and storing the column indices of each non-zero element.
For the reference matrix above, these arrays look like

======= == == == == == == == == == ==
data     9 -3  4  7  8 -1  8  4  5  6
col      0  4  0  1  1  2  3  0  2  3
row_ptr  0  2  4  7 10
======= == == == == == == == == == ==

The :code:`data` and :code:`col` arrays are of length :code:`nnz` and
are the same as the COO storage format. The :code:`row_ptr` array
has length :math:`M+1`, and :code:`row_ptr[i]` gives the index in
:code:`data` of the start of row :code:`i`. Therefore, the
:math:`i`-th row of the matrix is stored in
:code:`data[row_ptr[i]]`, :code:`data[row_ptr[i] + 1]`,
..., :code:`data[row_ptr[i+1] - 1]`.
The last element of :code:`row_ptr` is :code:`nnz`.

.. index::
   single: sparse matrices, overview

Overview
========

These routines provide support for constructing and manipulating
sparse matrices in GSL, using an API similar to :type:`gsl_matrix`.
The basic structure is called :type:`gsl_spmatrix`.

.. type:: gsl_spmatrix

   This structure is defined as::

      typedef struct
      {
        size_t size1;
        size_t size2;
        int *i;
        double *data;
        int *p;
        size_t nzmax;
        size_t nz;
        [ ... variables for binary tree and memory management ... ]
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

   There are additional variables in the :type:`gsl_spmatrix` structure related
   to binary tree storage and memory management. The GSL implementation of sparse
   matrices uses balanced AVL trees to sort matrix elements in the triplet representation.
   This speeds up element searches and duplicate detection during the matrix assembly process.
   The :type:`gsl_spmatrix` structure also contains additional workspace variables needed
   for various operations like converting from triplet to compressed storage.
   :data:`sptype` indicates the type of storage format being used (COO, CSC or CSR).

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

   .. macro:: GSL_SPMATRIX_COO

      This flag specifies coordinate (triplet) storage.

   .. macro:: GSL_SPMATRIX_CSC

      This flag specifies compressed sparse column storage.

   .. macro:: GSL_SPMATRIX_CSR

      This flag specifies compressed sparse row storage.

   The allocated :type:`gsl_spmatrix` structure is of size :math:`O(nzmax)`.

.. function:: int gsl_spmatrix_realloc (const size_t nzmax, gsl_spmatrix * m)

   This function reallocates the storage space for :data:`m` to accomodate
   :data:`nzmax` non-zero elements. It is typically called internally by
   :func:`gsl_spmatrix_set` if the user wants to add more elements to the
   sparse matrix than the previously specified :data:`nzmax`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: void gsl_spmatrix_free (gsl_spmatrix * m)

   This function frees the memory associated with the sparse matrix :data:`m`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. index::
   single: sparse matrices, accessing elements

Accessing Matrix Elements
=========================

.. function:: double gsl_spmatrix_get (const gsl_spmatrix * m, const size_t i, const size_t j)

   This function returns element (:data:`i`, :data:`j`) of the matrix :data:`m`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_set (gsl_spmatrix * m, const size_t i, const size_t j, const double x)

   This function sets element (:data:`i`, :data:`j`) of the matrix :data:`m` to
   the value :data:`x`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`

.. function:: double * gsl_spmatrix_ptr (gsl_spmatrix * m, const size_t i, const size_t j)

   This function returns a pointer to the (:data:`i`, :data:`j`) element of the matrix :data:`m`.
   If the (:data:`i`, :data:`j`) element is not explicitly stored in the matrix,
   a null pointer is returned.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. index::
   single: sparse matrices, initializing elements

Initializing Matrix Elements
============================

Since the sparse matrix format only stores the non-zero elements, it is automatically
initialized to zero upon allocation. The function :func:`gsl_spmatrix_set_zero` may
be used to re-initialize a matrix to zero after elements have been added to it.

.. function:: int gsl_spmatrix_set_zero (gsl_spmatrix * m)

   This function sets (or resets) all the elements of the matrix :data:`m` to zero.
   For CSC and CSR matrices, the cost of this operation is :math:`O(1)`. For
   COO matrices, the binary tree structure must be dismantled, so the cost is
   :math:`O(nz)`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

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

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_fread (FILE * stream, gsl_spmatrix * m)

   This function reads into the matrix :data:`m` from the open stream
   :data:`stream` in binary format.  The matrix :data:`m` must be preallocated
   with the correct storage format, dimensions and have a sufficiently large :data:`nzmax`
   in order to read in all matrix elements, otherwise :macro:`GSL_EBADLEN`
   is returned. The return value is 0 for success and
   :macro:`GSL_EFAILED` if there was a problem reading from the file.  The
   data is assumed to have been written in the native binary format on the
   same architecture.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_fprintf (FILE * stream, const gsl_spmatrix * m, const char * format)

   This function writes the elements of the matrix :data:`m` line-by-line to
   the stream :data:`stream` using the format specifier :data:`format`, which
   should be one of the :code:`%g`, :code:`%e` or :code:`%f` formats for
   floating point numbers.  The function returns 0 for success and
   :macro:`GSL_EFAILED` if there was a problem writing to the file. The
   input matrix :data:`m` may be in any storage format, and the output file
   will be written in MatrixMarket format.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: gsl_spmatrix * gsl_spmatrix_fscanf (FILE * stream)

   This function reads sparse matrix data in the MatrixMarket format
   from the stream :data:`stream` and stores it in a newly allocated matrix
   which is returned in :ref:`COO <sec_spmatrix-coo>` format.  The function returns 0 for success and
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

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. index::
   single: sparse matrices, exchanging rows and columns

Exchanging Rows and Columns
===========================

.. function:: int gsl_spmatrix_transpose_memcpy (gsl_spmatrix * dest, const gsl_spmatrix * src)

   This function copies the transpose of the sparse matrix :data:`src` into
   :data:`dest`. The dimensions of :data:`dest` must match the transpose of the
   matrix :data:`src`. Also, both matrices must use the same sparse storage
   format.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_transpose (gsl_spmatrix * m)

   This function replaces the matrix :data:`m` by its transpose, but
   changes the storage format for compressed matrix inputs. Since
   compressed column storage is the transpose of compressed row storage,
   this function simply converts a CSC matrix to CSR and vice versa.
   This is the most efficient way to transpose a compressed storage
   matrix, but the user should note that the storage format of their
   compressed matrix will change on output. For COO matrix inputs,
   the output matrix is also in COO storage.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. index::
   single: sparse matrices, operations

Matrix Operations
=================

.. function:: int gsl_spmatrix_scale (gsl_spmatrix * m, const double x)

   This function scales all elements of the matrix :data:`m` by the constant
   factor :data:`x`. The result :math:`m(i,j) \leftarrow x m(i,j)` is stored in :data:`m`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_scale_columns (gsl_spmatrix * A, const gsl_vector * x)

   This function scales the columns of the :math:`M`-by-:math:`N` sparse matrix
   :data:`A` by the elements of the vector :data:`x`, of length :math:`N`. The
   :math:`j`-th column of :data:`A` is multiplied by :code:`x[j]`. This is equivalent to
   forming

   .. math:: A \rightarrow A X

   where :math:`X = \textrm{diag}(x)`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_scale_rows (gsl_spmatrix * A, const gsl_vector * x)

   This function scales the rows of the :math:`M`-by-:math:`N` sparse matrix
   :data:`A` by the elements of the vector :data:`x`, of length :math:`M`. The
   :math:`i`-th row of :data:`A` is multiplied by :code:`x[i]`. This is equivalent to
   forming

   .. math:: A \rightarrow X A

   where :math:`X = \textrm{diag}(x)`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_add (gsl_spmatrix * c, const gsl_spmatrix * a, const gsl_spmatrix * b)

   This function computes the sum :math:`c = a + b`. The three matrices must
   have the same dimensions.

   Input matrix formats supported: :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_dense_add (gsl_matrix * a, const gsl_spmatrix * b)

   This function adds the elements of the sparse matrix :data:`b` to the elements of
   the dense matrix :data:`a`. The result :math:`a(i,j) \leftarrow a(i,j) + b(i,j)` is
   stored in :data:`a` and :data:`b` remains unchanged. The two matrices must have
   the same dimensions.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_dense_sub (gsl_matrix * a, const gsl_spmatrix * b)

   This function subtracts the elements of the sparse matrix :data:`b` from the elements of
   the dense matrix :data:`a`. The result :math:`a(i,j) \leftarrow a(i,j) - b(i,j)` is
   stored in :data:`a` and :data:`b` remains unchanged. The two matrices must have
   the same dimensions.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. index::
   single: sparse matrices, properties

Matrix Properties
=================

.. function:: const char * gsl_spmatrix_type (const gsl_spmatrix * m)

   This function returns a string describing the sparse storage format of
   the matrix :data:`m`. For example::

      printf ("matrix is '%s' format.\n", gsl_spmatrix_type (m));

   would print something like::

      matrix is 'CSR' format.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: size_t gsl_spmatrix_nnz (const gsl_spmatrix * m)

   This function returns the number of non-zero elements in :data:`m`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_equal (const gsl_spmatrix * a, const gsl_spmatrix * b)

   This function returns 1 if the matrices :data:`a` and :data:`b` are equal (by comparison of
   element values) and 0 otherwise. The matrices :data:`a` and :data:`b` must be in the same
   sparse storage format for comparison.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: double gsl_spmatrix_norm1 (const gsl_spmatrix * A)

   This function returns the 1-norm of the :math:`m`-by-:math:`n` matrix :data:`A`, defined as
   the maximum column sum,

   .. math:: ||A||_1 = \textrm{max}_{1 \le j \le n} \sum_{i=1}^m |A_{ij}|

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. index::
   single: sparse matrices, min/max elements

Finding Maximum and Minimum Elements
====================================

.. function:: int gsl_spmatrix_minmax (const gsl_spmatrix * m, double * min_out, double * max_out)

   This function returns the minimum and maximum elements of the matrix
   :data:`m`, storing them in :data:`min_out` and :data:`max_out`, and searching
   only the non-zero values.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. function:: int gsl_spmatrix_min_index (const gsl_spmatrix * m, size_t * imin, size_t * jmin)

   This function returns the indices of the minimum value in the matrix
   :data:`m`, searching only the non-zero values, and storing them in :data:`imin` and :data:`jmin`.
   When there are several equal minimum elements then the first element found is returned.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. index::
   single: sparse matrices, compression

Compressed Format
=================

These routines calculate a compressed matrix from a coordinate representation.

.. function:: int gsl_spmatrix_csc (gsl_spmatrix * dest, const gsl_spmatrix * src)

   This function creates a sparse matrix in :ref:`compressed sparse column <sec_spmatrix-csc>`
   format from the input sparse matrix :data:`src` which must be in COO format. The
   compressed matrix is stored in :data:`dest`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`

.. function:: int gsl_spmatrix_csr (gsl_spmatrix * dest, const gsl_spmatrix * src)

   This function creates a sparse matrix in :ref:`compressed sparse row <sec_spmatrix-csr>`
   format from the input sparse matrix :data:`src` which must be in COO format. The
   compressed matrix is stored in :data:`dest`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`

.. function:: gsl_spmatrix * gsl_spmatrix_compress (const gsl_spmatrix * src, const int sptype)

   This function allocates a new sparse matrix, and stores :data:`src` into it using the
   format specified by :data:`sptype`. The input :data:`sptype` can be one of
   :macro:`GSL_SPMATRIX_COO`, :macro:`GSL_SPMATRIX_CSC`, or :macro:`GSL_SPMATRIX_CSR`.
   A pointer to the newly allocated matrix is returned, and must be freed by the caller
   when no longer needed.

.. index::
   single: sparse matrices, conversion

Conversion Between Sparse and Dense Matrices
============================================

The :type:`gsl_spmatrix` structure can be converted into the dense :type:`gsl_matrix`
format and vice versa with the following routines.

.. function:: int gsl_spmatrix_d2sp (gsl_spmatrix * S, const gsl_matrix * A)

   This function converts the dense matrix :data:`A` into sparse COO format
   and stores the result in :data:`S`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`

.. function:: int gsl_spmatrix_sp2d (gsl_matrix * A, const gsl_spmatrix * S)

   This function converts the sparse matrix :data:`S` into a dense matrix and
   stores the result in :data:`A`.

   Input matrix formats supported: :ref:`COO <sec_spmatrix-coo>`, :ref:`CSC <sec_spmatrix-csc>`, :ref:`CSR <sec_spmatrix-csr>`

.. index::
   single: sparse matrices, examples

Examples
========

The following example program builds a 5-by-4 sparse matrix
and prints it in coordinate, compressed column, and compressed
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
