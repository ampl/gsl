.. index::
   single: linear algebra
   single: solution of linear systems, Ax=b
   single: matrix factorization
   single: factorization of matrices

**************
Linear Algebra
**************

.. include:: include.rst

This chapter describes functions for solving linear systems.  The
library provides linear algebra operations which operate directly on
the :type:`gsl_vector` and :type:`gsl_matrix` objects.  These routines
use the standard algorithms from Golub & Van Loan's *Matrix
Computations* with Level-1 and Level-2 BLAS calls for efficiency.

The functions described in this chapter are declared in the header file
:file:`gsl_linalg.h`.

.. index:: LU decomposition

.. _sec_lu-decomposition:

LU Decomposition
================

A general :math:`M`-by-:math:`N` matrix :math:`A` has an :math:`LU` decomposition

.. math:: P A = L U

where :math:`P` is an :math:`M`-by-:math:`M` permutation matrix, :math:`L` is
:math:`M`-by-:math:`\min(M,N)` and :math:`U` is :math:`\min(M,N)`-by-:math:`N`.
For square matrices, :math:`L` is a lower unit triangular matrix and
:math:`U` is upper triangular. For :math:`M > N`, :math:`L` is a unit lower
trapezoidal matrix of size :math:`M`-by-:math:`N`. For :math:`M < N`,
:math:`U` is upper trapezoidal of size :math:`M`-by-:math:`N`.
For square matrices this decomposition can be used to convert the linear system
:math:`A x = b` into a pair of triangular systems (:math:`L y = P b`,
:math:`U x = y`), which can be solved by forward and back-substitution.
Note that the :math:`LU` decomposition is valid for singular matrices.

.. function:: int gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int * signum)
              int gsl_linalg_complex_LU_decomp (gsl_matrix_complex * A, gsl_permutation * p, int * signum)

   These functions factorize the matrix :data:`A` into the :math:`LU`
   decomposition :math:`PA = LU`.  On output the diagonal and upper
   triangular (or trapezoidal) part of the input matrix :data:`A` contain the matrix
   :math:`U`. The lower triangular (or trapezoidal) part of the input matrix (excluding the
   diagonal) contains :math:`L`.  The diagonal elements of :math:`L` are
   unity, and are not stored.

   The permutation matrix :math:`P` is encoded in the permutation
   :data:`p` on output. The :math:`j`-th column of the matrix :math:`P`
   is given by the :math:`k`-th column of the identity matrix, where
   :math:`k = p_j` the
   :math:`j`-th element of the permutation vector. The sign of the
   permutation is given by :data:`signum`. It has the value :math:`(-1)^n`,
   where :math:`n` is the number of interchanges in the permutation.

   The algorithm used in the decomposition is Gaussian Elimination with
   partial pivoting (Golub & Van Loan, *Matrix Computations*,
   Algorithm 3.4.1), combined with a recursive algorithm based on
   Level 3 BLAS (Peise and Bientinesi, 2016).

   The functions return :macro:`GSL_SUCCESS` for non-singular matrices.
   If the matrix is singular, the factorization is still completed, and
   the functions return an integer :math:`k \in [1, MIN(M,N)]` such that
   :math:`U(k,k) = 0`. In this case, :math:`U` is singular and the
   decomposition should not be used to solve linear systems.

.. index:: linear systems, solution of

.. function:: int gsl_linalg_LU_solve (const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)
              int gsl_linalg_complex_LU_solve (const gsl_matrix_complex * LU, const gsl_permutation * p, const gsl_vector_complex * b, gsl_vector_complex * x)

   These functions solve the square system :math:`A x = b` using the :math:`LU`
   decomposition of :math:`A` into (:data:`LU`, :data:`p`) given by
   :func:`gsl_linalg_LU_decomp` or :func:`gsl_linalg_complex_LU_decomp` as input.

.. function:: int gsl_linalg_LU_svx (const gsl_matrix * LU, const gsl_permutation * p, gsl_vector * x)
              int gsl_linalg_complex_LU_svx (const gsl_matrix_complex * LU, const gsl_permutation * p, gsl_vector_complex * x)

   These functions solve the square system :math:`A x = b` in-place using the
   precomputed :math:`LU` decomposition of :math:`A` into (:data:`LU`, :data:`p`). On input
   :data:`x` should contain the right-hand side :math:`b`, which is replaced
   by the solution on output.

.. index::
   single: refinement of solutions in linear systems
   single: iterative refinement of solutions in linear systems
   single: linear systems, refinement of solutions

.. function:: int gsl_linalg_LU_refine (const gsl_matrix * A, const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x, gsl_vector * work)
              int gsl_linalg_complex_LU_refine (const gsl_matrix_complex * A, const gsl_matrix_complex * LU, const gsl_permutation * p, const gsl_vector_complex * b, gsl_vector_complex * x, gsl_vector_complex * work)

   These functions apply an iterative improvement to :data:`x`, the solution
   of :math:`A x = b`, from the precomputed :math:`LU` decomposition of :math:`A` into
   (:data:`LU`, :data:`p`). Additional workspace of length :data:`N` is required in :data:`work`.

.. index::
   single: inverse of a matrix, by LU decomposition
   single: matrix inverse

.. function:: int gsl_linalg_LU_invert (const gsl_matrix * LU, const gsl_permutation * p, gsl_matrix * inverse)
              int gsl_linalg_complex_LU_invert (const gsl_matrix_complex * LU, const gsl_permutation * p, gsl_matrix_complex * inverse)

   These functions compute the inverse of a matrix :math:`A` from its
   :math:`LU` decomposition (:data:`LU`, :data:`p`), storing the result in the
   matrix :data:`inverse`.  The inverse is computed by computing the inverses
   :math:`U^{-1}`, :math:`L^{-1}` and finally forming the product
   :math:`A^{-1} = U^{-1} L^{-1} P`. Each step is based on Level 3 BLAS calls.
   
   It is preferable
   to avoid direct use of the inverse whenever possible, as the linear
   solver functions can obtain the same result more efficiently and
   reliably (consult any introductory textbook on numerical linear algebra
   for details).

.. function:: int gsl_linalg_LU_invx (gsl_matrix * LU, const gsl_permutation * p)
              int gsl_linalg_complex_LU_invx (gsl_matrix_complex * LU, const gsl_permutation * p)

   These functions compute the inverse of a matrix :math:`A` from its
   :math:`LU` decomposition (:data:`LU`, :data:`p`), storing the result in-place
   in the matrix :data:`LU`. The inverse is computed by computing the inverses
   :math:`U^{-1}`, :math:`L^{-1}` and finally forming the product
   :math:`A^{-1} = U^{-1} L^{-1} P`. Each step is based on Level 3 BLAS calls.

   It is preferable to avoid direct use of the inverse whenever possible, as the linear
   solver functions can obtain the same result more efficiently and
   reliably (consult any introductory textbook on numerical linear algebra
   for details).

.. index::
   single: determinant of a matrix, by LU decomposition
   single: matrix determinant

.. function:: double gsl_linalg_LU_det (gsl_matrix * LU, int signum)
              gsl_complex gsl_linalg_complex_LU_det (gsl_matrix_complex * LU, int signum)

   These functions compute the determinant of a matrix :math:`A` from its
   :math:`LU` decomposition, :data:`LU`. The determinant is computed as the
   product of the diagonal elements of :math:`U` and the sign of the row
   permutation :data:`signum`.

.. index:: logarithm of the determinant of a matrix

.. function:: double gsl_linalg_LU_lndet (gsl_matrix * LU)
              double gsl_linalg_complex_LU_lndet (gsl_matrix_complex * LU)

   These functions compute the logarithm of the absolute value of the
   determinant of a matrix :math:`A`, :math:`\ln|\det(A)|`, from its :math:`LU`
   decomposition, :data:`LU`.  This function may be useful if the direct
   computation of the determinant would overflow or underflow.

.. index:: sign of the determinant of a matrix

.. function:: int gsl_linalg_LU_sgndet (gsl_matrix * LU, int signum)
              gsl_complex gsl_linalg_complex_LU_sgndet (gsl_matrix_complex * LU, int signum)

   These functions compute the sign or phase factor of the determinant of a
   matrix :math:`A`, :math:`\det(A)/|\det(A)|`, from its :math:`LU` decomposition,
   :data:`LU`.

.. index:: QR decomposition

.. _linalg-qr:

QR Decomposition
================

A general rectangular :math:`M`-by-:math:`N` matrix :math:`A` has a
:math:`QR` decomposition into the product of a unitary
:math:`M`-by-:math:`M` square matrix :math:`Q` (where :math:`Q^{\dagger} Q = I`) and
an :math:`M`-by-:math:`N` right-triangular matrix :math:`R`,

.. math:: A = Q R

This decomposition can be used to convert the square linear system :math:`A x = b`
into the triangular system :math:`R x = Q^{\dagger} b`, which can be solved by
back-substitution. Another use of the :math:`QR` decomposition is to
compute an orthonormal basis for a set of vectors. The first :math:`N`
columns of :math:`Q` form an orthonormal basis for the range of :math:`A`,
:math:`ran(A)`, when :math:`A` has full column rank.

When :math:`M > N`, the bottom :math:`M - N` rows of :math:`R` are zero,
and so :math:`A` can be naturally partioned as

.. only:: not texinfo

   .. math:: A = \begin{pmatrix} Q_1 & Q_2 \end{pmatrix} \begin{pmatrix} R_1 \\ 0 \end{pmatrix} = Q_1 R_1

.. only:: texinfo

   ::

     A = [ Q_1 Q_2 ] [ R_1  ] = Q_1 R_1
                     [  0   ]

where :math:`R_1` is :math:`N`-by-:math:`N` upper triangular, :math:`Q_1` is :math:`M`-by-:math:`N`, and
:math:`Q_2` is :math:`M`-by-:math:`(M-N)`. :math:`Q_1 R_1` is sometimes called the *thin* or *reduced*
QR decomposition. The solution of the least squares problem :math:`\min_x || b - A x ||^2` when :math:`A`
has full rank is:

.. math:: x = R_1^{-1} c_1

where :math:`c_1` is the first :math:`N` elements of :math:`Q^{\dagger} b`. If :math:`A` is rank deficient,
see :ref:`linalg-qrpt` and :ref:`linalg-cod`.

GSL offers two interfaces for the :math:`QR` decomposition. The first proceeds by zeroing
out columns below the diagonal of :math:`A`, one column at a time using Householder transforms.
In this method, the factor :math:`Q` is represented as a product of Householder reflectors:

.. math:: Q = H_n \cdots H_2 H_1

where each :math:`H_i = I - \tau_i v_i v_i^{\dagger}` for a scalar :math:`\tau_i` and column vector
:math:`v_i`. In this method, functions which compute the full matrix :math:`Q` or apply
:math:`Q^{\dagger}` to a right hand side vector operate by applying the Householder matrices one
at a time using Level 2 BLAS.

The second interface is based on a Level 3 BLAS block recursive algorithm developed by
Elmroth and Gustavson. In this case, :math:`Q` is written in block form as

.. math:: Q = I - V T V^{\dagger}

where :math:`V` is an :math:`M`-by-:math:`N` matrix of the column vectors :math:`v_i`
and :math:`T` is an :math:`N`-by-:math:`N` upper triangular matrix, whose diagonal elements
are the :math:`\tau_i`. Computing the full :math:`T`, while requiring more flops than
the Level 2 approach, offers the advantage that all standard operations can take advantage
of cache efficient Level 3 BLAS operations, and so this method often performs faster
than the Level 2 approach. The functions for the recursive block algorithm have a
:code:`_r` suffix, and it is recommended to use this interface for performance
critical applications.

.. function:: int gsl_linalg_QR_decomp_r (gsl_matrix * A, gsl_matrix * T)
              int gsl_linalg_complex_QR_decomp_r (gsl_matrix_complex * A, gsl_matrix_complex * T)

   These functions factor the :math:`M`-by-:math:`N` matrix :data:`A` into
   the :math:`QR` decomposition :math:`A = Q R` using the recursive Level 3 BLAS
   algorithm of Elmroth and Gustavson.  On output the diagonal and
   upper triangular part of :data:`A` contain the matrix
   :math:`R`.  The :math:`N`-by-:math:`N`
   matrix :data:`T` stores the upper triangular factor appearing in :math:`Q`.
   The matrix :math:`Q` is given by :math:`Q = I - V T V^{\dagger}`, where the elements
   below the diagonal of :data:`A` contain the columns of :math:`V` on output.

   This algorithm requires :math:`M \ge N` and performs best for
   "tall-skinny" matrices, i.e. :math:`M \gg N`.

.. function:: int gsl_linalg_QR_solve_r (const gsl_matrix * QR, const gsl_matrix * T, const gsl_vector * b, gsl_vector * x)
              int gsl_linalg_complex_QR_solve_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T, const gsl_vector_complex * b, gsl_vector_complex * x)

   These functions solve the square system :math:`A x = b` using the :math:`QR`
   decomposition of :math:`A` held in (:data:`QR`, :data:`T`).
   The least-squares solution for rectangular systems can be found using
   :func:`gsl_linalg_QR_lssolve_r` or :func:`gsl_linalg_complex_QR_lssolve_r`.

.. function:: int gsl_linalg_QR_lssolve_r (const gsl_matrix * QR, const gsl_matrix * T, const gsl_vector * b, gsl_vector * x, gsl_vector * work)
              int gsl_linalg_complex_QR_lssolve_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T, const gsl_vector_complex * b, gsl_vector_complex * x, gsl_vector_complex * work)

   These functions find the least squares solution to the overdetermined
   system :math:`A x = b`, where the matrix :data:`A` has more rows than
   columns.  The least squares solution minimizes the Euclidean norm of the
   residual, :math:`||b - Ax||`. The routine requires as input 
   the :math:`QR` decomposition of :math:`A` into (:data:`QR`, :data:`T`) given by
   :func:`gsl_linalg_QR_decomp_r` or :func:`gsl_linalg_complex_QR_decomp_r`.

   The parameter :data:`x` is of length :math:`M`.
   The solution :math:`x` is returned in the first :math:`N` rows of :data:`x`,
   i.e. :math:`x =` :code:`x[0], x[1], ..., x[N-1]`. The last :math:`M - N` rows
   of :data:`x` contain a vector whose norm is equal to the residual norm
   :math:`|| b - A x ||`. This similar to the behavior of LAPACK DGELS.
   Additional workspace of length :math:`N` is required in :data:`work`.

.. function:: int gsl_linalg_QR_lssolvem_r (const gsl_matrix * QR, const gsl_matrix * T, const gsl_matrix * B, gsl_matrix * X, gsl_matrix * work)
              int gsl_linalg_complex_QR_lssolvem_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T, const gsl_matrix_complex * B, gsl_matrix_complex * X, gsl_matrix_complex * work)

   These functions find the least squares solutions to the overdetermined
   systems :math:`A x_k = b_k` where the matrix :data:`A` has more rows than
   columns.  The least squares solution minimizes the Euclidean norm of the
   residual, :math:`||b_k - A x_k||`.  The routine requires as input 
   the :math:`QR` decomposition of :math:`A` into (:data:`QR`, :data:`T`) given by
   :func:`gsl_linalg_QR_decomp_r` or :func:`gsl_linalg_complex_QR_decomp_r`.
   The right hand side :math:`b_k` is provided in column :math:`k`
   of the input :data:`B`, while the solution :math:`x_k` is stored in the
   first :math:`N` rows of column :math:`k` of the output :data:`X`.

   The parameters :data:`X` and :data:`B` are of size :math:`M`-by-:math:`nrhs`.
   The last :math:`M - N` rows of :data:`X` contain vectors whose norm is equal to the residual norm
   :math:`|| b_k - A x_k ||`. This similar to the behavior of LAPACK DGELS.
   Additional workspace of length :math:`N`-by-:math:`nrhs` is required in :data:`work`.

.. function:: int gsl_linalg_QR_QTvec_r (const gsl_matrix * QR, const gsl_matrix * T, gsl_vector * v, gsl_vector * work)
              int gsl_linalg_complex_QR_QHvec_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T, gsl_vector_complex * v, gsl_vector_complex * work)

   These functions apply the matrix :math:`Q^T` (or :math:`Q^{\dagger}`) encoded in the decomposition
   (:data:`QR`, :data:`T`) to the vector :data:`v`,
   storing the result :math:`Q^T v` (or :math:`Q^{\dagger} v`) in :data:`v`.  The matrix multiplication is carried
   out directly using the encoding of the Householder vectors without needing to form the full
   matrix :math:`Q`. Additional workspace of size :math:`N` is required in :data:`work`.

.. function:: int gsl_linalg_QR_QTmat_r (const gsl_matrix * QR, const gsl_matrix * T, gsl_matrix * B, gsl_matrix * work)
              int gsl_linalg_complex_QR_QHmat_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T, gsl_matrix_complex * B, gsl_matrix_complex * work)

   This function applies the matrix :math:`Q^T` (or :math:`Q^{\dagger}`) encoded in the decomposition
   (:data:`QR`, :data:`T`) to the :math:`M`-by-:math:`K` matrix :data:`B`,
   storing the result :math:`Q^T B` (or :math:`Q^{\dagger} B`) in :data:`B`.  The matrix multiplication is carried
   out directly using the encoding of the Householder vectors without needing to form the full
   matrix :math:`Q`.  Additional workspace of size :math:`N`-by-:math:`K` is required in :data:`work`.

.. function:: int gsl_linalg_QR_unpack_r (const gsl_matrix * QR, const gsl_matrix * T, gsl_matrix * Q, gsl_matrix * R)
              int gsl_linalg_complex_QR_unpack_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T, gsl_matrix_complex * Q, gsl_matrix_complex * R)

   These functions unpack the encoded :math:`QR` decomposition
   (:data:`QR`, :data:`T`) as output from :func:`gsl_linalg_QR_decomp_r`
   or :func:`gsl_linalg_complex_QR_decomp_r` into the matrices
   :data:`Q` and :data:`R`, where :data:`Q` is :math:`M`-by-:math:`M` and :data:`R` is :math:`N`-by-:math:`N`.
   Note that the full :math:`R` matrix is :math:`M`-by-:math:`N`, however the lower trapezoidal portion
   is zero, so only the upper triangular factor is stored.

.. function:: int gsl_linalg_QR_rcond (const gsl_matrix * QR, double * rcond, gsl_vector * work)

   This function estimates the reciprocal condition number (using the 1-norm) of the :math:`R` factor,
   stored in the upper triangle of :data:`QR`. The reciprocal condition number estimate, defined as
   :math:`1 / (||R||_1 \cdot ||R^{-1}||_1)`, is stored in :data:`rcond`.
   Additional workspace of size :math:`3 N` is required in :data:`work`.

Level 2 Interface
-----------------

The functions below are for the slower Level 2 interface to the QR decomposition. It is
recommended to use these functions only for :math:`M < N`, since the Level 3 interface
above performs much faster for :math:`M \ge N`.

.. function:: int gsl_linalg_QR_decomp (gsl_matrix * A, gsl_vector * tau)
              int gsl_linalg_complex_QR_decomp (gsl_matrix_complex * A, gsl_vector_complex * tau)

   These functions factor the :math:`M`-by-:math:`N` matrix :data:`A` into
   the :math:`QR` decomposition :math:`A = Q R`.  On output the diagonal and
   upper triangular part of the input matrix contain the matrix
   :math:`R`. The vector :data:`tau` and the columns of the lower triangular
   part of the matrix :data:`A` contain the Householder coefficients and
   Householder vectors which encode the orthogonal matrix :data:`Q`. The
   vector :data:`tau` must be of length :math:`N`. The matrix
   :math:`Q` is related to these components by the product of
   :math:`k=min(M,N)` reflector matrices, :math:`Q = H_k ... H_2 H_1`
   where :math:`H_i = I - \tau_i v_i v_i^{\dagger}` and :math:`v_i` is the
   Householder vector :math:`v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i))`.
   This is the same storage scheme as used by |lapack|.

   The algorithm used to perform the decomposition is Householder QR (Golub
   & Van Loan, "Matrix Computations", Algorithm 5.2.1).

.. function:: int gsl_linalg_QR_solve (const gsl_matrix * QR, const gsl_vector * tau, const gsl_vector * b, gsl_vector * x)
              int gsl_linalg_complex_QR_solve (const gsl_matrix_complex * QR, const gsl_vector_complex * tau, const gsl_vector_complex * b, gsl_vector_complex * x)

   These functions solve the square system :math:`A x = b` using the :math:`QR`
   decomposition of :math:`A` held in (:data:`QR`, :data:`tau`).
   The least-squares solution for rectangular systems can be found using :func:`gsl_linalg_QR_lssolve`.

.. function:: int gsl_linalg_QR_svx (const gsl_matrix * QR, const gsl_vector * tau, gsl_vector * x)
              int gsl_linalg_complex_QR_svx (const gsl_matrix_complex * QR, const gsl_vector_complex * tau, gsl_vector_complex * x)

   These functions solve the square system :math:`A x = b` in-place using
   the :math:`QR` decomposition of :math:`A` held in (:data:`QR`, :data:`tau`).
   On input :data:`x` should contain the right-hand side :math:`b`, which is replaced by the solution on output.

.. function:: int gsl_linalg_QR_lssolve (const gsl_matrix * QR, const gsl_vector * tau, const gsl_vector * b, gsl_vector * x, gsl_vector * residual)
              int gsl_linalg_complex_QR_lssolve (const gsl_matrix_complex * QR, const gsl_vector_complex * tau, const gsl_vector_complex * b, gsl_vector_complex * x, gsl_vector_complex * residual)

   These functions find the least squares solution to the overdetermined
   system :math:`A x = b` where the matrix :data:`A` has more rows than
   columns.  The least squares solution minimizes the Euclidean norm of the
   residual, :math:`||Ax - b||`.The routine requires as input 
   the :math:`QR` decomposition
   of :math:`A` into (:data:`QR`, :data:`tau`) given by
   :func:`gsl_linalg_QR_decomp` or :func:`gsl_linalg_complex_QR_decomp`.  The solution is returned in :data:`x`.  The
   residual is computed as a by-product and stored in :data:`residual`.

.. function:: int gsl_linalg_QR_QTvec (const gsl_matrix * QR, const gsl_vector * tau, gsl_vector * v)
              int gsl_linalg_complex_QR_QHvec (const gsl_matrix_complex * QR, const gsl_vector_complex * tau, gsl_vector_complex * v)

   These functions apply the matrix :math:`Q^T` (or :math:`Q^{\dagger}`) encoded in the decomposition
   (:data:`QR`, :data:`tau`) to the vector :data:`v`,
   storing the result :math:`Q^T v` (or :math:`Q^{\dagger} v`) in :data:`v`.  The matrix multiplication is carried
   out directly using the encoding of the Householder vectors without needing to form the full
   matrix :math:`Q^T` (or :math:`Q^{\dagger}`).

.. function:: int gsl_linalg_QR_Qvec (const gsl_matrix * QR, const gsl_vector * tau, gsl_vector * v)
              int gsl_linalg_complex_QR_Qvec (const gsl_matrix_complex * QR, const gsl_vector_complex * tau, gsl_vector_complex * v)

   These functions apply the matrix :math:`Q` encoded in the decomposition
   (:data:`QR`, :data:`tau`) to the vector :data:`v`, storing the result :math:`Q v`
   in :data:`v`.  The matrix multiplication is carried out directly using
   the encoding of the Householder vectors without needing to form the full
   matrix :math:`Q`.

.. function:: int gsl_linalg_QR_QTmat (const gsl_matrix * QR, const gsl_vector * tau, gsl_matrix * B)

   This function applies the matrix :math:`Q^T` encoded in the decomposition
   (:data:`QR`, :data:`tau`) to the :math:`M`-by-:math:`K` matrix :data:`B`,
   storing the result :math:`Q^T B` in :data:`B`.  The matrix multiplication is carried
   out directly using the encoding of the Householder vectors without needing to form the full
   matrix :math:`Q^T`.

.. function:: int gsl_linalg_QR_Rsolve (const gsl_matrix * QR, const gsl_vector * b, gsl_vector * x)

   This function solves the triangular system :math:`R x = b` for
   :data:`x`. It may be useful if the product :math:`b' = Q^T b` has already
   been computed using :func:`gsl_linalg_QR_QTvec`.

.. function:: int gsl_linalg_QR_Rsvx (const gsl_matrix * QR, gsl_vector * x)

   This function solves the triangular system :math:`R x = b` for :data:`x`
   in-place. On input :data:`x` should contain the right-hand side :math:`b`
   and is replaced by the solution on output. This function may be useful if
   the product :math:`b' = Q^T b` has already been computed using
   :func:`gsl_linalg_QR_QTvec`.

.. function:: int gsl_linalg_QR_unpack (const gsl_matrix * QR, const gsl_vector * tau, gsl_matrix * Q, gsl_matrix * R)

   This function unpacks the encoded :math:`QR` decomposition
   (:data:`QR`, :data:`tau`) into the matrices :data:`Q` and :data:`R`, where
   :data:`Q` is :math:`M`-by-:math:`M` and :data:`R` is :math:`M`-by-:math:`N`. 

.. function:: int gsl_linalg_QR_QRsolve (gsl_matrix * Q, gsl_matrix * R, const gsl_vector * b, gsl_vector * x)

   This function solves the system :math:`R x = Q^T b` for :data:`x`. It can
   be used when the :math:`QR` decomposition of a matrix is available in
   unpacked form as (:data:`Q`, :data:`R`).

.. function:: int gsl_linalg_QR_update (gsl_matrix * Q, gsl_matrix * R, gsl_vector * w, const gsl_vector * v)

   This function performs a rank-1 update :math:`w v^T` of the :math:`QR`
   decomposition (:data:`Q`, :data:`R`). The update is given by :math:`Q'R' = Q (R + w v^T)`
   where the output matrices :math:`Q` and :math:`R` are also
   orthogonal and right triangular. Note that :data:`w` is destroyed by the
   update.

.. function:: int gsl_linalg_R_solve (const gsl_matrix * R, const gsl_vector * b, gsl_vector * x)

   This function solves the triangular system :math:`R x = b` for the
   :math:`N`-by-:math:`N` matrix :data:`R`.

.. function:: int gsl_linalg_R_svx (const gsl_matrix * R, gsl_vector * x)

   This function solves the triangular system :math:`R x = b` in-place. On
   input :data:`x` should contain the right-hand side :math:`b`, which is
   replaced by the solution on output.

Triangle on Top of Rectangle
----------------------------

This section provides routines for computing the :math:`QR` decomposition of the
specialized matrix

.. only:: not texinfo

   .. math:: \begin{pmatrix} U \\ A \end{pmatrix} = Q R

.. only:: texinfo

   ::

     [ U ] = Q R
     [ A ]

where :math:`U` is an :math:`N`-by-:math:`N` upper triangular matrix, and :math:`A` is
an :math:`M`-by-:math:`N` dense matrix. This type of matrix arises, for example,
in the sequential TSQR algorithm. The Elmroth and Gustavson algorithm is used to
efficiently factor this matrix. Due to the upper triangular factor, the :math:`Q`
matrix takes the form

.. math:: Q = I - V T V^T

with

.. only:: not texinfo

   .. math:: V = \begin{pmatrix} I \\ Y \end{pmatrix}

.. only:: texinfo

   ::

     V = [ I ]
         [ Y ]

and :math:`Y` is dense and of the same dimensions as :math:`A`.

.. function:: int gsl_linalg_QR_UR_decomp (gsl_matrix * U, gsl_matrix * A, gsl_matrix * T)

   This function computes the :math:`QR` decomposition of the matrix :math:`(U ; A)`, where
   :math:`U` is :math:`N`-by-:math:`N` upper triangular and :math:`A` is :math:`M`-by-:math:`N`
   dense. On output, :math:`U` is replaced by the :math:`R` factor, and :math:`A` is replaced
   by :math:`Y`.  The :math:`N`-by-:math:`N` upper triangular block reflector is
   stored in :data:`T` on output.

.. function:: int gsl_linalg_QR_UR_lssolve (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T, const gsl_vector * b, gsl_vector * x, gsl_vector * work)

   This function finds the least squares solution to the overdetermined
   system,
   
   .. math:: \min_x \left| \left| b - \begin{pmatrix} U \\ A \end{pmatrix} x \right| \right|^2
      
   where :math:`U` is a :math:`N`-by-:math:`N` upper triangular matrix, and
   :math:`A` is a :math:`M`-by-:math:`N` dense matrix.
   The routine requires as input the :math:`QR` decomposition
   of :math:`(U; A)` into (:data:`R`, :data:`Y`, :data:`T`) given by
   :func:`gsl_linalg_QR_UR_decomp`.
   The parameter :data:`x` is of length :math:`N+M`.
   The solution :math:`x` is returned in the first :math:`N` rows of :data:`x`,
   i.e. :math:`x =` :code:`x[0], x[1], ..., x[N-1]`. The last :math:`M` rows
   of :data:`x` contain a vector whose norm is equal to the residual norm
   :math:`|| b - (U; A) x ||`. This similar to the behavior of LAPACK DGELS.
   Additional workspace of length :math:`N` is required in :data:`work`.

.. function:: int gsl_linalg_QR_UR_lssvx (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T, gsl_vector * x, gsl_vector * work)

   This function finds the least squares solution to the overdetermined
   system,
   
   .. math:: \min_x \left| \left| b - \begin{pmatrix} U \\ A \end{pmatrix} x \right| \right|^2
      
   in-place, where :math:`U` is a :math:`N`-by-:math:`N` upper triangular matrix, and
   :math:`A` is a :math:`M`-by-:math:`N` dense matrix.
   The routine requires as input the :math:`QR` decomposition
   of :math:`(U; A)` into (:data:`R`, :data:`Y`, :data:`T`) given by
   :func:`gsl_linalg_QR_UR_decomp`.
   The parameter :data:`x` is of length :math:`N+M` and contains the right
   hand side vector :math:`b` on input.
   The solution :math:`x` is returned in the first :math:`N` rows of :data:`x`,
   i.e. :math:`x =` :code:`x[0], x[1], ..., x[N-1]`. The last :math:`M` rows
   of :data:`x` contain a vector whose norm is equal to the residual norm
   :math:`|| b - (U; A) x ||`. This similar to the behavior of LAPACK DGELS.
   Additional workspace of length :math:`N` is required in :data:`work`.

.. function:: int gsl_linalg_QR_UR_QTvec(const gsl_matrix * Y, const gsl_matrix * T, gsl_vector * b, gsl_vector * work)

   This function computes :math:`Q^T b` using the decomposition
   (:data:`Y`, :data:`T`) previously computed by :func:`gsl_linalg_QR_UR_decomp`.
   On input, :data:`b` contains the length :math:`N+M` vector :math:`b`, and on output it will contain
   :math:`Q^T b`. Additional workspace of length :math:`N` is required in :data:`work`.

Triangle on Top of Triangle
---------------------------

This section provides routines for computing the :math:`QR` decomposition of the
specialized matrix

.. only:: not texinfo

   .. math:: \begin{pmatrix} U_1 \\ U_2 \end{pmatrix} = Q R

.. only:: texinfo

   ::

     [ U_1 ] = Q R
     [ U_2 ]

where :math:`U_1,U_2` are :math:`N`-by-:math:`N` upper triangular matrices.
The Elmroth and Gustavson algorithm is used to efficiently factor this matrix.
The :math:`Q` matrix takes the form

.. math:: Q = I - V T V^T

with

.. only:: not texinfo

   .. math:: V = \begin{pmatrix} I \\ Y \end{pmatrix}

.. only:: texinfo

   ::

     V = [ I ]
         [ Y ]

and :math:`Y` is :math:`N`-by-:math:`N` upper triangular.

.. function:: int gsl_linalg_QR_UU_decomp (gsl_matrix * U1, gsl_matrix * U2, gsl_matrix * T)

   This function computes the :math:`QR` decomposition of the matrix :math:`(U_1 ; U_2)`, where
   :math:`U_1,U_2` are :math:`N`-by-:math:`N` upper triangular. On output, :data:`U1`
   is replaced by the :math:`R` factor, and :data:`U2` is replaced by :math:`Y`. The
   :math:`N`-by-:math:`N` upper triangular block reflector is stored in :data:`T` on output.

.. function:: int gsl_linalg_QR_UU_lssolve (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T, const gsl_vector * b, gsl_vector * x, gsl_vector * work)

   This function finds the least squares solution to the overdetermined
   system,
   
   .. math:: \min_x \left| \left| b - \begin{pmatrix} U_1 \\ U_2 \end{pmatrix} x \right| \right|^2
      
   where :math:`U_1,U_2` are :math:`N`-by-:math:`N` upper triangular matrices.
   The routine requires as input the :math:`QR` decomposition
   of :math:`(U_1; U_2)` into (:data:`R`, :data:`Y`, :data:`T`) given by
   :func:`gsl_linalg_QR_UU_decomp`.
   The parameter :data:`x` is of length :math:`2N`.
   The solution :math:`x` is returned in the first :math:`N` rows of :data:`x`,
   i.e. :math:`x =` :code:`x[0], x[1], ..., x[N-1]`. The last :math:`N` rows
   of :data:`x` contain a vector whose norm is equal to the residual norm
   :math:`|| b - (U_1; U_2) x ||`. This similar to the behavior of LAPACK DGELS.
   Additional workspace of length :math:`N` is required in :data:`work`.

.. function:: int gsl_linalg_QR_UU_lssvx (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T, gsl_vector * x, gsl_vector * work)

   This function finds the least squares solution to the overdetermined
   system,
   
   .. math:: \min_x \left| \left| b - \begin{pmatrix} U_1 \\ U_2 \end{pmatrix} x \right| \right|^2
      
   in-place, where :math:`U_1,U_2` are :math:`N`-by-:math:`N` upper triangular matrices.
   The routine requires as input the :math:`QR` decomposition
   of :math:`(U_1; U_2)` into (:data:`R`, :data:`Y`, :data:`T`) given by
   :func:`gsl_linalg_QR_UU_decomp`.
   The parameter :data:`x` is of length :math:`2N` and contains the right hand
   side vector :math:`b` on input.
   The solution :math:`x` is returned in the first :math:`N` rows of :data:`x`,
   i.e. :math:`x =` :code:`x[0], x[1], ..., x[N-1]`. The last :math:`N` rows
   of :data:`x` contain a vector whose norm is equal to the residual norm
   :math:`|| b - (U_1; U_2) x ||`. This similar to the behavior of LAPACK DGELS.
   Additional workspace of length :math:`N` is required in :data:`work`.

.. function:: int gsl_linalg_QR_UU_QTvec (const gsl_matrix * Y, const gsl_matrix * T, gsl_vector * b, gsl_vector * work)

   This function computes :math:`Q^T b` using the decomposition
   (:data:`Y`, :data:`T`) previously computed by :func:`gsl_linalg_QR_UU_decomp`.
   On input, :data:`b` contains the vector :math:`b`, and on output it will contain
   :math:`Q^T b`. Additional workspace of length :math:`N` is required in :data:`work`.

Triangle on Top of Trapezoidal
------------------------------

This section provides routines for computing the :math:`QR` decomposition of the
specialized matrix

.. only:: not texinfo

   .. math:: \begin{pmatrix} U \\ A \end{pmatrix} = Q R

.. only:: texinfo

   ::

     [ U ] = Q R
     [ A ]

where :math:`U` is an :math:`N`-by-:math:`N` upper triangular matrix, and :math:`A` is
an :math:`M`-by-:math:`N` upper trapezoidal matrix with :math:`M \ge N`. :math:`A` has
the structure,

.. only:: not texinfo

   .. math:: A = \begin{pmatrix} A_d \\ A_u \end{pmatrix}

.. only:: texinfo

   ::

     A = [ A_d ]
         [ A_u ]

where :math:`A_d` is :math:`(M-N)`-by-:math:`N` dense, and :math:`A_u` is
:math:`N`-by-:math:`N` upper triangular.
The Elmroth and Gustavson algorithm is used to efficiently factor this matrix.
The :math:`Q` matrix takes the form

.. math:: Q = I - V T V^T

with

.. only:: not texinfo

   .. math:: V = \begin{pmatrix} I \\ Y \end{pmatrix}

.. only:: texinfo

   ::

     V = [ I ]
         [ Y ]

and :math:`Y` is upper trapezoidal and of the same dimensions as :math:`A`.

.. function:: int gsl_linalg_QR_UZ_decomp (gsl_matrix * U, gsl_matrix * A, gsl_matrix * T)

   This function computes the :math:`QR` decomposition of the matrix :math:`(U ; A)`, where
   :math:`U` is :math:`N`-by-:math:`N` upper triangular and :math:`A` is :math:`M`-by-:math:`N`
   upper trapezoidal. On output, :math:`U` is replaced by the :math:`R` factor, and :math:`A`
   is replaced by :math:`Y`. The :math:`N`-by-:math:`N` upper triangular block reflector is
   stored in :data:`T` on output.

Triangle on Top of Diagonal
---------------------------

This section provides routines for computing the :math:`QR` decomposition of the
specialized matrix

.. only:: not texinfo

   .. math:: \begin{pmatrix} U \\ D \end{pmatrix} = Q R

.. only:: texinfo

   ::

     [ U ] = Q R
     [ D ]

where :math:`U` is an :math:`N`-by-:math:`N` upper triangular matrix and
:math:`D` is an :math:`N`-by-:math:`N` diagonal matrix. This type of matrix
arises in regularized least squares problems.
The Elmroth and Gustavson algorithm is used to efficiently factor this matrix.
The :math:`Q` matrix takes the form

.. math:: Q = I - V T V^T

with

.. only:: not texinfo

   .. math:: V = \begin{pmatrix} I \\ Y \end{pmatrix}

.. only:: texinfo

   ::

     V = [ I ]
         [ Y ]

and :math:`Y` is :math:`N`-by-:math:`N` upper triangular.

.. function:: int gsl_linalg_QR_UD_decomp (gsl_matrix * U, const gsl_vector * D, gsl_matrix * Y, gsl_matrix * T)

   This function computes the :math:`QR` decomposition of the matrix :math:`(U ; D)`, where
   :math:`U` is :math:`N`-by-:math:`N` upper triangular and :math:`D` is
   :math:`N`-by-:math:`N` diagonal. On output, :data:`U`
   is replaced by the :math:`R` factor and :math:`Y` is stored in :data:`Y`. The
   :math:`N`-by-:math:`N` upper triangular block reflector is stored in :data:`T` on output.

.. function:: int gsl_linalg_QR_UD_lssolve (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T, const gsl_vector * b, gsl_vector * x, gsl_vector * work)

   This function finds the least squares solution to the overdetermined
   system,
   
   .. math:: \min_x \left| \left| b - \begin{pmatrix} U \\ D \end{pmatrix} x \right| \right|^2
      
   where :math:`U` is :math:`N`-by-:math:`N` upper triangular and :math:`D` is
   :math:`N`-by-:math:`N` diagonal.  The routine requires as input 
   the :math:`QR` decomposition of :math:`(U; D)` into (:data:`R`, :data:`Y`, :data:`T`)
   given by :func:`gsl_linalg_QR_UD_decomp`.
   The parameter :data:`x` is of length :math:`2N`.
   The solution :math:`x` is returned in the first :math:`N` rows of :data:`x`,
   i.e. :math:`x =` :code:`x[0], x[1], ..., x[N-1]`. The last :math:`N` rows
   of :data:`x` contain a vector whose norm is equal to the residual norm
   :math:`|| b - (U; D) x ||`. This similar to the behavior of LAPACK DGELS.
   Additional workspace of length :math:`N` is required in :data:`work`.

.. function:: int gsl_linalg_QR_UD_lssvx (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T, gsl_vector * x, gsl_vector * work)

   This function finds the least squares solution to the overdetermined
   system,
   
   .. math:: \min_x \left| \left| b - \begin{pmatrix} U \\ D \end{pmatrix} x \right| \right|^2
      
   in-place, where :math:`U` is :math:`N`-by-:math:`N` upper triangular and :math:`D` is
   :math:`N`-by-:math:`N` diagonal.  The routine requires as input 
   the :math:`QR` decomposition of :math:`(U; D)` into (:data:`R`, :data:`Y`, :data:`T`)
   given by :func:`gsl_linalg_QR_UD_decomp`.
   The parameter :data:`x` is of length :math:`2N` and contains the right hand side
   vector :math:`b` on input.
   The solution :math:`x` is returned in the first :math:`N` rows of :data:`x`,
   i.e. :math:`x =` :code:`x[0], x[1], ..., x[N-1]`. The last :math:`N` rows
   of :data:`x` contain a vector whose norm is equal to the residual norm
   :math:`|| b - (U; D) x ||`. This similar to the behavior of LAPACK DGELS.
   Additional workspace of length :math:`N` is required in :data:`work`.

.. function:: int gsl_linalg_QR_UD_QTvec (const gsl_matrix * Y, const gsl_matrix * T, gsl_vector * b, gsl_vector * work)

   This function computes :math:`Q^T b` using the decomposition
   (:data:`Y`, :data:`T`) previously computed by :func:`gsl_linalg_QR_UD_decomp`.
   On input, :data:`b` contains the vector :math:`b`, and on output it will contain
   :math:`Q^T b`. Additional workspace of length :math:`N` is required in :data:`work`.

.. index:: QR decomposition with column pivoting

.. _linalg-qrpt:

QR Decomposition with Column Pivoting
=====================================

The :math:`QR` decomposition of an :math:`M`-by-:math:`N` matrix :math:`A`
can be extended to the rank deficient case by introducing a column permutation :math:`P`,

.. math:: A P = Q R

The first :math:`r` columns of :math:`Q` form an orthonormal basis
for the range of :math:`A` for a matrix with column rank :math:`r`.  This
decomposition can also be used to convert the square linear system :math:`A x = b`
into the triangular system :math:`R y = Q^T b, x = P y`, which can be
solved by back-substitution and permutation.  We denote the :math:`QR`
decomposition with column pivoting by :math:`QRP^T` since :math:`A = Q R P^T`.
When :math:`A` is rank deficient with :math:`r = {\rm rank}(A)`, the matrix
:math:`R` can be partitioned as

.. only:: not texinfo

   .. math::

      R = \left(
      \begin{matrix}
        R_{11} & R_{12} \\
        0 & R_{22}
      \end{matrix}
      \right) \approx
      \left(
      \begin{matrix}
        R_{11} & R_{12} \\
        0 & 0
      \end{matrix}
      \right)

.. only:: texinfo

   ::

      R = [ R11 R12 ] =~ [ R11 R12 ]
          [  0  R22 ]    [  0   0  ]

where :math:`R_{11}` is :math:`r`-by-:math:`r` and nonsingular. In this case,
a *basic* least squares solution for the overdetermined system :math:`A x = b`
can be obtained as

.. only:: not texinfo

   .. math::

      x = P \left(
      \begin{matrix}
        R_{11}^{-1} c_1 \\
        0
      \end{matrix}
      \right)

.. only:: texinfo

   ::

      x = P [ R11^-1 c1 ]
            [     0     ]

where :math:`c_1` consists of the first :math:`r` elements of :math:`Q^T b`.
This basic solution is not guaranteed to be the minimum norm solution unless
:math:`R_{12} = 0` (see :ref:`Complete Orthogonal Decomposition <linalg-cod>`).

.. function:: int gsl_linalg_QRPT_decomp (gsl_matrix * A, gsl_vector * tau, gsl_permutation * p, int * signum, gsl_vector * norm)

   This function factorizes the :math:`M`-by-:math:`N` matrix :data:`A` into
   the :math:`QRP^T` decomposition :math:`A = Q R P^T`.  On output the
   diagonal and upper triangular part of the input matrix contain the
   matrix :math:`R`. The permutation matrix :math:`P` is stored in the
   permutation :data:`p`.  The sign of the permutation is given by
   :data:`signum`. It has the value :math:`(-1)^n`, where :math:`n` is the
   number of interchanges in the permutation. The vector :data:`tau` and the
   columns of the lower triangular part of the matrix :data:`A` contain the
   Householder coefficients and vectors which encode the orthogonal matrix
   :data:`Q`.  The vector :data:`tau` must be of length :math:`k=\min(M,N)`. The
   matrix :math:`Q` is related to these components by, :math:`Q = Q_k ... Q_2 Q_1`
   where :math:`Q_i = I - \tau_i v_i v_i^T` and :math:`v_i` is the
   Householder vector
   
   .. math:: v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i))

   This is the same storage scheme
   as used by |lapack|.  The vector :data:`norm` is a workspace of length
   :data:`N` used for column pivoting.

   The algorithm used to perform the decomposition is Householder QR with
   column pivoting (Golub & Van Loan, "Matrix Computations", Algorithm
   5.4.1).

.. function:: int gsl_linalg_QRPT_decomp2 (const gsl_matrix * A, gsl_matrix * q, gsl_matrix * r, gsl_vector * tau, gsl_permutation * p, int * signum, gsl_vector * norm)

   This function factorizes the matrix :data:`A` into the decomposition
   :math:`A = Q R P^T` without modifying :data:`A` itself and storing the
   output in the separate matrices :data:`q` and :data:`r`.

.. function:: int gsl_linalg_QRPT_solve (const gsl_matrix * QR, const gsl_vector * tau, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)

   This function solves the square system :math:`A x = b` using the :math:`QRP^T`
   decomposition of :math:`A` held in (:data:`QR`, :data:`tau`, :data:`p`) which must 
   have been computed previously by :func:`gsl_linalg_QRPT_decomp`.

.. function:: int gsl_linalg_QRPT_svx (const gsl_matrix * QR, const gsl_vector * tau, const gsl_permutation * p, gsl_vector * x)

   This function solves the square system :math:`A x = b` in-place using the
   :math:`QRP^T` decomposition of :math:`A` held in
   (:data:`QR`, :data:`tau`, :data:`p`). On input :data:`x` should contain the
   right-hand side :math:`b`, which is replaced by the solution on output.

.. function:: int gsl_linalg_QRPT_lssolve (const gsl_matrix * QR, const gsl_vector * tau, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x, gsl_vector * residual)

   This function finds the least squares solution to the overdetermined
   system :math:`A x = b` where the matrix :data:`A` has more rows than
   columns and is assumed to have full rank. The least squares solution minimizes
   the Euclidean norm of the residual, :math:`||b - A x||`. The routine requires as input 
   the :math:`QR` decomposition of :math:`A` into (:data:`QR`, :data:`tau`, :data:`p`) given by
   :func:`gsl_linalg_QRPT_decomp`.  The solution is returned in :data:`x`.  The
   residual is computed as a by-product and stored in :data:`residual`. For rank
   deficient matrices, :func:`gsl_linalg_QRPT_lssolve2` should be used instead.

.. function:: int gsl_linalg_QRPT_lssolve2 (const gsl_matrix * QR, const gsl_vector * tau, const gsl_permutation * p, const gsl_vector * b, const size_t rank, gsl_vector * x, gsl_vector * residual)

   This function finds the least squares solution to the overdetermined
   system :math:`A x = b` where the matrix :data:`A` has more rows than
   columns and has rank given by the input :data:`rank`. If the user does not
   know the rank of :math:`A`, the routine :func:`gsl_linalg_QRPT_rank` can be
   called to estimate it. The least squares solution is
   the so-called "basic" solution discussed above and may not be the minimum
   norm solution. The routine requires as input 
   the :math:`QR` decomposition of :math:`A` into (:data:`QR`, :data:`tau`, :data:`p`) given by
   :func:`gsl_linalg_QRPT_decomp`.  The solution is returned in :data:`x`.  The
   residual is computed as a by-product and stored in :data:`residual`.

.. function:: int gsl_linalg_QRPT_QRsolve (const gsl_matrix * Q, const gsl_matrix * R, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)

   This function solves the square system :math:`R P^T x = Q^T b` for
   :data:`x`. It can be used when the :math:`QR` decomposition of a matrix is
   available in unpacked form as (:data:`Q`, :data:`R`).

.. function:: int gsl_linalg_QRPT_update (gsl_matrix * Q, gsl_matrix * R, const gsl_permutation * p, gsl_vector * w, const gsl_vector * v)

   This function performs a rank-1 update :math:`w v^T` of the :math:`QRP^T`
   decomposition (:data:`Q`, :data:`R`, :data:`p`). The update is given by
   :math:`Q'R' = Q (R + w v^T P)` where the output matrices :math:`Q'` and
   :math:`R'` are also orthogonal and right triangular. Note that :data:`w` is
   destroyed by the update. The permutation :data:`p` is not changed.

.. function:: int gsl_linalg_QRPT_Rsolve (const gsl_matrix * QR, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)

   This function solves the triangular system :math:`R P^T x = b` for the
   :math:`N`-by-:math:`N` matrix :math:`R` contained in :data:`QR`.

.. function:: int gsl_linalg_QRPT_Rsvx (const gsl_matrix * QR, const gsl_permutation * p, gsl_vector * x)

   This function solves the triangular system :math:`R P^T x = b` in-place
   for the :math:`N`-by-:math:`N` matrix :math:`R` contained in :data:`QR`. On
   input :data:`x` should contain the right-hand side :math:`b`, which is
   replaced by the solution on output.

.. function:: size_t gsl_linalg_QRPT_rank (const gsl_matrix * QR, const double tol)

   This function estimates the rank of the triangular matrix :math:`R` contained in :data:`QR`.
   The algorithm simply counts the number of diagonal elements of :math:`R` whose absolute value
   is greater than the specified tolerance :data:`tol`. If the input :data:`tol` is negative,
   a default value of :math:`20 (M + N) eps(max(|diag(R)|))` is used.

.. function:: int gsl_linalg_QRPT_rcond (const gsl_matrix * QR, double * rcond, gsl_vector * work)

   This function estimates the reciprocal condition number (using the 1-norm) of the :math:`R` factor,
   stored in the upper triangle of :data:`QR`. The reciprocal condition number estimate, defined as
   :math:`1 / (||R||_1 \cdot ||R^{-1}||_1)`, is stored in :data:`rcond`.
   Additional workspace of size :math:`3 N` is required in :data:`work`.

.. index:: LQ decomposition

LQ Decomposition
================

A general rectangular :math:`M`-by-:math:`N` matrix :math:`A` has a
:math:`LQ` decomposition into the product of a lower trapezoidal
:math:`M`-by-:math:`N` matrix :math:`L` and an orthogonal
:math:`N`-by-:math:`N` square matrix :math:`Q`:

.. math:: A = L Q

If :math:`M \le N`, then :math:`L` can be written as :math:`L = (L_1 \quad 0)` where
:math:`L_1` is :math:`M`-by-:math:`M` lower triangular,
and

.. math:: A =
          \begin{pmatrix}
            L_1 & 0
          \end{pmatrix}
          \begin{pmatrix}
            Q_1 \\ Q_2
          \end{pmatrix} = L_1 Q_1

where :math:`Q_1` consists of the first :math:`M` rows of :math:`Q`, and :math:`Q_2`
contains the remaining :math:`N - M` rows. The :math:`LQ` factorization of :math:`A`
is essentially the same as the :ref:`QR factorization <linalg-qr>` of :math:`A^T`.

The :math:`LQ` factorization may be used to find the minimum norm solution of
an underdetermined system of equations :math:`A x = b`, where :math:`A` is
:math:`M`-by-:math:`N` and :math:`M \le N`. The solution is

.. math:: x = Q^T
              \begin{pmatrix}
                L_1^{-1} b \\
                0
              \end{pmatrix}

.. function:: int gsl_linalg_LQ_decomp (gsl_matrix * A, gsl_vector * tau)

   This function factorizes the :math:`M`-by-:math:`N` matrix :data:`A` into
   the :math:`LQ` decomposition :math:`A = L Q`.  On output the diagonal and
   lower trapezoidal part of the input matrix contain the matrix
   :math:`L`. The vector :data:`tau` and the elements above the diagonal
   of the matrix :data:`A` contain the Householder coefficients and
   Householder vectors which encode the orthogonal matrix :data:`Q`.  The
   vector :data:`tau` must be of length :math:`k=\min(M,N)`. The matrix
   :math:`Q` is related to these components by, :math:`Q = Q_k ... Q_2 Q_1`
   where :math:`Q_i = I - \tau_i v_i v_i^T` and :math:`v_i` is the
   Householder vector :math:`v_i = (0,...,1,A(i,i+1),A(i,i+2),...,A(i,N))`.
   This is the same storage scheme as used by |lapack|.

.. function:: int gsl_linalg_LQ_lssolve (const gsl_matrix * LQ, const gsl_vector * tau, const gsl_vector * b, gsl_vector * x, gsl_vector * residual)

   This function finds the minimum norm least squares solution to the underdetermined system :math:`A x = b`,
   where the :math:`M`-by-:math:`N` matrix :data:`A` has :math:`M \le N`.
   The routine requires as input the :math:`LQ` decomposition of :math:`A` into (:data:`LQ`, :data:`tau`)
   given by :func:`gsl_linalg_LQ_decomp`.  The solution is returned in :data:`x`.  The
   residual, :math:`b - Ax`, is computed as a by-product and stored in :data:`residual`.

.. function:: int gsl_linalg_LQ_unpack (const gsl_matrix * LQ, const gsl_vector * tau, gsl_matrix * Q, gsl_matrix * L)

   This function unpacks the encoded :math:`LQ` decomposition
   (:data:`LQ`, :data:`tau`) into the matrices :data:`Q` and :data:`L`, where
   :data:`Q` is :math:`N`-by-:math:`N` and :data:`L` is :math:`M`-by-:math:`N`. 

.. function:: int gsl_linalg_LQ_QTvec (const gsl_matrix * LQ, const gsl_vector * tau, gsl_vector * v)

   This function applies :math:`Q^T` to the vector :data:`v`, storing the result :math:`Q^T v` in
   :data:`v` on output.

.. index:: QL decomposition

QL Decomposition
================

A general rectangular :math:`M`-by-:math:`N` matrix :math:`A` has a
:math:`QL` decomposition into the product of an orthogonal
:math:`M`-by-:math:`M` square matrix :math:`Q` (where :math:`Q^T Q = I`) and
an :math:`M`-by-:math:`N` left-triangular matrix :math:`L`.

When :math:`M \ge N`, the decomposition is given by

.. only:: not texinfo

   .. math:: A = Q \begin{pmatrix} 0 \\ L_1 \end{pmatrix}

.. only:: texinfo

   ::

     A = [ 0 ; L_1 ]

where :math:`L_1` is :math:`N`-by-:math:`N` lower triangular. When
:math:`M \le N`, the decomposition is given by

.. only:: not texinfo

   .. math:: A = Q \begin{pmatrix} L_1 & L_2 \end{pmatrix}

.. only:: texinfo

   ::

     A = [ L_1 L_2 ]

where :math:`L_1` is a dense :math:`M`-by-:math:`N-M` matrix and
:math:`L_2` is a lower triangular :math:`M`-by-:math:`M` matrix.

.. function:: int gsl_linalg_QL_decomp (gsl_matrix * A, gsl_vector * tau)

   This function factorizes the :math:`M`-by-:math:`N` matrix :data:`A` into
   the :math:`QL` decomposition :math:`A = Q L`.
   The vector :data:`tau` must be of length :math:`N` and contains the Householder
   coefficients on output.
   The matrix :math:`Q` is stored in packed form in :data:`A` on output, using the
   same storage scheme as |lapack|.

.. function:: int gsl_linalg_QL_unpack (const gsl_matrix * QL, const gsl_vector * tau, gsl_matrix * Q, gsl_matrix * L)

   This function unpacks the encoded :math:`QL` decomposition
   (:data:`QL`, :data:`tau`) into the matrices :data:`Q` and :data:`L`, where
   :data:`Q` is :math:`M`-by-:math:`M` and :data:`L` is :math:`M`-by-:math:`N`. 

.. index:: complete orthogonal decomposition

.. _linalg-cod:

Complete Orthogonal Decomposition
=================================

The complete orthogonal decomposition of a :math:`M`-by-:math:`N` matrix
:math:`A` is a generalization of the QR decomposition with column pivoting, given by

.. only:: not texinfo

   .. math::

      A P = Q
      \left(
      \begin{matrix}
        R_{11} & 0 \\
        0 & 0
      \end{matrix}
      \right) Z^T

.. only:: texinfo

   ::

      A P = Q [ R11 0 ] Z^T
              [  0  0 ]

where :math:`P` is a :math:`N`-by-:math:`N` permutation matrix,
:math:`Q` is :math:`M`-by-:math:`M` orthogonal, :math:`R_{11}` is
:math:`r`-by-:math:`r` upper triangular, with :math:`r = {\rm rank}(A)`,
and :math:`Z` is :math:`N`-by-:math:`N` orthogonal. If :math:`A`
has full rank, then :math:`R_{11} = R`, :math:`Z = I` and this reduces to the
QR decomposition with column pivoting.

For a rank deficient least squares problem, :math:`\min_x{|| b - Ax||^2}`, the solution vector
:math:`x` is not unique. However if we further require that :math:`||x||^2` is minimized,
then the complete orthogonal decomposition gives us the ability to compute
the unique minimum norm solution, which is given by

.. only:: not texinfo

   .. math::

      x = P Z
      \left(
      \begin{matrix}
        R_{11}^{-1} c_1 \\
        0
      \end{matrix}
      \right)

.. only:: texinfo

   ::

      x = P Z [ R11^-1 c1 ]
              [     0     ]

and the vector :math:`c_1` is the first :math:`r` elements of :math:`Q^T b`.

The COD also enables a straightforward solution of regularized least squares problems
in Tikhonov standard form, written as

.. math:: \min_x ||b - A x||^2 + \lambda^2 ||x||^2

where :math:`\lambda > 0` is a regularization parameter which represents a tradeoff between
minimizing the residual norm :math:`||b - A x||` and the solution norm :math:`||x||`. For this system,
the solution is given by

.. only:: not texinfo

   .. math::

      x = P Z
      \left(
      \begin{matrix}
        y_1 \\
        0
      \end{matrix}
      \right)

.. only:: texinfo

   ::

      x = P Z [ y1 ]
              [ 0  ]

where :math:`y_1` is a vector of length :math:`r` which is found by solving

.. only:: not texinfo

   .. math::

      \left(
      \begin{matrix}
        R_{11} \\
        \lambda I_r
      \end{matrix}
      \right) y_1 =
      \left(
      \begin{matrix}
        c_1 \\
        0
      \end{matrix}
      \right)

.. only:: texinfo

   ::

      [     R11     ] y_1 = [ c_1 ]
      [ \lambda I_r ]       [  0  ]

and :math:`c_1` is defined above. The equation above can be solved efficiently for different
values of :math:`\lambda` using QR factorizations of the left hand side matrix.

.. function:: int gsl_linalg_COD_decomp (gsl_matrix * A, gsl_vector * tau_Q, gsl_vector * tau_Z, gsl_permutation * p, size_t * rank, gsl_vector * work)
              int gsl_linalg_COD_decomp_e (gsl_matrix * A, gsl_vector * tau_Q, gsl_vector * tau_Z, gsl_permutation * p, double tol, size_t * rank, gsl_vector * work)

   These functions factor the :math:`M`-by-:math:`N` matrix :data:`A` into the decomposition :math:`A = Q R Z P^T`. The rank of :data:`A`
   is computed as the number of diagonal elements of :math:`R` greater than the tolerance :data:`tol` and output in :data:`rank`.
   If :data:`tol` is not specified, a default value is used (see :func:`gsl_linalg_QRPT_rank`). On output, the permutation
   matrix :math:`P` is stored in :data:`p`. The matrix :math:`R_{11}` is stored in the upper :data:`rank`-by-:data:`rank` block of :data:`A`.
   The matrices :math:`Q` and :math:`Z` are encoded in packed storage in :data:`A` on output. The vectors :data:`tau_Q` and :data:`tau_Z`
   contain the Householder scalars corresponding to the matrices :math:`Q` and :math:`Z` respectively and must be
   of length :math:`k = \min(M,N)`. The vector :data:`work` is additional workspace of length :math:`N`.

.. function:: int gsl_linalg_COD_lssolve (const gsl_matrix * QRZT, const gsl_vector * tau_Q, const gsl_vector * tau_Z, const gsl_permutation * p, const size_t rank, const gsl_vector * b, gsl_vector * x, gsl_vector * residual)

   This function finds the unique minimum norm least squares solution to the overdetermined
   system :math:`A x = b` where the matrix :data:`A` has more rows than
   columns.  The least squares solution minimizes the Euclidean norm of the
   residual, :math:`||b - A x||` as well as the norm of the solution :math:`||x||`.  The routine requires as input 
   the :math:`QRZT` decomposition of :math:`A` into (:data:`QRZT`, :data:`tau_Q`, :data:`tau_Z`, :data:`p`, :data:`rank`)
   given by :func:`gsl_linalg_COD_decomp`.  The solution is returned in :data:`x`.  The
   residual, :math:`b - Ax`, is computed as a by-product and stored in :data:`residual`.

.. function:: int gsl_linalg_COD_lssolve2 (const double lambda, const gsl_matrix * QRZT, const gsl_vector * tau_Q, const gsl_vector * tau_Z, const gsl_permutation * p, const size_t rank, const gsl_vector * b, gsl_vector * x, gsl_vector * residual, gsl_matrix * S, gsl_vector * work)

   This function finds the solution to the regularized least squares problem in Tikhonov
   standard form, :math:`\min_x ||b - Ax||^2 + \lambda^2 ||x||^2`. The routine requires as input 
   the :math:`QRZT` decomposition of :math:`A` into (:data:`QRZT`, :data:`tau_Q`, :data:`tau_Z`, :data:`p`, :data:`rank`)
   given by :func:`gsl_linalg_COD_decomp`. The parameter :math:`\lambda` is supplied in :data:`lambda`.  The solution
   is returned in :data:`x`. The residual, :math:`b - Ax`, is stored in :data:`residual` on output. :data:`S` is additional
   workspace of size :data:`rank`-by-:data:`rank`. :data:`work` is additional workspace of length :data:`rank`.

.. function:: int gsl_linalg_COD_unpack (const gsl_matrix * QRZT, const gsl_vector * tau_Q, const gsl_vector * tau_Z, const size_t rank, gsl_matrix * Q, gsl_matrix * R, gsl_matrix * Z)

   This function unpacks the encoded :math:`QRZT` decomposition
   (:data:`QRZT`, :data:`tau_Q`, :data:`tau_Z`, :data:`rank`) into the matrices
   :data:`Q`, :data:`R`, and :data:`Z`, where :data:`Q` is :math:`M`-by-:math:`M`,
   :data:`R` is :math:`M`-by-:math:`N`, and :data:`Z` is :math:`N`-by-:math:`N`.

.. function:: int gsl_linalg_COD_matZ (const gsl_matrix * QRZT, const gsl_vector * tau_Z, const size_t rank, gsl_matrix * A, gsl_vector * work)

   This function multiplies the input matrix :data:`A` on the right by :data:`Z`,
   :math:`A' = A Z` using the encoded :math:`QRZT` decomposition
   (:data:`QRZT`, :data:`tau_Z`, :data:`rank`). :data:`A` must have :math:`N` columns but may
   have any number of rows. Additional workspace of length :math:`M` is provided
   in :data:`work`.

.. index:: SVD, singular value decomposition

Singular Value Decomposition
============================

A general rectangular :math:`M`-by-:math:`N` matrix :math:`A` has a
singular value decomposition (SVD) into the product of an
:math:`M`-by-:math:`N` orthogonal matrix :math:`U`, an :math:`N`-by-:math:`N`
diagonal matrix of singular values :math:`S` and the transpose of an
:math:`N`-by-:math:`N` orthogonal square matrix :math:`V`,

.. math:: A = U S V^T

The singular values :math:`\sigma_i = S_{ii}`
are all non-negative and are
generally chosen to form a non-increasing sequence

.. only:: not texinfo

   .. math:: \sigma_1 \ge \sigma_2 \ge ... \ge \sigma_N \ge 0

.. only:: texinfo

   .. math:: \sigma_1 >= \sigma_2 >= ... >= \sigma_N >= 0

The singular value decomposition of a matrix has many practical uses.
The condition number of the matrix is given by the ratio of the largest
singular value to the smallest singular value. The presence of a zero
singular value indicates that the matrix is singular. The number of
non-zero singular values indicates the rank of the matrix.  In practice
singular value decomposition of a rank-deficient matrix will not produce
exact zeroes for singular values, due to finite numerical
precision.  Small singular values should be edited by choosing a suitable
tolerance.

For a rank-deficient matrix, the null space of :math:`A` is given by
the columns of :math:`V` corresponding to the zero singular values.
Similarly, the range of :math:`A` is given by columns of :math:`U`
corresponding to the non-zero singular values.

Note that the routines here compute the "thin" version of the SVD
with :math:`U` as :math:`M`-by-:math:`N` orthogonal matrix. This allows
in-place computation and is the most commonly-used form in practice.
Mathematically, the "full" SVD is defined with :math:`U` as an
:math:`M`-by-:math:`M` orthogonal matrix and :math:`S` as an
:math:`M`-by-:math:`N` diagonal matrix (with additional rows of zeros).

.. function:: int gsl_linalg_SV_decomp (gsl_matrix * A, gsl_matrix * V, gsl_vector * S, gsl_vector * work)

   This function factorizes the :math:`M`-by-:math:`N` matrix :data:`A` into
   the singular value decomposition :math:`A = U S V^T` for :math:`M \ge N`.
   On output the matrix :data:`A` is replaced by
   :math:`U`. The diagonal elements of the singular value matrix :math:`S`
   are stored in the vector :data:`S`. The singular values are non-negative
   and form a non-increasing sequence from :math:`S_1` to :math:`S_N`. The
   matrix :data:`V` contains the elements of :math:`V` in untransposed
   form. To form the product :math:`U S V^T` it is necessary to take the
   transpose of :data:`V`.  A workspace of length :data:`N` is required in
   :data:`work`.

   This routine uses the Golub-Reinsch SVD algorithm.  

.. function:: int gsl_linalg_SV_decomp_mod (gsl_matrix * A, gsl_matrix * X, gsl_matrix * V, gsl_vector * S, gsl_vector * work)

   This function computes the SVD using the modified Golub-Reinsch
   algorithm, which is faster for :math:`M \gg N`.
   It requires the vector :data:`work` of length :data:`N` and the
   :math:`N`-by-:math:`N` matrix :data:`X` as additional working space.

.. index:: Jacobi orthogonalization

.. function:: int gsl_linalg_SV_decomp_jacobi (gsl_matrix * A, gsl_matrix * V, gsl_vector * S)

   This function computes the SVD of the :math:`M`-by-:math:`N` matrix :data:`A`
   using one-sided Jacobi orthogonalization for :math:`M \ge N`.
   The Jacobi method can compute singular values to higher
   relative accuracy than Golub-Reinsch algorithms (see references for
   details).

.. function:: int gsl_linalg_SV_solve (const gsl_matrix * U, const gsl_matrix * V, const gsl_vector * S, const gsl_vector * b, gsl_vector * x)

   This function solves the system :math:`A x = b` using the singular value
   decomposition (:data:`U`, :data:`S`, :data:`V`) of :math:`A` which must 
   have been computed previously with :func:`gsl_linalg_SV_decomp`.

   Only non-zero singular values are used in computing the solution. The
   parts of the solution corresponding to singular values of zero are
   ignored.  Other singular values can be edited out by setting them to
   zero before calling this function. 

   In the over-determined case where :data:`A` has more rows than columns the
   system is solved in the least squares sense, returning the solution
   :data:`x` which minimizes :math:`||A x - b||_2`.

.. function:: int gsl_linalg_SV_solve2 (const double tol, const gsl_matrix * U, const gsl_matrix * V, const gsl_vector * S, const gsl_vector * b, gsl_vector * x, gsl_vector * work)

   This function solves the system :math:`A x = b` using the singular value
   decomposition (:data:`U`, :data:`S`, :data:`V`) of :math:`A` which must 
   have been computed previously with :func:`gsl_linalg_SV_decomp`.

   Singular values which satisfy, :math:`s_i \leq tol \times s_{max}` are excluded
   from the solution. Additional workspace of length :math:`N` must be provided in
   :data:`work`.

   In the over-determined case where :data:`A` has more rows than columns the
   system is solved in the least squares sense, returning the solution
   :data:`x` which minimizes :math:`||b - A x||_2`.

.. function:: int gsl_linalg_SV_lssolve (const double lambda, const gsl_matrix * U, const gsl_matrix * V, const gsl_vector * S, const gsl_vector * b, gsl_vector * x, double * rnorm, gsl_vector * work)

   This function solves the regularized least squares problem,

   .. math:: \min_x || b - A x ||^2 + \lambda^2 ||x||^2

   using the singular value decomposition (:data:`U`, :data:`S`, :data:`V`) of
   :math:`A` which must have been computed previously with :func:`gsl_linalg_SV_decomp`.
   The residual norm :math:`||b - Ax||` is stored in :data:`rnorm` on output.
   Additional workspace of size :math:`M+N` is required in :data:`work`.

.. function:: int gsl_linalg_SV_leverage (const gsl_matrix * U, gsl_vector * h)

   This function computes the statistical leverage values :math:`h_i` of a matrix :math:`A`
   using its singular value decomposition (:data:`U`, :data:`S`, :data:`V`) previously computed
   with :func:`gsl_linalg_SV_decomp`. :math:`h_i` are the diagonal values of the matrix
   :math:`A (A^T A)^{-1} A^T` and depend only on the matrix :data:`U` which is the input to
   this function.

.. index::
   single: Cholesky decomposition
   single: square root of a matrix, Cholesky decomposition
   single: matrix square root, Cholesky decomposition

.. _sec_cholesky-decomposition:

Cholesky Decomposition
======================

A symmetric, positive definite square matrix :math:`A` has a Cholesky
decomposition into a product of a lower triangular matrix :math:`L` and
its transpose :math:`L^T`,

.. math:: A = L L^T

This is sometimes referred to as taking the square-root of a matrix. The
Cholesky decomposition can only be carried out when all the eigenvalues
of the matrix are positive.  This decomposition can be used to convert
the linear system :math:`A x = b` into a pair of triangular systems
(:math:`L y = b`, :math:`L^T x = y`), which can be solved by forward and
back-substitution.

If the matrix :math:`A` is near singular, it is sometimes possible to reduce
the condition number and recover a more accurate solution vector :math:`x`
by scaling as

.. only:: not texinfo

   .. math:: \left( S A S \right) \left( S^{-1} x \right) = S b

.. only:: texinfo

   .. math:: ( S A S ) ( S^(-1) x ) = S b

where :math:`S` is a diagonal matrix whose elements are given by
:math:`S_{ii} = 1/\sqrt{A_{ii}}`. This scaling is also known as
Jacobi preconditioning. There are routines below to solve
both the scaled and unscaled systems.

.. function:: int gsl_linalg_cholesky_decomp1 (gsl_matrix * A)
              int gsl_linalg_complex_cholesky_decomp (gsl_matrix_complex * A)

   These functions factorize the symmetric, positive-definite square matrix
   :data:`A` into the Cholesky decomposition :math:`A = L L^T` (or
   :math:`A = L L^{\dagger}`
   for the complex case). On input, the values from the diagonal and lower-triangular
   part of the matrix :data:`A` are used (the upper triangular part is ignored).  On output
   the diagonal and lower triangular part of the input matrix :data:`A` contain the matrix
   :math:`L`, while the upper triangular part contains the original matrix.  If the matrix is not
   positive-definite then the decomposition will fail, returning the
   error code :macro:`GSL_EDOM`.

   When testing whether a matrix is positive-definite, disable the error
   handler first to avoid triggering an error. These functions use
   Level 3 BLAS to compute the Cholesky factorization (Peise and Bientinesi, 2016).

.. function:: int gsl_linalg_cholesky_decomp (gsl_matrix * A)

   This function is now deprecated and is provided only for backward compatibility.

.. function:: int gsl_linalg_cholesky_solve (const gsl_matrix * cholesky, const gsl_vector * b, gsl_vector * x)
              int gsl_linalg_complex_cholesky_solve (const gsl_matrix_complex * cholesky, const gsl_vector_complex * b, gsl_vector_complex * x)

   These functions solve the system :math:`A x = b` using the Cholesky
   decomposition of :math:`A` held in the matrix :data:`cholesky` which must
   have been previously computed by :func:`gsl_linalg_cholesky_decomp` or
   :func:`gsl_linalg_complex_cholesky_decomp`.

.. function:: int gsl_linalg_cholesky_svx (const gsl_matrix * cholesky, gsl_vector * x)
              int gsl_linalg_complex_cholesky_svx (const gsl_matrix_complex * cholesky, gsl_vector_complex * x)

   These functions solve the system :math:`A x = b` in-place using the
   Cholesky decomposition of :math:`A` held in the matrix :data:`cholesky`
   which must have been previously computed by
   :func:`gsl_linalg_cholesky_decomp` or
   :func:`gsl_linalg_complex_cholesky_decomp`. On input :data:`x` should
   contain the right-hand side :math:`b`, which is replaced by the
   solution on output.

.. function:: int gsl_linalg_cholesky_invert (gsl_matrix * cholesky)
              int gsl_linalg_complex_cholesky_invert (gsl_matrix_complex * cholesky)

   These functions compute the inverse of a matrix from its Cholesky
   decomposition :data:`cholesky`, which must have been previously computed
   by :func:`gsl_linalg_cholesky_decomp` or
   :func:`gsl_linalg_complex_cholesky_decomp`.  On output, the inverse is
   stored in-place in :data:`cholesky`.

.. function:: int gsl_linalg_cholesky_decomp2 (gsl_matrix * A, gsl_vector * S)
              int gsl_linalg_complex_cholesky_decomp2 (gsl_matrix_complex * A, gsl_vector * S)

   This function calculates a diagonal scaling transformation :math:`S` for
   the symmetric, positive-definite square matrix :data:`A`, and then
   computes the Cholesky decomposition :math:`S A S = L L^T`.
   On input, the values from the diagonal and lower-triangular part of the matrix :data:`A` are
   used (the upper triangular part is ignored).  On output the diagonal and lower triangular part
   of the input matrix :data:`A` contain the matrix :math:`L`.
   If the matrix is not positive-definite then the decomposition will fail, returning the
   error code :macro:`GSL_EDOM`. The diagonal scale factors are stored in :data:`S`
   on output.

   When testing whether a matrix is positive-definite, disable the error
   handler first to avoid triggering an error.

.. function:: int gsl_linalg_cholesky_solve2 (const gsl_matrix * cholesky, const gsl_vector * S, const gsl_vector * b, gsl_vector * x)
              int gsl_linalg_complex_cholesky_solve2 (const gsl_matrix_complex * cholesky, const gsl_vector * S, const gsl_vector_complex * b, gsl_vector_complex * x)

   This function solves the system :math:`(S A S) (S^{-1} x) = S b` using the Cholesky
   decomposition of :math:`S A S` held in the matrix :data:`cholesky` which must
   have been previously computed by :func:`gsl_linalg_cholesky_decomp2` or
   :func:`gsl_linalg_complex_cholesky_decomp2`.

.. function:: int gsl_linalg_cholesky_svx2 (const gsl_matrix * cholesky, const gsl_vector * S, gsl_vector * x)
              int gsl_linalg_complex_cholesky_svx2 (const gsl_matrix_complex * cholesky, const gsl_vector * S, gsl_vector_complex * x)

   This function solves the system :math:`(S A S) (S^{-1} x) = S b` in-place using the
   Cholesky decomposition of :math:`S A S` held in the matrix :data:`cholesky`
   which must have been previously computed by
   :func:`gsl_linalg_cholesky_decomp2` or :func:`gsl_linalg_complex_cholesky_decomp2`.
   On input :data:`x` should contain the right-hand side :math:`b`, which is replaced by the
   solution on output.

.. function:: int gsl_linalg_cholesky_scale (const gsl_matrix * A, gsl_vector * S)
              int gsl_linalg_complex_cholesky_scale (const gsl_matrix_complex * A, gsl_vector * S)

   This function calculates a diagonal scaling transformation of the
   symmetric, positive definite matrix :data:`A`, such that
   :math:`S A S` has a condition number within a factor of :math:`N`
   of the matrix of smallest possible condition number over all
   possible diagonal scalings. On output, :data:`S` contains the
   scale factors, given by :math:`S_i = 1/\sqrt{A_{ii}}`.
   For any :math:`A_{ii} \le 0`, the corresponding scale factor :math:`S_i`
   is set to :math:`1`.

.. function:: int gsl_linalg_cholesky_scale_apply (gsl_matrix * A, const gsl_vector * S)
              int gsl_linalg_complex_cholesky_scale_apply (gsl_matrix_complex * A, const gsl_vector * S)

   This function applies the scaling transformation :data:`S` to the matrix :data:`A`. On output,
   :data:`A` is replaced by :math:`S A S`.

.. function:: int gsl_linalg_cholesky_rcond (const gsl_matrix * cholesky, double * rcond, gsl_vector * work)

   This function estimates the reciprocal condition number (using the 1-norm) of the symmetric positive
   definite matrix :math:`A`, using its Cholesky decomposition provided in :data:`cholesky`.
   The reciprocal condition number estimate, defined as :math:`1 / (||A||_1 \cdot ||A^{-1}||_1)`, is stored
   in :data:`rcond`.  Additional workspace of size :math:`3 N` is required in :data:`work`.

.. index::
   single: Cholesky decomposition, pivoted
   single: Pivoted Cholesky Decomposition

Pivoted Cholesky Decomposition
==============================

A symmetric positive semi-definite square matrix :math:`A` has an alternate
Cholesky decomposition into a product of a lower unit triangular matrix :math:`L`,
a diagonal matrix :math:`D` and :math:`L^T`, given by :math:`L D L^T`. For
postive definite matrices, this is equivalent to the Cholesky formulation discussed
above, with the standard Cholesky lower triangular factor given by :math:`L D^{1 \over 2}`.
For ill-conditioned matrices, it can help to use a pivoting strategy to
prevent the entries of :math:`D` and :math:`L` from growing too large, and also
ensure :math:`D_1 \ge D_2 \ge \cdots \ge D_n > 0`, where :math:`D_i` are
the diagonal entries of :math:`D`. The final decomposition is given by

.. math:: P A P^T = L D L^T

where :math:`P` is a permutation matrix.

.. function:: int gsl_linalg_pcholesky_decomp (gsl_matrix * A, gsl_permutation * p)

   This function factors the symmetric, positive-definite square matrix
   :data:`A` into the Pivoted Cholesky decomposition :math:`P A P^T = L D L^T`.
   On input, the values from the diagonal and lower-triangular part of the matrix :data:`A` are
   used to construct the factorization. On output the diagonal of the input matrix :data:`A` stores
   the diagonal elements of :math:`D`, and the lower triangular portion of :data:`A`
   contains the matrix :math:`L`. Since :math:`L` has ones on its diagonal these
   do not need to be explicitely stored. The upper triangular portion of :data:`A` is
   unmodified. The permutation matrix :math:`P` is stored in :data:`p` on output.

.. function:: int gsl_linalg_pcholesky_solve (const gsl_matrix * LDLT, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)

   This function solves the system :math:`A x = b` using the Pivoted Cholesky
   decomposition of :math:`A` held in the matrix :data:`LDLT` and permutation
   :data:`p` which must have been previously computed by :func:`gsl_linalg_pcholesky_decomp`.

.. function:: int gsl_linalg_pcholesky_svx (const gsl_matrix * LDLT, const gsl_permutation * p, gsl_vector * x)

   This function solves the system :math:`A x = b` in-place using the Pivoted Cholesky
   decomposition of :math:`A` held in the matrix :data:`LDLT` and permutation
   :data:`p` which must have been previously computed by :func:`gsl_linalg_pcholesky_decomp`.
   On input, :data:`x` contains the right hand side vector :math:`b` which is
   replaced by the solution vector on output.

.. function:: int gsl_linalg_pcholesky_decomp2 (gsl_matrix * A, gsl_permutation * p, gsl_vector * S)

   This function computes the pivoted Cholesky factorization of the matrix
   :math:`S A S`, where the input matrix :data:`A` is symmetric and positive
   definite, and the diagonal scaling matrix :data:`S` is computed to reduce the
   condition number of :data:`A` as much as possible. See
   :ref:`Cholesky Decomposition <sec_cholesky-decomposition>` for more information on the matrix :data:`S`.
   The Pivoted Cholesky decomposition satisfies :math:`P S A S P^T = L D L^T`.
   On input, the values from the diagonal and lower-triangular part of the matrix :data:`A` are
   used to construct the factorization.  On output the diagonal of the input matrix :data:`A` stores the diagonal
   elements of :math:`D`, and the lower triangular portion of :data:`A`
   contains the matrix :math:`L`. Since :math:`L` has ones on its diagonal these
   do not need to be explicitely stored. The upper triangular portion of :data:`A`
   is unmodified. The permutation matrix :math:`P` is stored in :data:`p` on output.
   The diagonal scaling transformation is stored in :data:`S` on output.

.. function:: int gsl_linalg_pcholesky_solve2 (const gsl_matrix * LDLT, const gsl_permutation * p, const gsl_vector * S, const gsl_vector * b, gsl_vector * x)

   This function solves the system :math:`(S A S) (S^{-1} x) = S b` using the Pivoted Cholesky
   decomposition of :math:`S A S` held in the matrix :data:`LDLT`, permutation
   :data:`p`, and vector :data:`S`, which must have been previously computed by
   :func:`gsl_linalg_pcholesky_decomp2`.

.. function:: int gsl_linalg_pcholesky_svx2 (const gsl_matrix * LDLT, const gsl_permutation * p, const gsl_vector * S, gsl_vector * x)

   This function solves the system :math:`(S A S) (S^{-1} x) = S b` in-place using the Pivoted Cholesky
   decomposition of :math:`S A S` held in the matrix :data:`LDLT`, permutation
   :data:`p` and vector :data:`S`, which must have been previously computed by
   :func:`gsl_linalg_pcholesky_decomp2`.
   On input, :data:`x` contains the right hand side vector :math:`b` which is
   replaced by the solution vector on output.

.. function:: int gsl_linalg_pcholesky_invert (const gsl_matrix * LDLT, const gsl_permutation * p, gsl_matrix * Ainv)

   This function computes the inverse of the matrix :math:`A`, using the Pivoted
   Cholesky decomposition stored in :data:`LDLT` and :data:`p`. On output, the
   matrix :data:`Ainv` contains :math:`A^{-1}`.

.. function:: int gsl_linalg_pcholesky_rcond (const gsl_matrix * LDLT, const gsl_permutation * p, double * rcond, gsl_vector * work)

   This function estimates the reciprocal condition number (using the 1-norm) of the symmetric positive
   definite matrix :math:`A`, using its pivoted Cholesky decomposition provided in :data:`LDLT`.
   The reciprocal condition number estimate, defined as :math:`1 / (||A||_1 \cdot ||A^{-1}||_1)`, is stored
   in :data:`rcond`.  Additional workspace of size :math:`3 N` is required in :data:`work`.

.. index::
   single: Cholesky decomposition, modified
   single: Modified Cholesky Decomposition

Modified Cholesky Decomposition
===============================

The modified Cholesky decomposition is suitable for solving systems
:math:`A x = b` where :math:`A` is a symmetric indefinite matrix. Such
matrices arise in nonlinear optimization algorithms. The standard
Cholesky decomposition requires a positive definite matrix and would
fail in this case. Instead of resorting to a method like QR or SVD,
which do not take into account the symmetry of the matrix, we can
instead introduce a small perturbation to the matrix :math:`A` to
make it positive definite, and then use a Cholesky decomposition on
the perturbed matrix. The resulting decomposition satisfies

.. math:: P (A + E) P^T = L D L^T

where :math:`P` is a permutation matrix, :math:`E` is a diagonal
perturbation matrix, :math:`L` is unit lower triangular, and
:math:`D` is diagonal. If :math:`A` is sufficiently positive
definite, then the perturbation matrix :math:`E` will be zero
and this method is equivalent to the pivoted Cholesky algorithm.
For indefinite matrices, the perturbation matrix :math:`E` is
computed to ensure that :math:`A + E` is positive definite and
well conditioned.

.. function:: int gsl_linalg_mcholesky_decomp (gsl_matrix * A, gsl_permutation * p, gsl_vector * E)

   This function factors the symmetric, indefinite square matrix
   :data:`A` into the Modified Cholesky decomposition :math:`P (A + E) P^T = L D L^T`.
   On input, the values from the diagonal and lower-triangular part of the matrix :data:`A` are
   used to construct the factorization. On output the diagonal of the input matrix :data:`A` stores the diagonal
   elements of :math:`D`, and the lower triangular portion of :data:`A`
   contains the matrix :math:`L`. Since :math:`L` has ones on its diagonal these
   do not need to be explicitely stored. The upper triangular portion of :data:`A`
   is unmodified. The permutation matrix :math:`P` is
   stored in :data:`p` on output. The diagonal perturbation matrix is stored in
   :data:`E` on output. The parameter :data:`E` may be set to NULL if it is not
   required.

.. function:: int gsl_linalg_mcholesky_solve (const gsl_matrix * LDLT, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)

   This function solves the perturbed system :math:`(A + E) x = b` using the Cholesky
   decomposition of :math:`A + E` held in the matrix :data:`LDLT` and permutation
   :data:`p` which must have been previously computed by :func:`gsl_linalg_mcholesky_decomp`.

.. function:: int gsl_linalg_mcholesky_svx (const gsl_matrix * LDLT, const gsl_permutation * p, gsl_vector * x)

   This function solves the perturbed system :math:`(A + E) x = b` in-place using the Cholesky
   decomposition of :math:`A + E` held in the matrix :data:`LDLT` and permutation
   :data:`p` which must have been previously computed by :func:`gsl_linalg_mcholesky_decomp`.
   On input, :data:`x` contains the right hand side vector :math:`b` which is
   replaced by the solution vector on output.

.. function:: int gsl_linalg_mcholesky_rcond (const gsl_matrix * LDLT, const gsl_permutation * p, double * rcond, gsl_vector * work)

   This function estimates the reciprocal condition number (using the 1-norm) of the perturbed matrix
   :math:`A + E`, using its pivoted Cholesky decomposition provided in :data:`LDLT`.
   The reciprocal condition number estimate, defined as :math:`1 / (||A + E||_1 \cdot ||(A + E)^{-1}||_1)`, is stored
   in :data:`rcond`.  Additional workspace of size :math:`3 N` is required in :data:`work`.

.. index::
   single: LDL decomposition
   single: LDLT decomposition
   single: Cholesky decomposition, square root free

.. _sec_ldlt-decomposition:

LDLT Decomposition
==================

If :math:`A` is a symmetric, nonsingular square matrix, then it has a unique
factorization of the form

.. math:: A = L D L^T

where :math:`L` is a unit lower triangular matrix and :math:`D` is diagonal.
If :math:`A` is positive definite, then this factorization is equivalent
to the Cholesky factorization, where the lower triangular Cholesky factor
is :math:`L D^{\frac{1}{2}}`. Some indefinite matrices for which no
Cholesky decomposition exists have an :math:`L D L^T` decomposition
with negative entries in :math:`D`. The :math:`L D L^T` algorithm
is sometimes referred to as the *square root free* Cholesky
decomposition, as the algorithm does not require the computation of
square roots. The algorithm is stable for positive definite matrices,
but is not guaranteed to be stable for indefinite matrices.

.. function:: int gsl_linalg_ldlt_decomp (gsl_matrix * A)

   This function factorizes the symmetric, non-singular square matrix
   :data:`A` into the decomposition :math:`A = L D L^T`.
   On input, the values from the diagonal and lower-triangular
   part of the matrix :data:`A` are used. The upper triangle of :data:`A`
   is used as temporary workspace.  On output
   the diagonal of :data:`A` contains the matrix :math:`D` and the lower
   triangle of :data:`A` contains the unit lower triangular matrix :math:`L`.
   The matrix 1-norm, :math:`||A||_1` is stored in the upper right corner
   on output, for later use by :func:`gsl_linalg_ldlt_rcond`.

   If the matrix is detected to be singular, the function returns
   the error code :macro:`GSL_EDOM`.

.. function:: int gsl_linalg_ldlt_solve (const gsl_matrix * LDLT, const gsl_vector * b, gsl_vector * x)

   This function solves the system :math:`A x = b` using the :math:`L D L^T`
   decomposition of :math:`A` held in the matrix :data:`LDLT` which must
   have been previously computed by :func:`gsl_linalg_ldlt_decomp`.

.. function:: int gsl_linalg_ldlt_svx (const gsl_matrix * LDLT, gsl_vector * x)

   This function solves the system :math:`A x = b` in-place using the
   :math:`L D L^T` decomposition of :math:`A` held in the matrix :data:`LDLT`
   which must have been previously computed by
   :func:`gsl_linalg_ldlt_decomp`.  On input :data:`x` should
   contain the right-hand side :math:`b`, which is replaced by the
   solution on output.

.. function:: int gsl_linalg_ldlt_rcond (const gsl_matrix * LDLT, double * rcond, gsl_vector * work)

   This function estimates the reciprocal condition number (using the 1-norm) of the symmetric
   nonsingular matrix :math:`A`, using its :math:`L D L^T` decomposition provided in :data:`LDLT`.
   The reciprocal condition number estimate, defined as :math:`1 / (||A||_1 \cdot ||A^{-1}||_1)`, is stored
   in :data:`rcond`.  Additional workspace of size :math:`3 N` is required in :data:`work`.

.. index:: tridiagonal decomposition

Tridiagonal Decomposition of Real Symmetric Matrices
====================================================

A symmetric matrix :math:`A` can be factorized by similarity
transformations into the form,

.. math:: A = Q T Q^T

where :math:`Q` is an orthogonal matrix and :math:`T` is a symmetric
tridiagonal matrix.

.. function:: int gsl_linalg_symmtd_decomp (gsl_matrix * A, gsl_vector * tau)

   This function factorizes the symmetric square matrix :data:`A` into the
   symmetric tridiagonal decomposition :math:`Q T Q^T`.  On output the
   diagonal and subdiagonal part of the input matrix :data:`A` contain the
   tridiagonal matrix :math:`T`.  The remaining lower triangular part of the
   input matrix contains the Householder vectors which, together with the
   Householder coefficients :data:`tau`, encode the orthogonal matrix
   :math:`Q`. This storage scheme is the same as used by |lapack|.  The
   upper triangular part of :data:`A` is not referenced.

.. function:: int gsl_linalg_symmtd_unpack (const gsl_matrix * A, const gsl_vector * tau, gsl_matrix * Q, gsl_vector * diag, gsl_vector * subdiag)

   This function unpacks the encoded symmetric tridiagonal decomposition
   (:data:`A`, :data:`tau`) obtained from :func:`gsl_linalg_symmtd_decomp` into
   the orthogonal matrix :data:`Q`, the vector of diagonal elements :data:`diag`
   and the vector of subdiagonal elements :data:`subdiag`.  

.. function:: int gsl_linalg_symmtd_unpack_T (const gsl_matrix * A, gsl_vector * diag, gsl_vector * subdiag)

   This function unpacks the diagonal and subdiagonal of the encoded
   symmetric tridiagonal decomposition (:data:`A`, :data:`tau`) obtained from
   :func:`gsl_linalg_symmtd_decomp` into the vectors :data:`diag` and :data:`subdiag`.

.. index:: tridiagonal decomposition

Tridiagonal Decomposition of Hermitian Matrices
===============================================

A hermitian matrix :math:`A` can be factorized by similarity
transformations into the form,

.. math:: A = U T U^T

where :math:`U` is a unitary matrix and :math:`T` is a real symmetric
tridiagonal matrix.

.. function:: int gsl_linalg_hermtd_decomp (gsl_matrix_complex * A, gsl_vector_complex * tau)

   This function factorizes the hermitian matrix :data:`A` into the symmetric
   tridiagonal decomposition :math:`U T U^T`.  On output the real parts of
   the diagonal and subdiagonal part of the input matrix :data:`A` contain
   the tridiagonal matrix :math:`T`.  The remaining lower triangular part of
   the input matrix contains the Householder vectors which, together with
   the Householder coefficients :data:`tau`, encode the unitary matrix
   :math:`U`. This storage scheme is the same as used by |lapack|.  The
   upper triangular part of :data:`A` and imaginary parts of the diagonal are
   not referenced.

.. function:: int gsl_linalg_hermtd_unpack (const gsl_matrix_complex * A, const gsl_vector_complex * tau, gsl_matrix_complex * U, gsl_vector * diag, gsl_vector * subdiag)

   This function unpacks the encoded tridiagonal decomposition (:data:`A`,
   :data:`tau`) obtained from :func:`gsl_linalg_hermtd_decomp` into the
   unitary matrix :data:`U`, the real vector of diagonal elements :data:`diag` and
   the real vector of subdiagonal elements :data:`subdiag`. 

.. function:: int gsl_linalg_hermtd_unpack_T (const gsl_matrix_complex * A, gsl_vector * diag, gsl_vector * subdiag)

   This function unpacks the diagonal and subdiagonal of the encoded
   tridiagonal decomposition (:data:`A`, :data:`tau`) obtained from the
   :func:`gsl_linalg_hermtd_decomp` into the real vectors
   :data:`diag` and :data:`subdiag`.

.. index:: Hessenberg decomposition

Hessenberg Decomposition of Real Matrices
=========================================

A general real matrix :math:`A` can be decomposed by orthogonal
similarity transformations into the form

.. math:: A = U H U^T

where :math:`U` is orthogonal and :math:`H` is an upper Hessenberg matrix,
meaning that it has zeros below the first subdiagonal. The
Hessenberg reduction is the first step in the Schur decomposition
for the nonsymmetric eigenvalue problem, but has applications in
other areas as well.

.. function:: int gsl_linalg_hessenberg_decomp (gsl_matrix * A, gsl_vector * tau)

   This function computes the Hessenberg decomposition of the matrix
   :data:`A` by applying the similarity transformation :math:`H = U^T A U`.
   On output, :math:`H` is stored in the upper portion of :data:`A`. The
   information required to construct the matrix :math:`U` is stored in
   the lower triangular portion of :data:`A`. :math:`U` is a product
   of :math:`N - 2` Householder matrices. The Householder vectors
   are stored in the lower portion of :data:`A` (below the subdiagonal)
   and the Householder coefficients are stored in the vector :data:`tau`.
   :data:`tau` must be of length :data:`N`.

.. function:: int gsl_linalg_hessenberg_unpack (gsl_matrix * H, gsl_vector * tau, gsl_matrix * U)

   This function constructs the orthogonal matrix :math:`U` from the
   information stored in the Hessenberg matrix :data:`H` along with the
   vector :data:`tau`. :data:`H` and :data:`tau` are outputs from
   :func:`gsl_linalg_hessenberg_decomp`.

.. function:: int gsl_linalg_hessenberg_unpack_accum (gsl_matrix * H, gsl_vector * tau, gsl_matrix * V)

   This function is similar to :func:`gsl_linalg_hessenberg_unpack`, except
   it accumulates the matrix :data:`U` into :data:`V`, so that :math:`V' = VU`.
   The matrix :data:`V` must be initialized prior to calling this function.
   Setting :data:`V` to the identity matrix provides the same result as
   :func:`gsl_linalg_hessenberg_unpack`. If :data:`H` is order :data:`N`, then
   :data:`V` must have :data:`N` columns but may have any number of rows.

.. function:: int gsl_linalg_hessenberg_set_zero (gsl_matrix * H)

   This function sets the lower triangular portion of :data:`H`, below
   the subdiagonal, to zero. It is useful for clearing out the
   Householder vectors after calling :func:`gsl_linalg_hessenberg_decomp`.

.. index:: Hessenberg triangular decomposition

Hessenberg-Triangular Decomposition of Real Matrices
====================================================

A general real matrix pair (:math:`A`, :math:`B`) can be decomposed by
orthogonal similarity transformations into the form

.. only:: not texinfo

   .. math::

      A &= U H V^T \\
      B &= U R V^T

.. only:: texinfo

   ::

      A = U H V^T
      B = U R V^T

where :math:`U` and :math:`V` are orthogonal, :math:`H` is an upper
Hessenberg matrix, and :math:`R` is upper triangular. The
Hessenberg-Triangular reduction is the first step in the generalized
Schur decomposition for the generalized eigenvalue problem.

.. function:: int gsl_linalg_hesstri_decomp (gsl_matrix * A, gsl_matrix * B, gsl_matrix * U, gsl_matrix * V, gsl_vector * work)

   This function computes the Hessenberg-Triangular decomposition of the
   matrix pair (:data:`A`, :data:`B`). On output, :math:`H` is stored in :data:`A`,
   and :math:`R` is stored in :data:`B`. If :data:`U` and :data:`V` are provided
   (they may be null), the similarity transformations are stored in them.
   Additional workspace of length :math:`N` is needed in :data:`work`.

.. index:: bidiagonalization of real matrices

Bidiagonalization
=================

A general matrix :math:`A` can be factorized by similarity
transformations into the form,

.. math:: A = U B V^T

where :math:`U` and :math:`V` are orthogonal matrices and :math:`B` is a
:math:`N`-by-:math:`N` bidiagonal matrix with non-zero entries only on the
diagonal and superdiagonal.  The size of :data:`U` is :math:`M`-by-:math:`N`
and the size of :data:`V` is :math:`N`-by-:math:`N`.

.. function:: int gsl_linalg_bidiag_decomp (gsl_matrix * A, gsl_vector * tau_U, gsl_vector * tau_V)

   This function factorizes the :math:`M`-by-:math:`N` matrix :data:`A` into
   bidiagonal form :math:`U B V^T`.  The diagonal and superdiagonal of the
   matrix :math:`B` are stored in the diagonal and superdiagonal of :data:`A`.
   The orthogonal matrices :math:`U` and :data:`V` are stored as compressed
   Householder vectors in the remaining elements of :data:`A`.  The
   Householder coefficients are stored in the vectors :data:`tau_U` and
   :data:`tau_V`.  The length of :data:`tau_U` must equal the number of
   elements in the diagonal of :data:`A` and the length of :data:`tau_V` should
   be one element shorter.

.. function:: int gsl_linalg_bidiag_unpack (const gsl_matrix * A, const gsl_vector * tau_U, gsl_matrix * U, const gsl_vector * tau_V, gsl_matrix * V, gsl_vector * diag, gsl_vector * superdiag)

   This function unpacks the bidiagonal decomposition of :data:`A` produced by
   :func:`gsl_linalg_bidiag_decomp`, (:data:`A`, :data:`tau_U`, :data:`tau_V`)
   into the separate orthogonal matrices :data:`U`, :data:`V` and the diagonal
   vector :data:`diag` and superdiagonal :data:`superdiag`.  Note that :data:`U`
   is stored as a compact :math:`M`-by-:math:`N` orthogonal matrix satisfying
   :math:`U^T U = I` for efficiency.

.. function:: int gsl_linalg_bidiag_unpack2 (gsl_matrix * A, gsl_vector * tau_U, gsl_vector * tau_V, gsl_matrix * V)

   This function unpacks the bidiagonal decomposition of :data:`A` produced by
   :func:`gsl_linalg_bidiag_decomp`, (:data:`A`, :data:`tau_U`, :data:`tau_V`)
   into the separate orthogonal matrices :data:`U`, :data:`V` and the diagonal
   vector :data:`diag` and superdiagonal :data:`superdiag`.  The matrix :data:`U`
   is stored in-place in :data:`A`.

.. function:: int gsl_linalg_bidiag_unpack_B (const gsl_matrix * A, gsl_vector * diag, gsl_vector * superdiag)

   This function unpacks the diagonal and superdiagonal of the bidiagonal
   decomposition of :data:`A` from :func:`gsl_linalg_bidiag_decomp`, into
   the diagonal vector :data:`diag` and superdiagonal vector :data:`superdiag`.

.. index:: Givens rotation

Givens Rotations
================

A Givens rotation is a rotation in the plane acting on two elements
of a given vector. It can be represented in matrix form as

.. only:: not texinfo

   .. math::

      G(i,j,\theta) =
      \left(
      \begin{matrix}
        1 & \ldots & 0 & \ldots & 0 & \ldots & 0 \\
        \vdots & \ddots & \vdots &  & \vdots & & \vdots \\
        0 & \ldots & \cos{\theta} & \ldots & -\sin{\theta} & \ldots & 0 \\
        \vdots &  & \vdots & \ddots & \vdots & & \vdots \\
        0 & \ldots & \sin{\theta} & \ldots & \cos{\theta} & \ldots & 0 \\
        \vdots &  & \vdots &  & \vdots & \ddots & \vdots \\
        0 & \ldots & 0 & \ldots & 0 & \ldots & 1
      \end{matrix}
      \right)

where the :math:`\cos{\theta}` and :math:`\sin{\theta}` appear at
the intersection of the :math:`i`-th and :math:`j`-th rows and columns.
When acting on a vector :math:`x`, :math:`G(i,j,\theta) x` performs
a rotation of the :math:`(i,j)` elements of :math:`x`. Givens
rotations are typically used to introduce zeros in
vectors, such as during the QR decomposition of a matrix. In this
case, it is typically desired to find :math:`c` and :math:`s` such that

.. only:: not texinfo

   .. math::

      \left(
      \begin{matrix}
        c & -s \\
        s & c
      \end{matrix}
      \right)
      \left(
      \begin{matrix}
        a \\
        b
      \end{matrix}
      \right) =
      \left(
      \begin{matrix}
        r \\
        0
      \end{matrix}
      \right)

.. only:: texinfo

   ::

     [ c -s ] [ a ] = [ r ]
     [ s  c ] [ b ]   [ 0 ]

with :math:`r = \sqrt{a^2 + b^2}`.

.. function:: void gsl_linalg_givens (const double a, const double b, double * c, double * s)

   This function computes :math:`c = \cos{\theta}` and :math:`s = \sin{\theta}`
   so that the Givens matrix :math:`G(\theta)` acting on the
   vector :math:`(a,b)` produces :math:`(r, 0)`, with :math:`r = \sqrt{a^2 + b^2}`.

.. function:: void gsl_linalg_givens_gv (gsl_vector * v, const size_t i, const size_t j, const double c, const double s)

   This function applies the Givens rotation defined by
   :math:`c = \cos{\theta}` and :math:`s = \sin{\theta}` to the :data:`i`
   and :data:`j` elements of :data:`v`. On output,
   :math:`(v(i),v(j)) \leftarrow G(\theta) (v(i),v(j))`.

.. index::
   single: Householder matrix
   single: Householder transformation
   single: transformation, Householder

Householder Transformations
===========================

A Householder transformation is a rank-1 modification of the identity
matrix which can be used to zero out selected elements of a vector.  A
Householder matrix :math:`H` takes the form,

.. math:: H = I - \tau v v^T

where :math:`v` is a vector (called the *Householder vector*) and
:math:`\tau = 2/(v^T v)`. The functions described in this section use the
rank-1 structure of the Householder matrix to create and apply
Householder transformations efficiently.

.. function:: double gsl_linalg_householder_transform (gsl_vector * w)
              gsl_complex gsl_linalg_complex_householder_transform (gsl_vector_complex * w)

   This function prepares a Householder transformation :math:`H = I - \tau v v^T`
   which can be used to zero all the elements of the input vector :data:`w`
   except the first. On output the Householder vector :data:`v` is stored in
   :data:`w` and the scalar :math:`\tau` is returned. The householder vector
   :data:`v` is normalized so that :code:`v[0] = 1`, however this 1 is not
   stored in the output vector. Instead, :code:`w[0]` is set to
   the first element of the transformed vector, so that if
   :math:`u = H w`, :code:`w[0] = u[0]` on output and the remainder
   of :math:`u` is zero.

.. function:: int gsl_linalg_householder_hm (double tau, const gsl_vector * v, gsl_matrix * A)
              int gsl_linalg_complex_householder_hm (gsl_complex tau, const gsl_vector_complex * v, gsl_matrix_complex * A)

   This function applies the Householder matrix :math:`H` defined by the
   scalar :data:`tau` and the vector :data:`v` to the left-hand side of the
   matrix :data:`A`. On output the result :math:`H A` is stored in :data:`A`.

.. function:: int gsl_linalg_householder_mh (double tau, const gsl_vector * v, gsl_matrix * A)
              int gsl_linalg_complex_householder_mh (gsl_complex tau, const gsl_vector_complex * v, gsl_matrix_complex * A)

   This function applies the Householder matrix :math:`H` defined by the
   scalar :data:`tau` and the vector :data:`v` to the right-hand side of the
   matrix :data:`A`. On output the result :math:`A H` is stored in :data:`A`.

.. function:: int gsl_linalg_householder_hv (double tau, const gsl_vector * v, gsl_vector * w)
              int gsl_linalg_complex_householder_hv (gsl_complex tau, const gsl_vector_complex * v, gsl_vector_complex * w)

   This function applies the Householder transformation :math:`H` defined by
   the scalar :data:`tau` and the vector :data:`v` to the vector :data:`w`.  On
   output the result :math:`H w` is stored in :data:`w`.

.. @deftypefun int gsl_linalg_householder_hm1 (double tau, gsl_matrix * A)
.. This function applies the Householder transform, defined by the scalar
.. :data:`tau` and the vector :data:`v`, to a matrix being build up from the
.. identity matrix, using the first column of :data:`A` as a householder vector.
.. @end deftypefun

.. index::
   single: solution of linear system by Householder transformations
   single: Householder linear solver

Householder solver for linear systems
=====================================

.. function:: int gsl_linalg_HH_solve (gsl_matrix * A, const gsl_vector * b, gsl_vector * x)

   This function solves the system :math:`A x = b` directly using
   Householder transformations. On output the solution is stored in :data:`x`
   and :data:`b` is not modified. The matrix :data:`A` is destroyed by the
   Householder transformations.

.. function:: int gsl_linalg_HH_svx (gsl_matrix * A, gsl_vector * x)

   This function solves the system :math:`A x = b` in-place using
   Householder transformations.  On input :data:`x` should contain the
   right-hand side :math:`b`, which is replaced by the solution on output.  The
   matrix :data:`A` is destroyed by the Householder transformations.

.. index:: tridiagonal systems

Tridiagonal Systems
===================

The functions described in this section efficiently solve symmetric,
non-symmetric and cyclic tridiagonal systems with minimal storage.
Note that the current implementations of these functions use a variant
of Cholesky decomposition, so the tridiagonal matrix must be positive
definite.  For non-positive definite matrices, the functions return
the error code :macro:`GSL_ESING`.

.. function:: int gsl_linalg_solve_tridiag (const gsl_vector * diag, const gsl_vector * e, const gsl_vector * f, const gsl_vector * b, gsl_vector * x)

   This function solves the general :math:`N`-by-:math:`N` system :math:`A x = b`
   where :data:`A` is tridiagonal (:math:`N \geq 2`).
   The super-diagonal and
   sub-diagonal vectors :data:`e` and :data:`f` must be one element shorter
   than the diagonal vector :data:`diag`.  The form of :data:`A` for the 4-by-4
   case is shown below,

   .. only:: not texinfo

      .. math::

         A =
         \left(
         \begin{matrix}
           d_0&e_0&  0& 0\\
           f_0&d_1&e_1& 0\\
           0  &f_1&d_2&e_2\\ 
           0  &0  &f_2&d_3
         \end{matrix}
         \right)

   .. only:: texinfo

      ::

         A = ( d_0 e_0  0   0  )
             ( f_0 d_1 e_1  0  )
             (  0  f_1 d_2 e_2 )
             (  0   0  f_2 d_3 )

.. function:: int gsl_linalg_solve_symm_tridiag (const gsl_vector * diag, const gsl_vector * e, const gsl_vector * b, gsl_vector * x)

   This function solves the general :math:`N`-by-:math:`N` system :math:`A x = b`
   where :data:`A` is symmetric tridiagonal (:math:`N \geq 2`).
   The off-diagonal vector
   :data:`e` must be one element shorter than the diagonal vector :data:`diag`.
   The form of :data:`A` for the 4-by-4 case is shown below,

   .. only:: not texinfo

      .. math::

         A =
         \left(
         \begin{matrix}
           d_0&e_0&  0& 0\\
           e_0&d_1&e_1& 0\\
           0  &e_1&d_2&e_2\\ 
           0  &0  &e_2&d_3
         \end{matrix}
         \right)

   .. only:: texinfo

      ::

         A = ( d_0 e_0  0   0  )
             ( e_0 d_1 e_1  0  )
             (  0  e_1 d_2 e_2 )
             (  0   0  e_2 d_3 )

.. function:: int gsl_linalg_solve_cyc_tridiag (const gsl_vector * diag, const gsl_vector * e, const gsl_vector * f, const gsl_vector * b, gsl_vector * x)

   This function solves the general :math:`N`-by-:math:`N` system :math:`A x = b`
   where :data:`A` is cyclic tridiagonal (:math:`N \geq 3`).
   The cyclic super-diagonal and
   sub-diagonal vectors :data:`e` and :data:`f` must have the same number of
   elements as the diagonal vector :data:`diag`.  The form of :data:`A` for the
   4-by-4 case is shown below,

   .. only:: not texinfo

      .. math::

         A =
         \left(
         \begin{matrix}
           d_0&e_0& 0 &f_3\\
           f_0&d_1&e_1& 0 \\
           0 &f_1&d_2&e_2\\ 
           e_3& 0 &f_2&d_3
         \end{matrix}
         \right)

   .. only:: texinfo

      ::

         A = ( d_0 e_0  0  f_3 )
             ( f_0 d_1 e_1  0  )
             (  0  f_1 d_2 e_2 )
             ( e_3  0  f_2 d_3 )

.. function:: int gsl_linalg_solve_symm_cyc_tridiag (const gsl_vector * diag, const gsl_vector * e, const gsl_vector * b, gsl_vector * x)

   This function solves the general :math:`N`-by-:math:`N` system :math:`A x = b`
   where :data:`A` is symmetric cyclic tridiagonal (:math:`N \geq 3`).
   The cyclic
   off-diagonal vector :data:`e` must have the same number of elements as the
   diagonal vector :data:`diag`.  The form of :data:`A` for the 4-by-4 case is
   shown below,

   .. only:: not texinfo

      .. math::

         A =
         \left(
         \begin{matrix}
           d_0&e_0& 0 &e_3\\
           e_0&d_1&e_1& 0 \\
           0 &e_1&d_2&e_2\\ 
           e_3& 0 &e_2&d_3
         \end{matrix}
         \right)

   .. only:: texinfo

      ::

         A = ( d_0 e_0  0  e_3 )
             ( e_0 d_1 e_1  0  )
             (  0  e_1 d_2 e_2 )
             ( e_3  0  e_2 d_3 )

.. index:: triangular systems

Triangular Systems
==================

.. function:: int gsl_linalg_tri_invert (CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix * T)
              int gsl_linalg_complex_tri_invert (CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix_complex * T)

   These functions compute the in-place inverse of the triangular matrix :data:`T`, stored in
   the lower triangle when :data:`Uplo` = :code:`CblasLower` and upper triangle
   when :data:`Uplo` = :code:`CblasUpper`. The parameter :data:`Diag` = :code:`CblasUnit`, :code:`CblasNonUnit`
   specifies whether the matrix is unit triangular.

.. function:: int gsl_linalg_tri_LTL (gsl_matrix * L)
              int gsl_linalg_complex_tri_LHL (gsl_matrix_complex * L)

   These functions compute the product :math:`L^T L` (or :math:`L^{\dagger} L`)
   in-place and stores it in the lower triangle of :data:`L` on output.

.. function:: int gsl_linalg_tri_UL (gsl_matrix * LU)
              int gsl_linalg_complex_tri_UL (gsl_matrix_complex * LU)

   These functions compute the product :math:`U L` where :math:`U` is upper
   triangular and :math:`L` is unit lower triangular, stored in :data:`LU`, as
   computed by :func:`gsl_linalg_LU_decomp` or :func:`gsl_linalg_complex_LU_decomp`.
   The product is computed in-place using Level 3 BLAS.

.. function:: int gsl_linalg_tri_rcond (CBLAS_UPLO_t Uplo, const gsl_matrix * A, double * rcond, gsl_vector * work)

   This function estimates the 1-norm reciprocal condition number of the triangular matrix :data:`A`,
   using the lower triangle when :data:`Uplo` is :code:`CblasLower` and upper triangle when
   :data:`Uplo` is :code:`CblasUpper`. The reciprocal condition number
   :math:`1 / \left( \left|\left| A \right|\right|_1 \left|\left| A^{-1} \right|\right|_1 \right)`
   is stored in :data:`rcond` on output. Additional workspace of size :math:`3N` is required
   in :data:`work`.

.. index::
   single: banded matrices
   single: matrices, banded

Banded Systems
==============

Band matrices are sparse matrices whose non-zero entries are confined to a diagonal
*band*. From a storage point of view, significant savings can be achieved by
storing only the non-zero diagonals of a banded matrix. Algorithms such as LU and
Cholesky factorizations preserve the band structure of these matrices. Computationally,
working with compact banded matrices is preferable to working on the full dense
matrix with many zero entries.

.. index::
   single: banded general matrices

General Banded Format
---------------------

An example of a general banded matrix is given below.

.. math::

   A = \begin{pmatrix}
         \alpha_1 & \beta_1 & \gamma_1 & 0 & 0 & 0 \\
         \delta_1 & \alpha_2 & \beta_2 & \gamma_2 & 0 & 0 \\
         0 & \delta_2 & \alpha_3 & \beta_3 & \gamma_3 & 0 \\
         0 & 0 & \delta_3 & \alpha_4 & \beta_4 & \gamma_4 \\
         0 & 0 & 0 & \delta_4 & \alpha_5 & \beta_5 \\
         0 & 0 & 0 & 0 & \delta_5 & \alpha_6
       \end{pmatrix}

This matrix has a *lower bandwidth* of 1 and an *upper bandwidth* of 2.
The lower bandwidth is the number of non-zero subdiagonals, and the
upper bandwidth is the number of non-zero superdiagonals. A
:math:`(p,q)` banded matrix has a lower bandwidth :math:`p` and
upper bandwidth :math:`q`. For example, diagonal matrices are
:math:`(0,0)`, tridiagonal matrices are :math:`(1,1)`, and
upper triangular matrices are :math:`(0,N-1)` banded matrices.

The corresponding :math:`6`-by-:math:`4` packed banded matrix looks like

.. math::

   AB = \begin{pmatrix}
          * & * & \alpha_1 & \delta_1 \\
          * & \beta_1 & \alpha_2 & \delta_2 \\
          \gamma_1 & \beta_2 & \alpha_3 & \delta_3 \\
          \gamma_2 & \beta_3 & \alpha_4 & \delta_4 \\
          \gamma_3 & \beta_4 & \alpha_5 & \delta_5 \\
          \gamma_4 & \beta_5 & \alpha_6 & *
        \end{pmatrix}

where the superdiagonals are stored in columns, followed by the diagonal,
followed by the subdiagonals.
The entries marked by :math:`*` are not referenced by the banded
routines. With this format, each row of :math:`AB` corresponds to
the non-zero entries of the corresponding column of :math:`A`.
For an :math:`N`-by-:math:`N` matrix :math:`A`, the dimension
of :math:`AB` will be :math:`N`-by-:math:`(p+q+1)`.

.. index::
   single: banded symmetric matrices
   single: symmetric matrices, banded

.. _sec_symmetric-banded:

Symmetric Banded Format
-----------------------

Symmetric banded matrices allow for additional storage savings.
As an example, consider the following :math:`6 \times 6` symmetric banded matrix
with lower bandwidth :math:`p = 2`:

.. math::

   A = \begin{pmatrix}
         \alpha_1 & \beta_1 & \gamma_1 & 0 & 0 & 0 \\
         \beta_1 & \alpha_2 & \beta_2 & \gamma_2 & 0 & 0 \\
         \gamma_1 & \beta_2 & \alpha_3 & \beta_3 & \gamma_3 & 0 \\
         0 & \gamma_2 & \beta_3 & \alpha_4 & \beta_4 & \gamma_4 \\
         0 & 0 & \gamma_3 & \beta_4 & \alpha_5 & \beta_5 \\
         0 & 0 & 0 & \gamma_4 & \beta_5 & \alpha_6
       \end{pmatrix}

The packed symmetric banded :math:`6 \times 3` matrix will look like:

.. math::

   AB = \begin{pmatrix}
          \alpha_1 & \beta_1 & \gamma_1 \\
          \alpha_2 & \beta_2 & \gamma_2 \\
          \alpha_3 & \beta_3 & \gamma_3 \\
          \alpha_4 & \beta_4 & \gamma_4 \\
          \alpha_5 & \beta_5 & * \\
          \alpha_6 & * & *
        \end{pmatrix}

The entries marked by :math:`*` are not referenced by the symmetric banded
routines. The relationship between the packed format and original matrix is,

.. math:: AB(i,j) = A(i, i + j) = A(i + j, i)

for :math:`i = 0, \dots, N - 1, j = 0, \dots, p`. Conversely,

.. math:: A(i,j) = AB(j, i - j)

for :math:`i = 0, \dots, N - 1, j = \textrm{max}(0,i-p), \dots, i`.

.. warning::

   Note that this format is the transpose of the symmetric banded format used by
   |LAPACK|. In order to develop efficient routines for symmetric banded matrices,
   it helps to have the nonzero elements in each column in contiguous memory locations.
   Since C uses row-major order, GSL stores the columns in the rows of the packed
   banded format, while |LAPACK|, written in Fortran, uses the transposed format.

.. index::
   single: LU decomposition, banded
   single: banded LU Decomposition

Banded LU Decomposition
-----------------------

The routines in this section are designed to factor banded
:math:`M`-by-:math:`N` matrices with an LU factorization,
:math:`P A = L U`. The matrix :math:`A` is banded of type
:math:`(p,q)`, i.e. a lower bandwidth of :math:`p` and an upper
bandwidth of :math:`q`. See :ref:`LU Decomposition <sec_lu-decomposition>`
for more information on the factorization. For banded :math:`(p,q)`
matrices, the :math:`U` factor will have an upper bandwidth of
:math:`p + q`, while the :math:`L` factor will have a lower bandwidth of
at most :math:`p`. Therefore, additional storage is needed to store the
:math:`p` additional bands of :math:`U`.

As an example, consider the :math:`M = N = 7` matrix with
lower bandwidth :math:`p = 3` and upper bandwidth :math:`q = 2`,

.. math::

   A = \begin{pmatrix}
         \alpha_1   & \beta_1    & \gamma_1   & 0          & 0          & 0        & 0 \\
         \delta_1   & \alpha_2   & \beta_2    & \gamma_2   & 0          & 0        & 0 \\
         \epsilon_1 & \delta_2   & \alpha_3   & \beta_3    & \gamma_3   & 0        & 0 \\
         \zeta_1    & \epsilon_2 & \delta_3   & \alpha_4   & \beta_4    & \gamma_4 & 0 \\
         0          & \zeta_2    & \epsilon_3 & \delta_4   & \alpha_5   & \beta_5  & \gamma_5 \\
         0          & 0          & \zeta_3    & \epsilon_4 & \delta_5   & \alpha_6 & \beta_6 \\
         0          & 0          & 0          & \zeta_4    & \epsilon_5 & \delta_6 & \alpha_7
       \end{pmatrix}

The corresponding :math:`N`-by-:math:`2p + q + 1` packed banded matrix looks like

.. math::

   AB = \begin{pmatrix}
          * & * & * & * & * & \alpha_1 & \delta_1 & \epsilon_1 & \zeta_1 \\
          * & * & * & * & \beta_1 & \alpha_2 & \delta_2 & \epsilon_2 & \zeta_2 \\
          * & * & * & \gamma_1 & \beta_2 & \alpha_3 & \delta_3 & \epsilon_3 & \zeta_3 \\
          * & * & - & \gamma_2 & \beta_3 & \alpha_4 & \delta_4 & \epsilon_4 & \zeta_4 \\
          * & - & - & \gamma_3 & \beta_4 & \alpha_5 & \delta_5 & \epsilon_5 & * \\
          - & - & - & \gamma_4 & \beta_5 & \alpha_6 & \delta_6 & * & * \\
          \undermat{p}{- & - & -} & \undermat{q}{\gamma_5 & \beta_6} & \alpha_7 & \undermat{p}{* & * & * & }
        \end{pmatrix}

Entries marked with :math:`-` are used to store the additional :math:`p` diagonals of the
:math:`U` factor. Entries marked with :math:`*` are not referenced by the banded routines.

.. function:: int gsl_linalg_LU_band_decomp (const size_t M, const size_t lb, const size_t ub, gsl_matrix * AB, gsl_vector_uint * piv)

   This function computes the LU factorization of the banded matrix :data:`AB` which
   is stored in packed band format (see above) and has dimension :math:`N`-by-:math:`2p + q + 1`. The
   number of rows :math:`M` of the original matrix is provided in :data:`M`.
   The lower bandwidth :math:`p` is provided in :data:`lb` and
   the upper bandwidth :math:`q` is provided in :data:`ub`. The vector :data:`piv` has length :math:`\textrm{min}(M,N)`
   and stores the pivot indices on output (for :math:`0 \le i < \textrm{min}(M,N)`, row :math:`i` of the matrix
   was interchanged with row :code:`piv[i]`). On output, :data:`AB` contains both the :math:`L` and :math:`U` factors
   in packed format.

.. function:: int gsl_linalg_LU_band_solve (const size_t lb, const size_t ub, const gsl_matrix * LUB, const gsl_vector_uint * piv, const gsl_vector * b, gsl_vector * x)

   This function solves the square system :math:`Ax = b` using the banded LU factorization
   (:data:`LUB`, :data:`piv`) computed by :func:`gsl_linalg_LU_band_decomp`. The lower and
   upper bandwidths are provided in :data:`lb` and :data:`ub` respectively. The right
   hand side vector is provided in :data:`b`. The solution vector is stored in
   :data:`x` on output.

.. function:: int gsl_linalg_LU_band_svx (const size_t lb, const size_t ub, const gsl_matrix * LUB, const gsl_vector_uint * piv, gsl_vector * x)

   This function solves the square system :math:`Ax = b` in-place, using the banded LU factorization
   (:data:`LUB`, :data:`piv`) computed by :func:`gsl_linalg_LU_band_decomp`. The lower and
   upper bandwidths are provided in :data:`lb` and :data:`ub` respectively. On input,
   the right hand side vector :math:`b` is provided in :data:`x`, which is replaced by the
   solution vector :math:`x` on output.

.. function:: int gsl_linalg_LU_band_unpack (const size_t M, const size_t lb, const size_t ub, const gsl_matrix * LUB, const gsl_vector_uint * piv, gsl_matrix * L, gsl_matrix * U)

   This function unpacks the banded LU factorization (:data:`LUB`, :data:`piv`) previously
   computed by :func:`gsl_linalg_LU_band_decomp` into the matrices :data:`L` and :data:`U`.
   The matrix :data:`U` has dimension :math:`\textrm{min}(M,N)`-by-:math:`N` and stores
   the upper triangular factor on output. The matrix :data:`L` has dimension
   :math:`M`-by-:math:`\textrm{min}(M,N)` and stores the matrix :math:`P^T L` on output.

.. index::
   single: Cholesky decomposition, banded
   single: banded Cholesky Decomposition

Banded Cholesky Decomposition
-----------------------------

The routines in this section are designed to factor and solve
:math:`N`-by-:math:`N` linear systems of the form :math:`A x = b`
where :math:`A` is a banded, symmetric, and positive definite matrix
with lower bandwidth :math:`p`. See
:ref:`Cholesky Decomposition <sec_cholesky-decomposition>` for more
information on the factorization.  The lower triangular
factor of the Cholesky decomposition preserves the same
banded structure as the matrix :math:`A`, enabling an efficient
algorithm which overwrites the original matrix with the
:math:`L` factor.

.. function:: int gsl_linalg_cholesky_band_decomp(gsl_matrix * A)

   This function factorizes the symmetric, positive-definite square matrix
   :data:`A` into the Cholesky decomposition :math:`A = L L^T`. The input matrix
   :data:`A` is given in :ref:`symmetric banded format <sec_symmetric-banded>`, and has dimensions
   :math:`N`-by-:math:`(p + 1)`, where :math:`p` is the lower bandwidth of the matrix.
   On output, the entries of :data:`A` are replaced by the entries of the matrix
   :math:`L` in the same format. In addition, the lower right element of :data:`A` is
   used to store the matrix 1-norm, used later by :func:`gsl_linalg_cholesky_band_rcond` to
   calculate the reciprocal condition number.

   If the matrix is not positive-definite then the decomposition will fail, returning the
   error code :macro:`GSL_EDOM`. When testing whether a matrix is positive-definite, disable the error
   handler first to avoid triggering an error.

.. function:: int gsl_linalg_cholesky_band_solve (const gsl_matrix * LLT, const gsl_vector * b, gsl_vector * x)
              int gsl_linalg_cholesky_band_solvem (const gsl_matrix * LLT, const gsl_matrix * B, gsl_matrix * X)

   This function solves the symmetric banded system :math:`A x = b` (or :math:`A X = B`) using the Cholesky
   decomposition of :math:`A` held in the matrix :data:`LLT` which must
   have been previously computed by :func:`gsl_linalg_cholesky_band_decomp`.

.. function:: int gsl_linalg_cholesky_band_svx (const gsl_matrix * LLT, gsl_vector * x)
              int gsl_linalg_cholesky_band_svxm (const gsl_matrix * LLT, gsl_matrix * X)

   This function solves the symmetric banded system :math:`A x = b` (or :math:`A X = B`) in-place using the
   Cholesky decomposition of :math:`A` held in the matrix :data:`LLT`
   which must have been previously computed by :func:`gsl_linalg_cholesky_band_decomp`.
   On input :data:`x` (or :data:`X`) should contain the right-hand side :math:`b` (or :math:`B`),
   which is replaced by the solution on output.

.. function:: int gsl_linalg_cholesky_band_invert (const gsl_matrix * LLT, gsl_matrix * Ainv)

   This function computes the inverse of a symmetric banded matrix from its Cholesky
   decomposition :data:`LLT`, which must have been previously computed
   by :func:`gsl_linalg_cholesky_band_decomp`. On output, the inverse is
   stored in :data:`Ainv`, using both the lower and upper portions.

.. function:: int gsl_linalg_cholesky_band_unpack (const gsl_matrix * LLT, gsl_matrix * L)

   This function unpacks the lower triangular Cholesky factor from :data:`LLT` and stores
   it in the lower triangular portion of the :math:`N`-by-:math:`N` matrix :data:`L`. The
   upper triangular portion of :data:`L` is not referenced.

.. function:: int gsl_linalg_cholesky_band_scale (const gsl_matrix * A, gsl_vector * S)

   This function calculates a diagonal scaling transformation of the
   symmetric, positive definite banded matrix :data:`A`, such that
   :math:`S A S` has a condition number within a factor of :math:`N`
   of the matrix of smallest possible condition number over all
   possible diagonal scalings. On output, :data:`S` contains the
   scale factors, given by :math:`S_i = 1/\sqrt{A_{ii}}`.
   For any :math:`A_{ii} \le 0`, the corresponding scale factor :math:`S_i`
   is set to :math:`1`.

.. function:: int gsl_linalg_cholesky_band_scale_apply (gsl_matrix * A, const gsl_vector * S)

   This function applies the scaling transformation :data:`S` to the banded symmetric
   positive definite matrix :data:`A`. On output,
   :data:`A` is replaced by :math:`S A S`.

.. function:: int gsl_linalg_cholesky_band_rcond (const gsl_matrix * LLT, double * rcond, gsl_vector * work)

   This function estimates the reciprocal condition number (using the 1-norm) of the symmetric banded positive
   definite matrix :math:`A`, using its Cholesky decomposition provided in :data:`LLT`.
   The reciprocal condition number estimate, defined as :math:`1 / (||A||_1 \cdot ||A^{-1}||_1)`, is stored
   in :data:`rcond`. Additional workspace of size :math:`3 N` is required in :data:`work`.

.. index::
   single: banded LDLT decomposition
   single: LDLT decomposition, banded

Banded LDLT Decomposition
-------------------------

The routines in this section are designed to factor and solve
:math:`N`-by-:math:`N` linear systems of the form :math:`A x = b`
where :math:`A` is a banded, symmetric, and non-singular matrix
with lower bandwidth :math:`p`. See :ref:`sec_ldlt-decomposition`
for more information on the factorization. The lower triangular
factor of the :math:`L D L^T` decomposition preserves the same
banded structure as the matrix :math:`A`, enabling an efficient
algorithm which overwrites the original matrix with the
:math:`L` and :math:`D` factors.

.. function:: int gsl_linalg_ldlt_band_decomp (gsl_matrix * A)

   This function factorizes the symmetric, non-singular square matrix
   :data:`A` into the decomposition :math:`A = L D L^T`. The input matrix
   :data:`A` is given in :ref:`symmetric banded format <sec_symmetric-banded>`, and has dimensions
   :math:`N`-by-:math:`(p + 1)`, where :math:`p` is the lower bandwidth of the matrix.
   On output, the entries of :data:`A` are replaced by the entries of the matrices
   :math:`D` and :math:`L` in the same format.

   If the matrix is singular then the decomposition will fail, returning the
   error code :macro:`GSL_EDOM`.

.. function:: int gsl_linalg_ldlt_band_solve (const gsl_matrix * LDLT, const gsl_vector * b, gsl_vector * x)

   This function solves the symmetric banded system :math:`A x = b` using the
   :math:`L D L^T` decomposition of :math:`A` held in the matrix :data:`LDLT` which must
   have been previously computed by :func:`gsl_linalg_ldlt_band_decomp`.

.. function:: int gsl_linalg_ldlt_band_svx (const gsl_matrix * LDLT, gsl_vector * x)

   This function solves the symmetric banded system :math:`A x = b` in-place using the
   :math:`L D L^T` decomposition of :math:`A` held in the matrix :data:`LDLT`
   which must have been previously computed by :func:`gsl_linalg_ldlt_band_decomp`.
   On input :data:`x` should contain the right-hand side :math:`b`, which is replaced by the
   solution on output.

.. function:: int gsl_linalg_ldlt_band_unpack (const gsl_matrix * LDLT, gsl_matrix * L, gsl_vector * D)

   This function unpacks the unit lower triangular factor :math:`L` from :data:`LDLT` and stores
   it in the lower triangular portion of the :math:`N`-by-:math:`N` matrix :data:`L`. The
   upper triangular portion of :data:`L` is not referenced. The diagonal matrix
   :math:`D` is stored in the vector :data:`D`.

.. function:: int gsl_linalg_ldlt_band_rcond (const gsl_matrix * LDLT, double * rcond, gsl_vector * work)

   This function estimates the reciprocal condition number (using the 1-norm) of the symmetric banded
   nonsingular matrix :math:`A`, using its :math:`L D L^T` decomposition provided in :data:`LDLT`.
   The reciprocal condition number estimate, defined as :math:`1 / (||A||_1 \cdot ||A^{-1}||_1)`, is stored
   in :data:`rcond`. Additional workspace of size :math:`3 N` is required in :data:`work`.

.. index:: balancing matrices

.. _balancing:

Balancing
=========

The process of balancing a matrix applies similarity transformations
to make the rows and columns have comparable norms. This is
useful, for example, to reduce roundoff errors in the solution
of eigenvalue problems. Balancing a matrix :math:`A` consists
of replacing :math:`A` with a similar matrix

.. math:: A' = D^{-1} A D

where :math:`D` is a diagonal matrix whose entries are powers
of the floating point radix.

.. function:: int gsl_linalg_balance_matrix (gsl_matrix * A, gsl_vector * D)

   This function replaces the matrix :data:`A` with its balanced counterpart
   and stores the diagonal elements of the similarity transformation
   into the vector :data:`D`.

Examples
========

The following program solves the linear system :math:`A x = b`. The
system to be solved is,

.. only:: not texinfo

   .. math::

      \left(
      \begin{matrix}
        0.18& 0.60& 0.57& 0.96\\
        0.41& 0.24& 0.99& 0.58\\
        0.14& 0.30& 0.97& 0.66\\
        0.51& 0.13& 0.19& 0.85
      \end{matrix}
      \right)
      \left(
      \begin{matrix}
        x_0\\
        x_1\\
        x_2\\
        x_3
      \end{matrix}
      \right)
      =
      \left(
      \begin{matrix}
        1.0\\
        2.0\\
        3.0\\
        4.0
      \end{matrix}
      \right)

.. only:: texinfo

   ::

      [ 0.18 0.60 0.57 0.96 ] [x0]   [1.0]
      [ 0.41 0.24 0.99 0.58 ] [x1] = [2.0]
      [ 0.14 0.30 0.97 0.66 ] [x2]   [3.0]
      [ 0.51 0.13 0.19 0.85 ] [x3]   [4.0]

and the solution is found using LU decomposition of the matrix :math:`A`.

.. include:: examples/linalglu.c
   :code:

Here is the output from the program,

.. include:: examples/linalglu.txt
   :code:

This can be verified by multiplying the solution :math:`x` by the
original matrix :math:`A` using |octave|,

::

  octave> A = [ 0.18, 0.60, 0.57, 0.96;
                0.41, 0.24, 0.99, 0.58; 
                0.14, 0.30, 0.97, 0.66; 
                0.51, 0.13, 0.19, 0.85 ];

  octave> x = [ -4.05205; -12.6056; 1.66091; 8.69377];

  octave> A * x
  ans =
    1.0000
    2.0000
    3.0000
    4.0000

This reproduces the original right-hand side vector, :math:`b`, in
accordance with the equation :math:`A x = b`.

References and Further Reading
==============================

Further information on the algorithms described in this section can be
found in the following book,

* G. H. Golub, C. F. Van Loan, "Matrix Computations" (3rd Ed, 1996),
  Johns Hopkins University Press, ISBN 0-8018-5414-8.

The |lapack| library is described in the following manual,

* *LAPACK Users' Guide* (Third Edition, 1999), Published by SIAM,
  ISBN 0-89871-447-8

The |lapack| source code can be found at http://www.netlib.org/lapack,
along with an online copy of the users guide.

Further information on recursive Level 3 BLAS algorithms may be
found in the following paper,

* E. Peise and P. Bientinesi, "Recursive algorithms for dense linear algebra: the ReLAPACK collection",
  http://arxiv.org/abs/1602.06763, 2016.

The recursive Level 3 BLAS QR decomposition is described in the following paper,

* E. Elmroth and F. G. Gustavson, 2000. Applying recursion to serial and parallel QR
  factorization leads to better performance. IBM Journal of Research and Development,
  44(4), pp.605-624.

The Modified Golub-Reinsch algorithm is described in the following paper,

* T.F. Chan, "An Improved Algorithm for Computing the Singular Value
  Decomposition", ACM Transactions on Mathematical Software, 8
  (1982), pp 72--83.

The Jacobi algorithm for singular value decomposition is described in
the following papers,

* J.C. Nash, "A one-sided transformation method for the singular value
  decomposition and algebraic eigenproblem", Computer Journal,
  Volume 18, Number 1 (1975), p 74--76

* J.C. Nash and S. Shlien "Simple algorithms for the partial singular
  value decomposition", Computer Journal, Volume 30 (1987), p
  268--275.

* J. Demmel, K. Veselic, "Jacobi's Method is more accurate than
  QR", Lapack Working Note 15 (LAWN-15), October 1989. Available
  from netlib, http://www.netlib.org/lapack/ in the :code:`lawns` or
  :code:`lawnspdf` directories.

The algorithm for estimating a matrix condition number is described in
the following paper,

* N. J. Higham, "FORTRAN codes for estimating the one-norm of
  a real or complex matrix, with applications to condition estimation",
  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
