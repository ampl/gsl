.. index:: eigenvalues and eigenvectors

************
Eigensystems
************

.. include:: include.rst

This chapter describes functions for computing eigenvalues and
eigenvectors of matrices.  There are routines for real symmetric,
real nonsymmetric, complex hermitian, real generalized symmetric-definite,
complex generalized hermitian-definite, and real generalized nonsymmetric
eigensystems. Eigenvalues can be computed with or without eigenvectors.
The hermitian and real symmetric matrix algorithms are symmetric bidiagonalization
followed by QR reduction. The nonsymmetric algorithm is the Francis QR
double-shift.  The generalized nonsymmetric algorithm is the QZ method due
to Moler and Stewart.

The functions described in this chapter are declared in the header file
:file:`gsl_eigen.h`.

Real Symmetric Matrices
=======================
.. index::
   single: symmetric matrix, real, eigensystem
   single: real symmetric matrix, eigensystem

For real symmetric matrices, the library uses the symmetric
bidiagonalization and QR reduction method.  This is described in Golub
& van Loan, section 8.3.  The computed eigenvalues are accurate to an
absolute accuracy of :math:`\epsilon ||A||_2`, where :math:`\epsilon` is
the machine precision.

.. type:: gsl_eigen_symm_workspace

   This workspace contains internal parameters used for solving symmetric eigenvalue
   problems.

.. function:: gsl_eigen_symm_workspace * gsl_eigen_symm_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues of
   :data:`n`-by-:data:`n` real symmetric matrices.  The size of the workspace
   is :math:`O(2n)`.

.. function:: void gsl_eigen_symm_free (gsl_eigen_symm_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_symm (gsl_matrix * A, gsl_vector * eval, gsl_eigen_symm_workspace * w)

   This function computes the eigenvalues of the real symmetric matrix
   :data:`A`.  Additional workspace of the appropriate size must be provided
   in :data:`w`.  The diagonal and lower triangular part of :data:`A` are
   destroyed during the computation, but the strict upper triangular part
   is not referenced.  The eigenvalues are stored in the vector :data:`eval`
   and are unordered.

.. type:: gsl_eigen_symmv_workspace

   This workspace contains internal parameters used for solving symmetric eigenvalue
   and eigenvector problems.

.. function:: gsl_eigen_symmv_workspace * gsl_eigen_symmv_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues and
   eigenvectors of :data:`n`-by-:data:`n` real symmetric matrices.  The size of
   the workspace is :math:`O(4n)`.

.. function:: void gsl_eigen_symmv_free (gsl_eigen_symmv_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_symmv (gsl_matrix * A, gsl_vector * eval, gsl_matrix * evec, gsl_eigen_symmv_workspace * w)

   This function computes the eigenvalues and eigenvectors of the real
   symmetric matrix :data:`A`.  Additional workspace of the appropriate size
   must be provided in :data:`w`.  The diagonal and lower triangular part of
   :data:`A` are destroyed during the computation, but the strict upper
   triangular part is not referenced.  The eigenvalues are stored in the
   vector :data:`eval` and are unordered.  The corresponding eigenvectors are
   stored in the columns of the matrix :data:`evec`.  For example, the
   eigenvector in the first column corresponds to the first eigenvalue.
   The eigenvectors are guaranteed to be mutually orthogonal and normalised
   to unit magnitude.

Complex Hermitian Matrices
==========================

For hermitian matrices, the library uses the complex form of
the symmetric bidiagonalization and QR reduction method.

.. index::
   single: hermitian matrix, complex, eigensystem
   single: complex hermitian matrix, eigensystem

.. type:: gsl_eigen_herm_workspace

   This workspace contains internal parameters used for solving hermitian eigenvalue
   problems.

.. function:: gsl_eigen_herm_workspace * gsl_eigen_herm_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues of
   :data:`n`-by-:data:`n` complex hermitian matrices.  The size of the workspace
   is :math:`O(3n)`.

.. function:: void gsl_eigen_herm_free (gsl_eigen_herm_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_herm (gsl_matrix_complex * A, gsl_vector * eval, gsl_eigen_herm_workspace * w)

   This function computes the eigenvalues of the complex hermitian matrix
   :data:`A`.  Additional workspace of the appropriate size must be provided
   in :data:`w`.  The diagonal and lower triangular part of :data:`A` are
   destroyed during the computation, but the strict upper triangular part
   is not referenced.  The imaginary parts of the diagonal are assumed to be
   zero and are not referenced. The eigenvalues are stored in the vector
   :data:`eval` and are unordered.

.. type:: gsl_eigen_hermv_workspace

   This workspace contains internal parameters used for solving hermitian eigenvalue
   and eigenvector problems.

.. function:: gsl_eigen_hermv_workspace * gsl_eigen_hermv_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues and
   eigenvectors of :data:`n`-by-:data:`n` complex hermitian matrices.  The size of
   the workspace is :math:`O(5n)`.

.. function:: void gsl_eigen_hermv_free (gsl_eigen_hermv_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_hermv (gsl_matrix_complex * A, gsl_vector * eval, gsl_matrix_complex * evec, gsl_eigen_hermv_workspace * w)

   This function computes the eigenvalues and eigenvectors of the complex
   hermitian matrix :data:`A`.  Additional workspace of the appropriate size
   must be provided in :data:`w`.  The diagonal and lower triangular part of
   :data:`A` are destroyed during the computation, but the strict upper
   triangular part is not referenced. The imaginary parts of the diagonal
   are assumed to be zero and are not referenced.  The eigenvalues are
   stored in the vector :data:`eval` and are unordered.  The corresponding
   complex eigenvectors are stored in the columns of the matrix :data:`evec`.
   For example, the eigenvector in the first column corresponds to the
   first eigenvalue.  The eigenvectors are guaranteed to be mutually
   orthogonal and normalised to unit magnitude.

Real Nonsymmetric Matrices
==========================
.. index::
   single: nonsymmetric matrix, real, eigensystem
   single: real nonsymmetric matrix, eigensystem

The solution of the real nonsymmetric eigensystem problem for a
matrix :math:`A` involves computing the Schur decomposition

.. math:: A = Z T Z^T

where :math:`Z` is an orthogonal matrix of Schur vectors and :math:`T`,
the Schur form, is quasi upper triangular with diagonal
:math:`1`-by-:math:`1` blocks which are real eigenvalues of :math:`A`, and
diagonal :math:`2`-by-:math:`2` blocks whose eigenvalues are complex
conjugate eigenvalues of :math:`A`. The algorithm used is the double-shift 
Francis method.

.. type:: gsl_eigen_nonsymm_workspace

   This workspace contains internal parameters used for solving nonsymmetric eigenvalue
   problems.

.. function:: gsl_eigen_nonsymm_workspace * gsl_eigen_nonsymm_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues of
   :data:`n`-by-:data:`n` real nonsymmetric matrices. The size of the workspace
   is :math:`O(2n)`.

.. function:: void gsl_eigen_nonsymm_free (gsl_eigen_nonsymm_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: void gsl_eigen_nonsymm_params (const int compute_t, const int balance, gsl_eigen_nonsymm_workspace * w)

   This function sets some parameters which determine how the eigenvalue
   problem is solved in subsequent calls to :func:`gsl_eigen_nonsymm`.

   If :data:`compute_t` is set to 1, the full Schur form :math:`T` will be
   computed by :func:`gsl_eigen_nonsymm`. If it is set to 0,
   :math:`T` will not be computed (this is the default setting). Computing
   the full Schur form :math:`T` requires approximately 1.5--2 times the
   number of flops.

   If :data:`balance` is set to 1, a balancing transformation is applied
   to the matrix prior to computing eigenvalues. This transformation is
   designed to make the rows and columns of the matrix have comparable
   norms, and can result in more accurate eigenvalues for matrices
   whose entries vary widely in magnitude. See :ref:`Balancing <balancing>` for more
   information. Note that the balancing transformation does not preserve
   the orthogonality of the Schur vectors, so if you wish to compute the
   Schur vectors with :func:`gsl_eigen_nonsymm_Z` you will obtain the Schur
   vectors of the balanced matrix instead of the original matrix. The
   relationship will be

   .. math:: T = Q^T D^{-1} A D Q

   where :data:`Q` is the matrix of Schur vectors for the balanced matrix, and
   :data:`D` is the balancing transformation. Then :func:`gsl_eigen_nonsymm_Z`
   will compute a matrix :data:`Z` which satisfies

   .. math:: T = Z^{-1} A Z

   with :math:`Z = D Q`. Note that :data:`Z` will not be orthogonal. For
   this reason, balancing is not performed by default.

.. function:: int gsl_eigen_nonsymm (gsl_matrix * A, gsl_vector_complex * eval, gsl_eigen_nonsymm_workspace * w)

   This function computes the eigenvalues of the real nonsymmetric matrix
   :data:`A` and stores them in the vector :data:`eval`. If :math:`T` is
   desired, it is stored in the upper portion of :data:`A` on output.
   Otherwise, on output, the diagonal of :data:`A` will contain the
   :math:`1`-by-:math:`1` real eigenvalues and :math:`2`-by-:math:`2`
   complex conjugate eigenvalue systems, and the rest of :data:`A` is
   destroyed. In rare cases, this function may fail to find all
   eigenvalues. If this happens, an error code is returned
   and the number of converged eigenvalues is stored in :code:`w->n_evals`.
   The converged eigenvalues are stored in the beginning of :data:`eval`.

.. function:: int gsl_eigen_nonsymm_Z (gsl_matrix * A, gsl_vector_complex * eval, gsl_matrix * Z, gsl_eigen_nonsymm_workspace * w)

   This function is identical to :func:`gsl_eigen_nonsymm` except that it also
   computes the Schur vectors and stores them into :data:`Z`.

.. type:: gsl_eigen_nonsymmv_workspace

   This workspace contains internal parameters used for solving nonsymmetric eigenvalue
   and eigenvector problems.

.. function:: gsl_eigen_nonsymmv_workspace * gsl_eigen_nonsymmv_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues and
   eigenvectors of :data:`n`-by-:data:`n` real nonsymmetric matrices. The
   size of the workspace is :math:`O(5n)`.

.. function:: void gsl_eigen_nonsymmv_free (gsl_eigen_nonsymmv_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: void gsl_eigen_nonsymmv_params (const int balance, gsl_eigen_nonsymm_workspace * w)

   This function sets parameters which determine how the eigenvalue
   problem is solved in subsequent calls to :func:`gsl_eigen_nonsymmv`.
   If :data:`balance` is set to 1, a balancing transformation is applied
   to the matrix. See :func:`gsl_eigen_nonsymm_params` for more information.
   Balancing is turned off by default since it does not preserve the
   orthogonality of the Schur vectors.

.. function:: int gsl_eigen_nonsymmv (gsl_matrix * A, gsl_vector_complex * eval, gsl_matrix_complex * evec, gsl_eigen_nonsymmv_workspace * w)

   This function computes eigenvalues and right eigenvectors of the
   :data:`n`-by-:data:`n` real nonsymmetric matrix :data:`A`. It first calls
   :func:`gsl_eigen_nonsymm` to compute the eigenvalues, Schur form :math:`T`, and
   Schur vectors. Then it finds eigenvectors of :math:`T` and backtransforms
   them using the Schur vectors. The Schur vectors are destroyed in the
   process, but can be saved by using :func:`gsl_eigen_nonsymmv_Z`. The
   computed eigenvectors are normalized to have unit magnitude. On
   output, the upper portion of :data:`A` contains the Schur form
   :math:`T`. If :func:`gsl_eigen_nonsymm` fails, no eigenvectors are
   computed, and an error code is returned.

.. function:: int gsl_eigen_nonsymmv_Z (gsl_matrix * A, gsl_vector_complex * eval, gsl_matrix_complex * evec, gsl_matrix * Z, gsl_eigen_nonsymmv_workspace * w)

   This function is identical to :func:`gsl_eigen_nonsymmv` except that it also saves
   the Schur vectors into :data:`Z`.

Real Generalized Symmetric-Definite Eigensystems
================================================
.. index:: generalized symmetric eigensystems

The real generalized symmetric-definite eigenvalue problem is to find
eigenvalues :math:`\lambda` and eigenvectors :math:`x` such that

.. math:: A x = \lambda B x

where :math:`A` and :math:`B` are symmetric matrices, and :math:`B` is
positive-definite. This problem reduces to the standard symmetric
eigenvalue problem by applying the Cholesky decomposition to :math:`B`:

.. only:: not texinfo

   .. math::

      A x & = \lambda B x \\
      A x & = \lambda L L^T x \\
      \left( L^{-1} A L^{-T} \right) L^T x & = \lambda L^T x

.. only:: texinfo

   ::

      A x = \lambda B x
      A x = \lambda L L^T x
      ( L^{-1} A L^{-T} ) L^T x = \lambda L^T x

Therefore, the problem becomes :math:`C y = \lambda y` where
:math:`C = L^{-1} A L^{-T}`
is symmetric, and :math:`y = L^T x`. The standard
symmetric eigensolver can be applied to the matrix :math:`C`.
The resulting eigenvectors are backtransformed to find the
vectors of the original problem. The eigenvalues and eigenvectors
of the generalized symmetric-definite eigenproblem are always real.

.. type:: gsl_eigen_gensymm_workspace

   This workspace contains internal parameters used for solving generalized symmetric eigenvalue
   problems.

.. function:: gsl_eigen_gensymm_workspace * gsl_eigen_gensymm_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues of
   :data:`n`-by-:data:`n` real generalized symmetric-definite eigensystems. The
   size of the workspace is :math:`O(2n)`.

.. function:: void gsl_eigen_gensymm_free (gsl_eigen_gensymm_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_gensymm (gsl_matrix * A, gsl_matrix * B, gsl_vector * eval, gsl_eigen_gensymm_workspace * w)

   This function computes the eigenvalues of the real generalized
   symmetric-definite matrix pair (:data:`A`, :data:`B`), and stores them 
   in :data:`eval`, using the method outlined above. On output, :data:`B`
   contains its Cholesky decomposition and :data:`A` is destroyed.

.. type:: gsl_eigen_gensymmv_workspace

   This workspace contains internal parameters used for solving generalized symmetric eigenvalue
   and eigenvector problems.

.. function:: gsl_eigen_gensymmv_workspace * gsl_eigen_gensymmv_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues and
   eigenvectors of :data:`n`-by-:data:`n` real generalized symmetric-definite
   eigensystems. The size of the workspace is :math:`O(4n)`.

.. function:: void gsl_eigen_gensymmv_free (gsl_eigen_gensymmv_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_gensymmv (gsl_matrix * A, gsl_matrix * B, gsl_vector * eval, gsl_matrix * evec, gsl_eigen_gensymmv_workspace * w)

   This function computes the eigenvalues and eigenvectors of the real
   generalized symmetric-definite matrix pair (:data:`A`, :data:`B`), and
   stores them in :data:`eval` and :data:`evec` respectively. The computed
   eigenvectors are normalized to have unit magnitude. On output,
   :data:`B` contains its Cholesky decomposition and :data:`A` is destroyed.

Complex Generalized Hermitian-Definite Eigensystems
===================================================
.. index:: generalized hermitian definite eigensystems

The complex generalized hermitian-definite eigenvalue problem is to find
eigenvalues :math:`\lambda` and eigenvectors :math:`x` such that

.. math:: A x = \lambda B x

where :math:`A` and :math:`B` are hermitian matrices, and :math:`B` is
positive-definite. Similarly to the real case, this can be reduced
to :math:`C y = \lambda y` where
:math:`C = L^{-1} A L^{-\dagger}`
is hermitian, and
:math:`y = L^{\dagger} x`.  The standard
hermitian eigensolver can be applied to the matrix :math:`C`.
The resulting eigenvectors are backtransformed to find the
vectors of the original problem. The eigenvalues
of the generalized hermitian-definite eigenproblem are always real.

.. type:: gsl_eigen_genherm_workspace

   This workspace contains internal parameters used for solving generalized hermitian eigenvalue
   problems.

.. function:: gsl_eigen_genherm_workspace * gsl_eigen_genherm_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues of
   :data:`n`-by-:data:`n` complex generalized hermitian-definite eigensystems. The
   size of the workspace is :math:`O(3n)`.

.. function:: void gsl_eigen_genherm_free (gsl_eigen_genherm_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_genherm (gsl_matrix_complex * A, gsl_matrix_complex * B, gsl_vector * eval, gsl_eigen_genherm_workspace * w)

   This function computes the eigenvalues of the complex generalized
   hermitian-definite matrix pair (:data:`A`, :data:`B`), and stores them 
   in :data:`eval`, using the method outlined above. On output, :data:`B`
   contains its Cholesky decomposition and :data:`A` is destroyed.

.. type:: gsl_eigen_genhermv_workspace

   This workspace contains internal parameters used for solving generalized hermitian eigenvalue
   and eigenvector problems.

.. function:: gsl_eigen_genhermv_workspace * gsl_eigen_genhermv_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues and
   eigenvectors of :data:`n`-by-:data:`n` complex generalized hermitian-definite
   eigensystems. The size of the workspace is :math:`O(5n)`.

.. function:: void gsl_eigen_genhermv_free (gsl_eigen_genhermv_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_genhermv (gsl_matrix_complex * A, gsl_matrix_complex * B, gsl_vector * eval, gsl_matrix_complex * evec, gsl_eigen_genhermv_workspace * w)

   This function computes the eigenvalues and eigenvectors of the complex
   generalized hermitian-definite matrix pair (:data:`A`, :data:`B`), and
   stores them in :data:`eval` and :data:`evec` respectively. The computed
   eigenvectors are normalized to have unit magnitude. On output,
   :data:`B` contains its Cholesky decomposition and :data:`A` is destroyed.

Real Generalized Nonsymmetric Eigensystems
==========================================
.. index:: generalized eigensystems

Given two square matrices (:math:`A`, :math:`B`), the generalized
nonsymmetric eigenvalue problem is to find eigenvalues :math:`\lambda` and
eigenvectors :math:`x` such that

.. math:: A x = \lambda B x

We may also define the problem as finding eigenvalues :math:`\mu` and
eigenvectors :math:`y` such that

.. math:: \mu A y = B y

Note that these two problems are equivalent (with :math:`\lambda = 1/\mu`)
if neither :math:`\lambda` nor :math:`\mu` is zero. If say, :math:`\lambda`
is zero, then it is still a well defined eigenproblem, but its alternate
problem involving :math:`\mu` is not. Therefore, to allow for zero
(and infinite) eigenvalues, the problem which is actually solved is

.. math:: \beta A x = \alpha B x

The eigensolver routines below will return two values :math:`\alpha`
and :math:`\beta` and leave it to the user to perform the divisions
:math:`\lambda = \alpha / \beta` and :math:`\mu = \beta / \alpha`.

If the determinant of the matrix pencil :math:`A - \lambda B` is zero
for all :math:`\lambda`, the problem is said to be singular; otherwise
it is called regular.  Singularity normally leads to some
:math:`\alpha = \beta = 0` which means the eigenproblem is ill-conditioned
and generally does not have well defined eigenvalue solutions. The
routines below are intended for regular matrix pencils and could yield
unpredictable results when applied to singular pencils.

The solution of the real generalized nonsymmetric eigensystem problem for a
matrix pair :math:`(A, B)` involves computing the generalized Schur
decomposition

.. only:: not texinfo

   .. math::

      A &= Q S Z^T \\
      B &= Q T Z^T

.. only:: texinfo

   ::

      A = Q S Z^T
      B = Q T Z^T

where :math:`Q` and :math:`Z` are orthogonal matrices of left and right
Schur vectors respectively, and :math:`(S, T)` is the generalized Schur
form whose diagonal elements give the :math:`\alpha` and :math:`\beta`
values. The algorithm used is the QZ method due to Moler and Stewart
(see references).

.. type:: gsl_eigen_gen_workspace

   This workspace contains internal parameters used for solving generalized eigenvalue
   problems.

.. function:: gsl_eigen_gen_workspace * gsl_eigen_gen_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues of
   :data:`n`-by-:data:`n` real generalized nonsymmetric eigensystems. The
   size of the workspace is :math:`O(n)`.

.. function:: void gsl_eigen_gen_free (gsl_eigen_gen_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: void gsl_eigen_gen_params (const int compute_s, const int compute_t, const int balance, gsl_eigen_gen_workspace * w)

   This function sets some parameters which determine how the eigenvalue
   problem is solved in subsequent calls to :func:`gsl_eigen_gen`.

   If :data:`compute_s` is set to 1, the full Schur form :math:`S` will be
   computed by :func:`gsl_eigen_gen`. If it is set to 0,
   :math:`S` will not be computed (this is the default setting). :math:`S`
   is a quasi upper triangular matrix with 1-by-1 and 2-by-2 blocks
   on its diagonal. 1-by-1 blocks correspond to real eigenvalues, and
   2-by-2 blocks correspond to complex eigenvalues.

   If :data:`compute_t` is set to 1, the full Schur form :math:`T` will be
   computed by :func:`gsl_eigen_gen`. If it is set to 0,
   :math:`T` will not be computed (this is the default setting). :math:`T`
   is an upper triangular matrix with non-negative elements on its diagonal.
   Any 2-by-2 blocks in :math:`S` will correspond to a 2-by-2 diagonal
   block in :math:`T`.

   The :data:`balance` parameter is currently ignored, since generalized
   balancing is not yet implemented.

.. function:: int gsl_eigen_gen (gsl_matrix * A, gsl_matrix * B, gsl_vector_complex * alpha, gsl_vector * beta, gsl_eigen_gen_workspace * w)

   This function computes the eigenvalues of the real generalized nonsymmetric
   matrix pair (:data:`A`, :data:`B`), and stores them as pairs in
   (:data:`alpha`, :data:`beta`), where :data:`alpha` is complex and :data:`beta` is
   real. If :math:`\beta_i` is non-zero, then
   :math:`\lambda = \alpha_i / \beta_i` is an eigenvalue. Likewise,
   if :math:`\alpha_i` is non-zero, then
   :math:`\mu = \beta_i / \alpha_i` is an eigenvalue of the alternate
   problem :math:`\mu A y = B y`. The elements of :data:`beta` are normalized
   to be non-negative.

   If :math:`S` is desired, it is stored in :data:`A` on output. If :math:`T`
   is desired, it is stored in :data:`B` on output. The ordering of
   eigenvalues in (:data:`alpha`, :data:`beta`) follows the ordering
   of the diagonal blocks in the Schur forms :math:`S` and :math:`T`. In rare
   cases, this function may fail to find all eigenvalues. If this occurs, an
   error code is returned.

.. function:: int gsl_eigen_gen_QZ (gsl_matrix * A, gsl_matrix * B, gsl_vector_complex * alpha, gsl_vector * beta, gsl_matrix * Q, gsl_matrix * Z, gsl_eigen_gen_workspace * w)

   This function is identical to :func:`gsl_eigen_gen` except that it also
   computes the left and right Schur vectors and stores them into :data:`Q`
   and :data:`Z` respectively.

.. type:: gsl_eigen_genv_workspace

   This workspace contains internal parameters used for solving generalized eigenvalue
   and eigenvector problems.

.. function:: gsl_eigen_genv_workspace * gsl_eigen_genv_alloc (const size_t n)

   This function allocates a workspace for computing eigenvalues and
   eigenvectors of :data:`n`-by-:data:`n` real generalized nonsymmetric
   eigensystems. The size of the workspace is :math:`O(7n)`.

.. function:: void gsl_eigen_genv_free (gsl_eigen_genv_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_eigen_genv (gsl_matrix * A, gsl_matrix * B, gsl_vector_complex * alpha, gsl_vector * beta, gsl_matrix_complex * evec, gsl_eigen_genv_workspace * w)

   This function computes eigenvalues and right eigenvectors of the
   :data:`n`-by-:data:`n` real generalized nonsymmetric matrix pair
   (:data:`A`, :data:`B`). The eigenvalues are stored in (:data:`alpha`, :data:`beta`)
   and the eigenvectors are stored in :data:`evec`. It first calls
   :func:`gsl_eigen_gen` to compute the eigenvalues, Schur forms, and
   Schur vectors. Then it finds eigenvectors of the Schur forms and
   backtransforms them using the Schur vectors. The Schur vectors are
   destroyed in the process, but can be saved by using
   :func:`gsl_eigen_genv_QZ`. The computed eigenvectors are normalized
   to have unit magnitude. On output, (:data:`A`, :data:`B`) contains
   the generalized Schur form (:math:`S`, :math:`T`). If :func:`gsl_eigen_gen`
   fails, no eigenvectors are computed, and an error code is returned.

.. function:: int gsl_eigen_genv_QZ (gsl_matrix * A, gsl_matrix * B, gsl_vector_complex * alpha, gsl_vector * beta, gsl_matrix_complex * evec, gsl_matrix * Q, gsl_matrix * Z, gsl_eigen_genv_workspace * w)

   This function is identical to :func:`gsl_eigen_genv` except that it also
   computes the left and right Schur vectors and stores them into :data:`Q`
   and :data:`Z` respectively.

Sorting Eigenvalues and Eigenvectors
====================================
.. index:: sorting eigenvalues and eigenvectors

.. function:: int gsl_eigen_symmv_sort (gsl_vector * eval, gsl_matrix * evec, gsl_eigen_sort_t sort_type)

   This function simultaneously sorts the eigenvalues stored in the vector
   :data:`eval` and the corresponding real eigenvectors stored in the columns
   of the matrix :data:`evec` into ascending or descending order according to
   the value of the parameter :data:`sort_type`,

   .. type:: gsl_eigen_sort_t

      ================================= ===================================
      :macro:`GSL_EIGEN_SORT_VAL_ASC`   ascending order in numerical value
      :macro:`GSL_EIGEN_SORT_VAL_DESC`  descending order in numerical value
      :macro:`GSL_EIGEN_SORT_ABS_ASC`   ascending order in magnitude
      :macro:`GSL_EIGEN_SORT_ABS_DESC`  descending order in magnitude
      ================================= ===================================

.. function:: int gsl_eigen_hermv_sort (gsl_vector * eval, gsl_matrix_complex * evec, gsl_eigen_sort_t sort_type)

   This function simultaneously sorts the eigenvalues stored in the vector
   :data:`eval` and the corresponding complex eigenvectors stored in the
   columns of the matrix :data:`evec` into ascending or descending order
   according to the value of the parameter :data:`sort_type` as shown above.

.. function:: int gsl_eigen_nonsymmv_sort (gsl_vector_complex * eval, gsl_matrix_complex * evec, gsl_eigen_sort_t sort_type)

   This function simultaneously sorts the eigenvalues stored in the vector
   :data:`eval` and the corresponding complex eigenvectors stored in the
   columns of the matrix :data:`evec` into ascending or descending order
   according to the value of the parameter :data:`sort_type` as shown above.
   Only :macro:`GSL_EIGEN_SORT_ABS_ASC` and :macro:`GSL_EIGEN_SORT_ABS_DESC` are
   supported due to the eigenvalues being complex.

.. function:: int gsl_eigen_gensymmv_sort (gsl_vector * eval, gsl_matrix * evec, gsl_eigen_sort_t sort_type)

   This function simultaneously sorts the eigenvalues stored in the vector
   :data:`eval` and the corresponding real eigenvectors stored in the columns
   of the matrix :data:`evec` into ascending or descending order according to
   the value of the parameter :data:`sort_type` as shown above.

.. function:: int gsl_eigen_genhermv_sort (gsl_vector * eval, gsl_matrix_complex * evec, gsl_eigen_sort_t sort_type)

   This function simultaneously sorts the eigenvalues stored in the vector
   :data:`eval` and the corresponding complex eigenvectors stored in the columns
   of the matrix :data:`evec` into ascending or descending order according to
   the value of the parameter :data:`sort_type` as shown above.

.. function:: int gsl_eigen_genv_sort (gsl_vector_complex * alpha, gsl_vector * beta, gsl_matrix_complex * evec, gsl_eigen_sort_t sort_type)

   This function simultaneously sorts the eigenvalues stored in the vectors
   (:data:`alpha`, :data:`beta`) and the corresponding complex eigenvectors
   stored in the columns of the matrix :data:`evec` into ascending or
   descending order according to the value of the parameter :data:`sort_type`
   as shown above. Only :macro:`GSL_EIGEN_SORT_ABS_ASC` and
   :macro:`GSL_EIGEN_SORT_ABS_DESC` are supported due to the eigenvalues being
   complex.

Examples
========

The following program computes the eigenvalues and eigenvectors of the 4-th order Hilbert matrix, :math:`H(i,j) = 1/(i + j + 1)`.

.. include:: examples/eigen.c
   :code:

Here is the beginning of the output from the program::

   $ ./a.out 
   eigenvalue = 9.67023e-05
   eigenvector = 
   -0.0291933
   0.328712
   -0.791411
   0.514553
   ...

This can be compared with the corresponding output from |octave|::

   octave> [v,d] = eig(hilb(4));
   octave> diag(d)  
   ans =
   
      9.6702e-05
      6.7383e-03
      1.6914e-01
      1.5002e+00

   octave> v 
   v =
   
      0.029193   0.179186  -0.582076   0.792608
     -0.328712  -0.741918   0.370502   0.451923
      0.791411   0.100228   0.509579   0.322416
     -0.514553   0.638283   0.514048   0.252161

Note that the eigenvectors can differ by a change of sign, since the
sign of an eigenvector is arbitrary.

The following program illustrates the use of the nonsymmetric
eigensolver, by computing the eigenvalues and eigenvectors of
the Vandermonde matrix
:math:`V(x;i,j) = x_i^{n - j}`
with :math:`x = (-1,-2,3,4)`.

.. include:: examples/eigen_nonsymm.c
   :code:

Here is the beginning of the output from the program::

   $ ./a.out 
   eigenvalue = -6.41391 + 0i
   eigenvector = 
   -0.0998822 + 0i
   -0.111251 + 0i
   0.292501 + 0i
   0.944505 + 0i
   eigenvalue = 5.54555 + 3.08545i
   eigenvector = 
   -0.043487 + -0.0076308i
   0.0642377 + -0.142127i
   -0.515253 + 0.0405118i
   -0.840592 + -0.00148565i
   ...

This can be compared with the corresponding output from |octave|::

   octave> [v,d] = eig(vander([-1 -2 3 4]));
   octave> diag(d)
   ans =
   
     -6.4139 + 0.0000i
      5.5456 + 3.0854i
      5.5456 - 3.0854i
      2.3228 + 0.0000i
   
   octave> v
   v =
   
    Columns 1 through 3:
   
     -0.09988 + 0.00000i  -0.04350 - 0.00755i  -0.04350 + 0.00755i
     -0.11125 + 0.00000i   0.06399 - 0.14224i   0.06399 + 0.14224i
      0.29250 + 0.00000i  -0.51518 + 0.04142i  -0.51518 - 0.04142i
      0.94451 + 0.00000i  -0.84059 + 0.00000i  -0.84059 - 0.00000i
   
    Column 4:
   
     -0.14493 + 0.00000i
      0.35660 + 0.00000i
      0.91937 + 0.00000i
      0.08118 + 0.00000i

Note that the eigenvectors corresponding to the eigenvalue
:math:`5.54555 + 3.08545i` differ by the multiplicative constant
:math:`0.9999984 + 0.0017674i` which is an arbitrary phase factor 
of magnitude 1.

References and Further Reading
==============================

Further information on the algorithms described in this section can be
found in the following book,

* G. H. Golub, C. F. Van Loan, "Matrix Computations" (3rd Ed, 1996),
  Johns Hopkins University Press, ISBN 0-8018-5414-8.

Further information on the generalized eigensystems QZ algorithm
can be found in this paper,

* C. Moler, G. Stewart, "An Algorithm for Generalized Matrix Eigenvalue
  Problems", SIAM J. Numer. Anal., Vol 10, No 2, 1973.

.. index:: LAPACK

Eigensystem routines for very large matrices can be found in the
Fortran library |lapack|. The |lapack| library is described in,

* LAPACK Users' Guide (Third Edition, 1999), Published by SIAM,
  ISBN 0-89871-447-8.

The |lapack| source code can be found at the website http://www.netlib.org/lapack
along with an online copy of the users guide.
