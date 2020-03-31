/* L D L^T Decomposition
 *
 * Copyright (C) 2016 Patrick Alken
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 */

/*
 * L D L^T decomposition of a symmetric positive definite matrix.
 *
 * This algorithm does:
 *   P A P' = L D L'
 * with
 *   L  := unit lower left triangle matrix
 *   D  := diagonal matrix
 *   L' := the transposed form of L.
 *   P  := permutation matrix
 *
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>

#include "cholesky_common.c"

static double cholesky_LDLT_norm1(const gsl_matrix * LDLT, const gsl_permutation * p,
                                  gsl_vector * work);
static int cholesky_LDLT_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params);

typedef struct
{
  const gsl_matrix * LDLT;
  const gsl_permutation * perm;
} pcholesky_params;

/*
pcholesky_decomp()
  Perform Pivoted Cholesky LDLT decomposition of a symmetric positive
semidefinite matrix

Inputs: copy_uplo - copy lower triangle to upper to save original matrix
                    for rcond calculation later
        A         - (input) symmetric, positive semidefinite matrix,
                    stored in lower triangle
                    (output) lower triangle contains L; diagonal contains D
        p         - permutation vector

Return: success/error

Notes:
1) Based on algorithm 4.2.2 (Outer Product LDLT with Pivoting) of
Golub and Van Loan, Matrix Computations (4th ed).
*/

static int
pcholesky_decomp (const int copy_uplo, gsl_matrix * A, gsl_permutation * p)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR("LDLT decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (p->size != N)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_view diag = gsl_matrix_diagonal(A);
      size_t k;

      if (copy_uplo)
        {
          /* save a copy of A in upper triangle (for later rcond calculation) */
          gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, A, A);
        }

      gsl_permutation_init(p);

      for (k = 0; k < N; ++k)
        {
          gsl_vector_view w;
          size_t j;

          /* compute j = max_idx { A_kk, ..., A_nn } */
          w = gsl_vector_subvector(&diag.vector, k, N - k);
          j = gsl_vector_max_index(&w.vector) + k;
          gsl_permutation_swap(p, k, j);

          cholesky_swap_rowcol(A, k, j);

          if (k < N - 1)
            {
              double alpha = gsl_matrix_get(A, k, k);
              double alphainv = 1.0 / alpha;

              /* v = A(k+1:n, k) */
              gsl_vector_view v = gsl_matrix_subcolumn(A, k, k + 1, N - k - 1);

              /* m = A(k+1:n, k+1:n) */
              gsl_matrix_view m = gsl_matrix_submatrix(A, k + 1, k + 1, N - k - 1, N - k - 1);

              /* m = m - v v^T / alpha */
              gsl_blas_dsyr(CblasLower, -alphainv, &v.vector, &m.matrix);

              /* v = v / alpha */
              gsl_vector_scale(&v.vector, alphainv);
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_pcholesky_decomp()
  Perform Pivoted Cholesky LDLT decomposition of a symmetric positive
semidefinite matrix

Inputs: A - (input) symmetric, positive semidefinite matrix,
                    stored in lower triangle
            (output) lower triangle contains L; diagonal contains D
        p - permutation vector

Return: success/error

Notes:
1) Based on algorithm 4.2.2 (Outer Product LDLT with Pivoting) of
Golub and Van Loan, Matrix Computations (4th ed).
*/

int
gsl_linalg_pcholesky_decomp (gsl_matrix * A, gsl_permutation * p)
{
  int status = pcholesky_decomp(1, A, p);
  return status;
}

int
gsl_linalg_pcholesky_solve(const gsl_matrix * LDLT,
                           const gsl_permutation * p,
                           const gsl_vector * b,
                           gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (LDLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      gsl_vector_memcpy (x, b);

      status = gsl_linalg_pcholesky_svx (LDLT, p, x);
      
      return status;
    }
}

int
gsl_linalg_pcholesky_svx(const gsl_matrix * LDLT,
                         const gsl_permutation * p,
                         gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_const_view D = gsl_matrix_const_diagonal(LDLT);

      /* x := P b */
      gsl_permute_vector(p, x);

      /* solve: L w = P b */
      gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasUnit, LDLT, x);

      /* solve: D y = w */
      gsl_vector_div(x, &D.vector);

      /* solve: L^T z = y */
      gsl_blas_dtrsv(CblasLower, CblasTrans, CblasUnit, LDLT, x);

      /* compute: x = P^T z */
      gsl_permute_vector_inverse(p, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_pcholesky_decomp2(gsl_matrix * A, gsl_permutation * p,
                             gsl_vector * S)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (N != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (N != S->size)
    {
      GSL_ERROR("S must have length N", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* save a copy of A in upper triangle (for later rcond calculation) */
      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, A, A);

      /* compute scaling factors to reduce cond(A) */
      status = gsl_linalg_cholesky_scale(A, S);
      if (status)
        return status;

      /* apply scaling factors */
      status = gsl_linalg_cholesky_scale_apply(A, S);
      if (status)
        return status;

      /* compute Cholesky decomposition of diag(S) A diag(S) */
      status = pcholesky_decomp(0, A, p);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_pcholesky_solve2(const gsl_matrix * LDLT,
                            const gsl_permutation * p,
                            const gsl_vector * S,
                            const gsl_vector * b,
                            gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (LDLT->size1 != S->size)
    {
      GSL_ERROR ("matrix size must match S", GSL_EBADLEN);
    }
  else if (LDLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      gsl_vector_memcpy (x, b);

      status = gsl_linalg_pcholesky_svx2 (LDLT, p, S, x);
      
      return status;
    }
}

int
gsl_linalg_pcholesky_svx2(const gsl_matrix * LDLT,
                          const gsl_permutation * p,
                          const gsl_vector * S,
                          gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (LDLT->size1 != S->size)
    {
      GSL_ERROR ("matrix size must match S", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* x := S b */
      gsl_vector_mul(x, S);

      /* solve: A~ x~ = b~, with A~ = S A S, b~ = S b */
      status = gsl_linalg_pcholesky_svx(LDLT, p, x);
      if (status)
        return status;

      /* compute: x = S x~ */
      gsl_vector_mul(x, S);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_pcholesky_invert()
  Compute the inverse of a symmetric positive definite matrix in
Cholesky form.

Inputs: LDLT - matrix in cholesky form
        p    - permutation
        Ainv - (output) A^{-1}

Return: success or error
*/

int
gsl_linalg_pcholesky_invert(const gsl_matrix * LDLT, const gsl_permutation * p,
                            gsl_matrix * Ainv)
{
  const size_t M = LDLT->size1;
  const size_t N = LDLT->size2;

  if (M != N)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (Ainv->size1 != Ainv->size2)
    {
      GSL_ERROR ("Ainv matrix must be square", GSL_ENOTSQR);
    }
  else if (Ainv->size1 != M)
    {
      GSL_ERROR ("Ainv matrix has wrong dimensions", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* invert the lower triangle of LDLT */
      gsl_matrix_memcpy(Ainv, LDLT);
      gsl_linalg_tri_invert(CblasLower, CblasUnit, Ainv);

      /* compute sqrt(D^{-1}) L^{-1} in the lower triangle of Ainv */
      for (i = 0; i < N; ++i)
        {
          double di = gsl_matrix_get(LDLT, i, i);
          double invsqrt_di = 1.0 / sqrt(di);

          if (i > 0)
            {
              gsl_vector_view v = gsl_matrix_subrow(Ainv, i, 0, i);
              gsl_blas_dscal(invsqrt_di, &v.vector);
            }

          gsl_matrix_set(Ainv, i, i, invsqrt_di);
        }

      /* compute A^{-1} = L^{-T} D^{-1} L^{-1} */
      gsl_linalg_tri_LTL(Ainv);

      /* copy lower triangle to upper */
      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, Ainv, Ainv);

      /* now apply permutation p to the matrix */

      /* compute L^{-T} D^{-1} L^{-1} P^T */
      for (i = 0; i < N; ++i)
        {
          gsl_vector_view v = gsl_matrix_row(Ainv, i);
          gsl_permute_vector_inverse(p, &v.vector);
        }

      /* compute P L^{-T} D^{-1} L^{-1} P^T */
      for (i = 0; i < N; ++i)
        {
          gsl_vector_view v = gsl_matrix_column(Ainv, i);
          gsl_permute_vector_inverse(p, &v.vector);
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_pcholesky_rcond (const gsl_matrix * LDLT, const gsl_permutation * p,
                            double * rcond, gsl_vector * work)
{
  const size_t M = LDLT->size1;
  const size_t N = LDLT->size2;

  if (M != N)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (work->size != 3 * N)
    {
      GSL_ERROR ("work vector must have length 3*N", GSL_EBADLEN);
    }
  else
    {
      int status;
      double Anorm = cholesky_LDLT_norm1(LDLT, p, work); /* ||A||_1 */
      double Ainvnorm;                                   /* ||A^{-1}||_1 */
      pcholesky_params params;

      *rcond = 0.0;

      /* don't continue if matrix is singular */
      if (Anorm == 0.0)
        return GSL_SUCCESS;

      params.LDLT = LDLT;
      params.perm = p;

      /* estimate ||A^{-1}||_1 */
      status = gsl_linalg_invnorm1(N, cholesky_LDLT_Ainv, &params, &Ainvnorm, work);

      if (status)
        return status;

      if (Ainvnorm != 0.0)
        *rcond = (1.0 / Anorm) / Ainvnorm;

      return GSL_SUCCESS;
    }
}

/*
cholesky_LDLT_norm1
  Compute 1-norm of original matrix A, stored in upper triangle of LDLT;
diagonal entries have to be reconstructed

Inputs: LDLT - Cholesky L D L^T decomposition (lower triangle) with
               original matrix in upper triangle
        p    - permutation vector
        work - workspace, length 2*N
*/

static double
cholesky_LDLT_norm1(const gsl_matrix * LDLT, const gsl_permutation * p, gsl_vector * work)
{
  const size_t N = LDLT->size1;
  gsl_vector_const_view D = gsl_matrix_const_diagonal(LDLT);
  gsl_vector_view diagA = gsl_vector_subvector(work, N, N);
  double max = 0.0;
  size_t i, j;

  /* reconstruct diagonal entries of original matrix A */
  for (j = 0; j < N; ++j)
    {
      double Ajj;

      /* compute diagonal (j,j) entry of A */
      Ajj = gsl_vector_get(&D.vector, j);
      for (i = 0; i < j; ++i)
        {
          double Di = gsl_vector_get(&D.vector, i);
          double Lji = gsl_matrix_get(LDLT, j, i);

          Ajj += Di * Lji * Lji;
        }

      gsl_vector_set(&diagA.vector, j, Ajj);
    }

  gsl_permute_vector_inverse(p, &diagA.vector);

  for (j = 0; j < N; ++j)
    {
      double sum = 0.0;
      double Ajj = gsl_vector_get(&diagA.vector, j);

      for (i = 0; i < j; ++i)
        {
          double *wi = gsl_vector_ptr(work, i);
          double Aij = gsl_matrix_get(LDLT, i, j);
          double absAij = fabs(Aij);

          sum += absAij;
          *wi += absAij;
        }

      gsl_vector_set(work, j, sum + fabs(Ajj));
    }

  for (i = 0; i < N; ++i)
    {
      double wi = gsl_vector_get(work, i);
      max = GSL_MAX(max, wi);
    }

  return max;
}

/* x := A^{-1} x = A^{-t} x, A = L D L^T */
static int
cholesky_LDLT_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params)
{
  int status;
  pcholesky_params *par = (pcholesky_params *) params;

  (void) TransA; /* unused parameter warning */

  status = gsl_linalg_pcholesky_svx(par->LDLT, par->perm, x);

  return status;
}
