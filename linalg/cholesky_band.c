/* linalg/cholesky_band.c
 *
 * Copyright (C) 2018, 2019, 2020 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

static double cholesky_band_norm1(const gsl_matrix * A);
static int cholesky_band_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params);

/*
gsl_linalg_cholesky_band_decomp()
  Cholesky decomposition of a square symmetric positive definite banded
matrix

Inputs: A - matrix in banded format, N-by-ndiag where N is the size of
            the matrix and ndiag is the number of nonzero diagonals.

Notes:
1) The main diagonal of the Cholesky factor is stored in the first column
of A; the first subdiagonal in the second column and so on.

2) If ndiag > 1, the 1-norm of A is stored in A(N,ndiag) on output

3) At each diagonal element, the matrix is factored as

A(j:end,j:end) = [ A11 A21^T ] = [ alpha  0 ] [ alpha v^T ]
                 [ A21 A22   ]   [   v    L ] [   0   L^T ]

where:

alpha = sqrt(A(j,j))
v = A(j+1:end, j) / alpha
A22 = L L^T + v v^T

So we start at A(1,1) and work right. Pseudo-code is:

loop j = 1, ..., N
  alpha = sqrt(A(j,j))
  A(j+1:end, j) := A(j+1:end, j) / alpha     (DSCAL)
  A(j+1:end, j+1:end) -= v v^T               (DSYR)

Due to the banded structure, v has at most p non-zero elements, where
p is the lower bandwidth
*/

int
gsl_linalg_cholesky_band_decomp(gsl_matrix * A)
{
  const size_t N = A->size1;     /* size of matrix */
  const size_t ndiag = A->size2; /* number of diagonals in band, including main diagonal */

  if (ndiag > N)
    {
      GSL_ERROR ("invalid matrix dimensions", GSL_EBADLEN);
    }
  else
    {
      const size_t p = ndiag - 1; /* lower bandwidth */
      const int kld = (int) GSL_MAX(1, p);
      size_t j;

      if (ndiag > 1)
        {
          /*
           * calculate 1-norm of A and store in lower right of matrix, which is not accessed
           * by rest of routine. gsl_linalg_cholesky_band_rcond() will use this later. If
           * A is diagonal, there is no empty slot to store the 1-norm, so the rcond routine
           * will have to reconstruct it from the Cholesky factor.
           */
          double Anorm = cholesky_band_norm1(A);
          gsl_matrix_set(A, N - 1, p, Anorm);
        }

      for (j = 0; j < N; ++j)
        {
          double ajj = gsl_matrix_get(A, j, 0);
          size_t lenv;

          if (ajj <= 0.0)
            {
              GSL_ERROR("matrix is not positive definite", GSL_EDOM);
            }

          ajj = sqrt(ajj);
          gsl_matrix_set(A, j, 0, ajj);

          /* number of elements in v, which will normally be p, unless we
           * are in lower right corner of matrix */
          lenv = GSL_MIN(p, N - j - 1);

          if (lenv > 0)
            {
              gsl_vector_view v = gsl_matrix_subrow(A, j, 1, lenv);
              gsl_matrix_view m = gsl_matrix_submatrix(A, j + 1, 0, lenv, lenv);

              gsl_blas_dscal(1.0 / ajj, &v.vector);

              m.matrix.tda = kld;
              gsl_blas_dsyr(CblasUpper, -1.0, &v.vector, &m.matrix);
            }
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_band_solve (const gsl_matrix * LLT,
                                const gsl_vector * b,
                                gsl_vector * x)
{
  if (LLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LLT->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* copy x <- b */
      gsl_vector_memcpy (x, b);

      status = gsl_linalg_cholesky_band_svx(LLT, x);

      return status;
    }
}

int
gsl_linalg_cholesky_band_svx (const gsl_matrix * LLT, gsl_vector * x)
{
  if (LLT->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* solve for c using forward-substitution, L c = b */
      cblas_dtbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                  (int) LLT->size1, (int) (LLT->size2 - 1), LLT->data, LLT->tda,
                  x->data, x->stride);

      /* perform back-substitution, L^T x = c */
      cblas_dtbsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit,
                  (int) LLT->size1, (int) (LLT->size2 - 1), LLT->data, LLT->tda,
                  x->data, x->stride);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_band_solvem (const gsl_matrix * LLT,
                                 const gsl_matrix * B,
                                 gsl_matrix * X)
{
  if (LLT->size1 != B->size1)
    {
      GSL_ERROR ("LLT size1 must match B size1", GSL_EBADLEN);
    }
  else if (LLT->size1 != X->size1)
    {
      GSL_ERROR ("LLT size1 must match solution size1", GSL_EBADLEN);
    }
  else if (B->size2 != X->size2)
    {
      GSL_ERROR ("B size2 must match X size2", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* copy X <- B */
      gsl_matrix_memcpy (X, B);

      status = gsl_linalg_cholesky_band_svxm(LLT, X);

      return status;
    }
}

int
gsl_linalg_cholesky_band_svxm (const gsl_matrix * LLT, gsl_matrix * X)
{
  if (LLT->size1 != X->size1)
    {
      GSL_ERROR ("LLT size1 must match solution size1", GSL_EBADLEN);
    }
  else
    {
      int status;
      const size_t nrhs = X->size2;
      size_t j;

      for (j = 0; j < nrhs; ++j)
        {
          gsl_vector_view xj = gsl_matrix_column(X, j);

          status = gsl_linalg_cholesky_band_svx (LLT, &xj.vector);
          if (status)
            return status;
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_band_invert (const gsl_matrix * LLT, gsl_matrix * Ainv)
{
  if (Ainv->size1 != Ainv->size2)
    {
      GSL_ERROR("Ainv must be square", GSL_ENOTSQR);
    }
  else if (LLT->size1 != Ainv->size1)
    {
      GSL_ERROR("cholesky matrix has different dimensions from Ainv", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* unpack Cholesky factor into lower triangle of Ainv */
      status = gsl_linalg_cholesky_band_unpack(LLT, Ainv);
      if (status)
        return status;

      /* call the standard Cholesky inversion routine */
      status = gsl_linalg_cholesky_invert(Ainv);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_band_unpack (const gsl_matrix * LLT, gsl_matrix * L)
{
  const size_t N = LLT->size1;

  if (N != L->size1)
    {
      GSL_ERROR("L matrix does not match LLT dimensions", GSL_EBADLEN);
    }
  else if (L->size1 != L->size2)
    {
      GSL_ERROR("L matrix is not square", GSL_ENOTSQR);
    }
  else
    {
      const size_t p = LLT->size2 - 1; /* lower bandwidth */
      size_t i;

      for (i = 0; i < p + 1; ++i)
        {
          gsl_vector_const_view v = gsl_matrix_const_subcolumn(LLT, i, 0, N - i);
          gsl_vector_view w = gsl_matrix_subdiagonal(L, i);

          gsl_vector_memcpy(&w.vector, &v.vector);
        }

      /* zero out remaining subdiagonals */
      for (i = p + 1; i < N; ++i)
        {
          gsl_vector_view w = gsl_matrix_subdiagonal(L, i);
          gsl_vector_set_zero(&w.vector);
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_band_rcond (const gsl_matrix * LLT, double * rcond, gsl_vector * work)
{
  const size_t N = LLT->size1;

  if (work->size != 3 * N)
    {
      GSL_ERROR ("work vector must have length 3*N", GSL_EBADLEN);
    }
  else
    {
      int status;
      const size_t ndiag = LLT->size2;
      double Anorm;    /* ||A||_1 */
      double Ainvnorm; /* ||A^{-1}||_1 */

      if (ndiag == 1)
        {
          /* diagonal matrix, compute 1-norm since it has not been stored */
          gsl_vector_const_view v = gsl_matrix_const_column(LLT, 0);
          Anorm = gsl_vector_max(&v.vector);
          Anorm = Anorm * Anorm;
        }
      else
        {
          /* 1-norm is stored in A(N, ndiag) by gsl_linalg_cholesky_band_decomp() */
          Anorm = gsl_matrix_get(LLT, N - 1, ndiag - 1);
        }

      *rcond = 0.0;

      /* return if matrix is singular */
      if (Anorm == 0.0)
        return GSL_SUCCESS;

      status = gsl_linalg_invnorm1(N, cholesky_band_Ainv, (void *) LLT, &Ainvnorm, work);
      if (status)
        return status;

      if (Ainvnorm != 0.0)
        *rcond = (1.0 / Anorm) / Ainvnorm;

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_cholesky_band_scale()
  This function computes scale factors diag(S), such that

diag(S) A diag(S)

has a condition number within a factor N of the matrix
with the smallest condition number over all possible
diagonal scalings. See Corollary 7.6 of:

N. J. Higham, Accuracy and Stability of Numerical Algorithms (2nd Edition),
SIAM, 2002.

Inputs: A - symmetric positive definite matrix
        S - (output) scale factors, S_i = 1 / sqrt(A_ii)
*/

int
gsl_linalg_cholesky_band_scale(const gsl_matrix * A, gsl_vector * S)
{
  const size_t N = A->size1;     /* size of matrix */
  const size_t ndiag = A->size2; /* number of diagonals in band, including main diagonal */

  if (ndiag > N)
    {
      GSL_ERROR ("invalid matrix dimensions", GSL_EBADLEN);
    }
  else if (N != S->size)
    {
      GSL_ERROR("S must have length N", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* compute S_i = 1/sqrt(A_{ii}) */
      for (i = 0; i < N; ++i)
        {
          double Aii = gsl_matrix_get(A, i, 0);

          if (Aii <= 0.0)
            gsl_vector_set(S, i, 1.0); /* matrix not positive definite */
          else
            gsl_vector_set(S, i, 1.0 / sqrt(Aii));
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_cholesky_band_scale_apply()
  This function applies scale transformation to A:

A <- diag(S) A diag(S)

Inputs: A - (input/output)
            on input, symmetric positive definite matrix in banded format
            on output, diag(S) * A * diag(S) in banded format
        S - (input) scale factors
*/

int
gsl_linalg_cholesky_band_scale_apply(gsl_matrix * A, const gsl_vector * S)
{
  const size_t N = A->size1;     /* size of matrix */
  const size_t ndiag = A->size2; /* number of diagonals in band, including main diagonal */

  if (ndiag > N)
    {
      GSL_ERROR ("invalid matrix dimensions", GSL_EBADLEN);
    }
  else if (N != S->size)
    {
      GSL_ERROR("S must have length N", GSL_EBADLEN);
    }
  else
    {
      size_t i, j;

      for (j = 0; j < N; ++j)
        {
          double sj = gsl_vector_get(S, j);

          for (i = j; i < GSL_MIN(N, j + ndiag); ++i)
            {
              double si = gsl_vector_get(S, i);
              double * ptr = gsl_matrix_ptr(A, j, i - j);
              *ptr *= sj * si;
            }
        }

      return GSL_SUCCESS;
    }
}

/* compute 1-norm of symmetric banded matrix */
static double
cholesky_band_norm1(const gsl_matrix * A)
{
  const size_t N = A->size1;
  const size_t ndiag = A->size2; /* number of diagonals in band, including main diagonal */
  double value;

  if (ndiag == 1)
    {
      /* diagonal matrix */
      gsl_vector_const_view v = gsl_matrix_const_column(A, 0);
      CBLAS_INDEX_t idx = gsl_blas_idamax(&v.vector);
      value = gsl_vector_get(&v.vector, idx);
    }
  else
    {
      size_t j;

      value = 0.0;
      for (j = 0; j < N; ++j)
        {
          size_t ncol = GSL_MIN(ndiag, N - j); /* number of elements in column j below and including main diagonal */
          gsl_vector_const_view v = gsl_matrix_const_subrow(A, j, 0, ncol);
          double sum = gsl_blas_dasum(&v.vector);
          size_t k, l;

          /* sum now contains the absolute sum of elements below and including main diagonal for column j; we
           * have to add the symmetric elements above the diagonal */
          k = j;
          l = 1;
          while (k > 0 && l < ndiag)
            {
              double Akl = gsl_matrix_get(A, --k, l++);
              sum += fabs(Akl);
            }

          value = GSL_MAX(value, sum);
        }
    }

  return value;
}

/* x := A^{-1} x = A^{-t} x, A = L L^T */
static int
cholesky_band_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params)
{
  gsl_matrix * LLT = (gsl_matrix * ) params;

  (void) TransA; /* unused parameter warning */

  /* compute x := L^{-1} x */
  cblas_dtbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
              (int) LLT->size1, (int) (LLT->size2 - 1), LLT->data, LLT->tda,
              x->data, x->stride);

  /* compute x := L^{-T} x */
  cblas_dtbsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit,
              (int) LLT->size1, (int) (LLT->size2 - 1), LLT->data, LLT->tda,
              x->data, x->stride);

  return GSL_SUCCESS;
}
