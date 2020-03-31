/* linalg/ldlt_band.c
 * 
 * Copyright (C) 2018 Patrick Alken
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

/* L D L^T decomposition of a symmetric banded positive semi-definite matrix */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static double symband_norm1(const gsl_matrix * A);
static int ldlt_band_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params);

/*
gsl_linalg_ldlt_band_decomp()
  L D L^T decomposition of a square symmetric positive semi-definite banded
matrix

Inputs: A - matrix in symmetric banded format, N-by-ndiag where N is the size of
            the matrix and ndiag is the number of nonzero diagonals.

Notes:
1) The matrix D is stored in the first column of A;
the first subdiagonal of L in the second column and so on.

2) If ndiag > 1, the 1-norm of A is stored in A(N,ndiag) on output

3) At each diagonal element, the matrix is factored as

A(j:end,j:end) = [ A11 A21^T ] = [ 1 0 ] [ alpha 0 ] [ 1 v^T ]
                 [ A21 A22   ]   [ v L ] [   0   D ] [ 0 L^T ]

where:

alpha = A(j,j)
v = A(j+1:end, j) / alpha
A22 = L D L^T + alpha v v^T

So we start at A(1,1) and work right. Pseudo-code is:

loop j = 1, ..., N
  alpha = A(j,j)
  A(j+1:end, j) := A(j+1:end, j) / alpha     (DSCAL)
  A(j+1:end, j+1:end) -= alpha v v^T         (DSYR)

Due to the banded structure, v has at most p non-zero elements, where
p is the lower bandwidth
*/

int
gsl_linalg_ldlt_band_decomp(gsl_matrix * A)
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
      double Anorm;
      size_t j;

      /* check for quick return */
      if (ndiag == 1)
        return GSL_SUCCESS;

      /*
       * calculate 1-norm of A and store in lower right of matrix, which is not accessed
       * by rest of routine. gsl_linalg_ldlt_band_rcond() will use this later. If
       * A is diagonal, there is no empty slot to store the 1-norm, so the rcond routine
       * will have to compute it.
       */
      Anorm = symband_norm1(A);
      gsl_matrix_set(A, N - 1, p, Anorm);

      for (j = 0; j < N - 1; ++j)
        {
          double ajj = gsl_matrix_get(A, j, 0);
          size_t lenv;

          if (ajj == 0.0)
            {
              GSL_ERROR("matrix is singular", GSL_EDOM);
            }

          /* number of elements in v, which will normally be p, unless we
           * are in lower right corner of matrix */
          lenv = GSL_MIN(p, N - j - 1);

          if (lenv > 0)
            {
              gsl_vector_view v = gsl_matrix_subrow(A, j, 1, lenv);
              gsl_matrix_view m = gsl_matrix_submatrix(A, j + 1, 0, lenv, lenv);

              gsl_blas_dscal(1.0 / ajj, &v.vector);

              m.matrix.tda = kld;
              gsl_blas_dsyr(CblasUpper, -ajj, &v.vector, &m.matrix);
            }
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_ldlt_band_solve (const gsl_matrix * LDLT,
                            const gsl_vector * b,
                            gsl_vector * x)
{
  if (LDLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LDLT->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* copy x <- b */
      gsl_vector_memcpy (x, b);

      status = gsl_linalg_ldlt_band_svx(LDLT, x);

      return status;
    }
}

int
gsl_linalg_ldlt_band_svx (const gsl_matrix * LDLT, gsl_vector * x)
{
  if (LDLT->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_const_view diag = gsl_matrix_const_column(LDLT, 0);

      /* solve for z using forward-substitution, L z = b */
      cblas_dtbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasUnit,
                  (int) LDLT->size1, (int) (LDLT->size2 - 1), LDLT->data, LDLT->tda,
                  x->data, x->stride);

      /* solve for y, D y = z */
      gsl_vector_div(x, &diag.vector);

      /* perform back-substitution, L^T x = y */
      cblas_dtbsv(CblasColMajor, CblasLower, CblasTrans, CblasUnit,
                  (int) LDLT->size1, (int) (LDLT->size2 - 1), LDLT->data, LDLT->tda,
                  x->data, x->stride);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_ldlt_band_unpack()
  Unpack symmetric banded format matrix LDLT into
larger matrix L and diagonal vector D
*/

int
gsl_linalg_ldlt_band_unpack (const gsl_matrix * LDLT, gsl_matrix * L, gsl_vector * D)
{
  const size_t N = LDLT->size1;

  if (N != L->size1)
    {
      GSL_ERROR("L matrix does not match LDLT dimensions", GSL_EBADLEN);
    }
  else if (L->size1 != L->size2)
    {
      GSL_ERROR("L matrix is not square", GSL_ENOTSQR);
    }
  else if (N != D->size)
    {
      GSL_ERROR("D vector does not match LDLT dimensions", GSL_EBADLEN);
    }
  else
    {
      const size_t p = LDLT->size2 - 1; /* lower bandwidth */
      gsl_vector_const_view diag = gsl_matrix_const_column(LDLT, 0);
      gsl_vector_view diagL = gsl_matrix_diagonal(L);
      size_t i;

      /* copy diagonal entries */
      gsl_vector_memcpy(D, &diag.vector);

      /* copy subdiagonals into L */
      for (i = 1; i <= p; ++i)
        {
          gsl_vector_const_view v = gsl_matrix_const_subcolumn(LDLT, i, 0, N - i);
          gsl_vector_view w = gsl_matrix_subdiagonal(L, i);
          gsl_vector_memcpy(&w.vector, &v.vector);
        }

      /* set main diagonal of L */
      gsl_vector_set_all(&diagL.vector, 1.0);

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
gsl_linalg_ldlt_band_rcond (const gsl_matrix * LDLT, double * rcond, gsl_vector * work)
{
  const size_t N = LDLT->size1;

  if (work->size != 3 * N)
    {
      GSL_ERROR ("work vector must have length 3*N", GSL_EBADLEN);
    }
  else
    {
      int status;
      const size_t ndiag = LDLT->size2;
      double Anorm;    /* ||A||_1 */
      double Ainvnorm; /* ||A^{-1}||_1 */

      if (ndiag == 1)
        {
          /* diagonal matrix, compute 1-norm since it has not been stored */
          Anorm = symband_norm1(LDLT);
        }
      else
        {
          /* 1-norm is stored in A(N, ndiag) by gsl_linalg_ldlt_band_decomp() */
          Anorm = gsl_matrix_get(LDLT, N - 1, ndiag - 1);
        }

      *rcond = 0.0;

      /* return if matrix is singular */
      if (Anorm == 0.0)
        return GSL_SUCCESS;

      status = gsl_linalg_invnorm1(N, ldlt_band_Ainv, (void *) LDLT, &Ainvnorm, work);
      if (status)
        return status;

      if (Ainvnorm != 0.0)
        *rcond = (1.0 / Anorm) / Ainvnorm;

      return GSL_SUCCESS;
    }
}

/* compute 1-norm of symmetric banded matrix */
static double
symband_norm1(const gsl_matrix * A)
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

/* x := A^{-1} x = A^{-t} x, A = L D L^T */
static int
ldlt_band_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params)
{
  gsl_matrix * LDLT = (gsl_matrix * ) params;
  gsl_vector_const_view diag = gsl_matrix_const_column(LDLT, 0);

  (void) TransA; /* unused parameter warning */

  /* compute x := L^{-1} x */
  cblas_dtbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasUnit,
              (int) LDLT->size1, (int) (LDLT->size2 - 1), LDLT->data, LDLT->tda,
              x->data, x->stride);

  /* compute x := D^{-1} x */
  gsl_vector_div(x, &diag.vector);

  /* compute x := L^{-T} x */
  cblas_dtbsv(CblasColMajor, CblasLower, CblasTrans, CblasUnit,
              (int) LDLT->size1, (int) (LDLT->size2 - 1), LDLT->data, LDLT->tda,
              x->data, x->stride);

  return GSL_SUCCESS;
}
