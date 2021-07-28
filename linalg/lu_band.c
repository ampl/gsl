/* linalg/lu.c
 * 
 * Copyright (C) 2020 Patrick Alken
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
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static int LU_band_decomp_L2 (const size_t M, const size_t lb, const size_t ub,
                              gsl_matrix * AB, gsl_vector_uint * ipiv);

/* Factorise a (p,q) banded M x N matrix A into,
 *
 *   P A = L U
 *
 * where P is a permutation matrix, L is unit lower triangular and U
 * is upper triangular.
 *
 * L is stored in the strict lower triangular part of the input
 * matrix. The diagonal elements of L are unity and are not stored.
 *
 * U is stored in the diagonal and upper triangular part of the
 * input matrix.  
 * 
 * P is stored in the permutation p. Column j of P is column k of the
 * identity matrix, where k = permutation->data[j]
 *
 * Inputs: M   - number of rows in matrix
 *         lb  - lower bandwidth
 *         ub  - upper bandwidth
 *         AB  - matrix in band storage format, N-by-(2*p + q + 1)
 *         piv - pivot vector, size MIN(M, N)
 */

int
gsl_linalg_LU_band_decomp (const size_t M, const size_t lb, const size_t ub, gsl_matrix * AB, gsl_vector_uint * piv)
{
  const size_t N = AB->size1;
  const size_t minMN = GSL_MIN(M, N);

  if (lb >= M)
    {
      GSL_ERROR ("lower bandwidth must be less than M", GSL_EDOM);
    }
  else if (ub >= N)
    {
      GSL_ERROR ("upper bandwidth must be less than N", GSL_EDOM);
    }
  else if (AB->size2 != 2*lb + ub + 1)
    {
      GSL_ERROR ("matrix size inconsistent with bandwidths", GSL_EBADLEN);
    }
  else if (piv->size != minMN)
    {
      GSL_ERROR ("pivot vector must have length MIN(M,N)", GSL_EBADLEN);
    }
  else
    {
      int status;

      status = LU_band_decomp_L2 (M, lb, ub, AB, piv);

      return status;
    }
}

int
gsl_linalg_LU_band_solve (const size_t lb, const size_t ub, const gsl_matrix * LUB,
                          const gsl_vector_uint * piv, const gsl_vector * b, gsl_vector * x)
{
  const size_t N = LUB->size1;

  if (N != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (N != b->size)
    {
      GSL_ERROR ("matrix size must match rhs size", GSL_EBADLEN);
    }
  else if (lb >= N)
    {
      GSL_ERROR ("lower bandwidth must be less than N", GSL_EDOM);
    }
  else if (ub >= N)
    {
      GSL_ERROR ("upper bandwidth must be less than N", GSL_EDOM);
    }
  else if (LUB->size2 != 2*lb + ub + 1)
    {
      GSL_ERROR ("matrix size inconsistent with bandwidths", GSL_EBADLEN);
    }
  else if (piv->size != N)
    {
      GSL_ERROR ("pivot vector must have length N", GSL_EBADLEN);
    }
  else
    {
      int status;

      gsl_vector_memcpy(x, b);

      status = gsl_linalg_LU_band_svx(lb, ub, LUB, piv, x);

      return status;
    }
}

int
gsl_linalg_LU_band_svx (const size_t lb, const size_t ub, const gsl_matrix * LUB,
                        const gsl_vector_uint * piv, gsl_vector * x)
{
  const size_t N = LUB->size1;

  if (N != x->size)
    {
      GSL_ERROR ("matrix size must match solution/rhs size", GSL_EBADLEN);
    }
  else if (lb >= N)
    {
      GSL_ERROR ("lower bandwidth must be less than N", GSL_EDOM);
    }
  else if (ub >= N)
    {
      GSL_ERROR ("upper bandwidth must be less than N", GSL_EDOM);
    }
  else if (LUB->size2 != 2*lb + ub + 1)
    {
      GSL_ERROR ("matrix size inconsistent with bandwidths", GSL_EBADLEN);
    }
  else if (piv->size != N)
    {
      GSL_ERROR ("pivot vector must have length N", GSL_EBADLEN);
    }
  else
    {
      if (lb > 0)
        {
          size_t j;

          for (j = 0; j < N - 1; ++j)
            {
              size_t pj = gsl_vector_uint_get(piv, j);
              double * xj = gsl_vector_ptr(x, j);
              size_t lm = GSL_MIN(lb, N - j - 1);
              gsl_vector_view xv = gsl_vector_subvector(x, j + 1, lm);
              gsl_vector_const_view yv = gsl_matrix_const_subrow(LUB, j, lb + ub + 1, lm);

              if (j != pj)
                {
                  double xl = gsl_vector_get(x, pj);
                  gsl_vector_set(x, pj, *xj);
                  *xj = xl;
                }

              gsl_blas_daxpy(-(*xj), &yv.vector, &xv.vector);
            }
        }

      /* solve U x = b */
      cblas_dtbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                  (int) N, (int) (lb + ub), LUB->data, LUB->tda,
                  x->data, x->stride);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_LU_band_unpack (const size_t M, const size_t lb, const size_t ub, const gsl_matrix * LUB,
                           const gsl_vector_uint * piv, gsl_matrix * L, gsl_matrix * U)
{
  const size_t N = LUB->size1;
  const size_t minMN = GSL_MIN(M, N);

  if (ub >= N)
    {
      GSL_ERROR ("upper bandwidth must be < N", GSL_EDOM);
    }
  else if (lb >= M)
    {
      GSL_ERROR ("lower bandwidth must be < M", GSL_EDOM);
    }
  else if (LUB->size2 != 2*lb + ub + 1)
    {
      GSL_ERROR ("matrix size inconsistent with bandwidths", GSL_EBADLEN);
    }
  else if (piv->size != minMN)
    {
      GSL_ERROR ("pivot vector must have length MIN(M,N)", GSL_EBADLEN);
    }
  else if (L->size1 != M || L->size2 != minMN)
    {
      GSL_ERROR ("L matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (U->size1 != minMN || U->size2 != N)
    {
      GSL_ERROR ("U matrix has wrong dimensions", GSL_EBADLEN);
    }
  else
    {
      const size_t ub_U = lb + ub;
      size_t j;

      gsl_matrix_set_identity(L);
      gsl_matrix_set_zero(U);

      /* compute L */
      if (lb > 0)
        {
          const size_t jstart = (M > N) ? minMN : minMN - 1;
          size_t j;

          for (j = jstart; j > 0 && j--; )
            {
              size_t pj = gsl_vector_uint_get(piv, j);
              size_t lm = GSL_MIN(lb, M - j - 1);
              gsl_vector_const_view xv = gsl_matrix_const_subrow(LUB, j, lb + ub + 1, lm);
              gsl_vector_const_view yv = gsl_matrix_const_subrow(L, j, 0, minMN);
              gsl_matrix_view m = gsl_matrix_submatrix(L, j + 1, 0, lm, minMN);

              gsl_blas_dger(1.0, &xv.vector, &yv.vector, &m.matrix);

              if (j != pj)
                {
                  gsl_vector_view Lj = gsl_matrix_row(L, j);
                  gsl_vector_view Lpj = gsl_matrix_row(L, pj);
                  gsl_blas_dswap(&Lj.vector, &Lpj.vector);
                }
            }
        }

      /* fill in U */
      for (j = 0; j <= GSL_MIN(ub_U, N - 1); ++j)
        {
          gsl_vector_const_view src = gsl_matrix_const_subcolumn(LUB, ub_U - j, j, GSL_MIN(M, N - j));
          gsl_vector_view dest = gsl_matrix_superdiagonal(U, j);
          gsl_vector_memcpy(&dest.vector, &src.vector);
        }

      return GSL_SUCCESS;
    }
}

/*
LU_band_decomp_L2
  LU decomposition with partial pivoting using Level 2 BLAS

Inputs: M    - number of rows in matrix
        lb   - lower bandwidth
        ub   - upper bandwidth
        AB   - on input, matrix to be factored; on output, L and U factors
               N-by-(2*p + q + 1)
        ipiv - (output) array containing row swaps

Notes:
1) Based on LAPACK DGBTF2
*/

static int
LU_band_decomp_L2 (const size_t M, const size_t lb, const size_t ub,
                   gsl_matrix * AB, gsl_vector_uint * ipiv)
{
  const size_t N = AB->size1;
  const size_t minMN = GSL_MIN(M, N);

  if (ipiv->size != minMN)
    {
      GSL_ERROR ("ipiv length must equal MIN(M,N)", GSL_EBADLEN);
    }
  else if (lb >= M)
    {
      GSL_ERROR ("lower bandwidth must be less than M", GSL_EDOM);
    }
  else if (ub >= N)
    {
      GSL_ERROR ("upper bandwidth must be less than N", GSL_EDOM);
    }
  else if (AB->size2 != 2*lb + ub + 1)
    {
      GSL_ERROR ("matrix size inconsistent with bandwidths", GSL_EBADLEN);
    }
  else
    {
      const size_t ub_U = lb + ub; /* upper bandwidth of U factor */
      const size_t ldab = AB->size2;
      size_t ju = 0;
      size_t j;

      if (lb > 0)
        {
          /* initialize fill-in elements to zero */
          gsl_matrix_view m = gsl_matrix_submatrix(AB, 0, 0, N, lb);
          gsl_matrix_set_zero(&m.matrix);
        }
      
      for (j = 0; j < minMN; ++j)
        {
          size_t lbj = GSL_MIN(lb, M - j - 1); /* subdiagonal elements in column j */
          gsl_vector_view x = gsl_matrix_subrow(AB, j, ub_U, lbj + 1);
          gsl_vector_view y;
          CBLAS_INDEX_t j_pivot = gsl_blas_idamax(&x.vector);
          double * ptr;

          gsl_vector_uint_set(ipiv, j, j + j_pivot);

          ptr = gsl_matrix_ptr(AB, j, ub_U + j_pivot);
          if (*ptr != 0.0)
            ju = GSL_MAX(ju, GSL_MIN(j + ub + j_pivot, N - 1));

          if (j_pivot != 0)
            {
              double *ptr2;

              /* swap columns */

              x = gsl_vector_view_array_with_stride(ptr, ldab - 1, ju - j + 1);

              ptr2 = gsl_matrix_ptr(AB, j, ub_U);
              y = gsl_vector_view_array_with_stride(ptr2, ldab - 1, ju - j + 1);

              gsl_blas_dswap(&x.vector, &y.vector);
            }

          if (lbj > 0)
            {
              double tmp = gsl_matrix_get(AB, j, ub_U);

              x = gsl_matrix_subrow(AB, j, ub_U + 1, lbj);
              gsl_blas_dscal(1.0 / tmp, &x.vector);

              if (ju > j)
                {
                  gsl_matrix_view m = gsl_matrix_submatrix(AB, j + 1, ub_U, ju - j, lbj);

                  ptr = gsl_matrix_ptr(AB, j + 1, ub_U - 1);
                  y = gsl_vector_view_array_with_stride(ptr, ldab - 1, ju - j);

                  m.matrix.tda = ldab - 1; /* unskew matrix */
                  gsl_blas_dger(-1.0, &y.vector, &x.vector, &m.matrix);
                }
            }
        }

      return GSL_SUCCESS;
    }
}
