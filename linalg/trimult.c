/* linalg/trimult.c
 *
 * Copyright (C) 2019 Patrick Alken
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * This module contains code to compute L^T L where L is a lower triangular matrix
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "recurse.h"

static int triangular_multsymm_L2(CBLAS_UPLO_t Uplo, gsl_matrix * T);
static int triangular_multsymm_L3(CBLAS_UPLO_t Uplo, gsl_matrix * T);
static int triangular_mult_L2(CBLAS_UPLO_t Uplo, gsl_matrix * A);
static int triangular_mult_L3(CBLAS_UPLO_t Uplo, gsl_matrix * A);

int
gsl_linalg_tri_LTL(gsl_matrix * L)
{
  return triangular_multsymm_L3(CblasLower, L);
}

int
gsl_linalg_tri_UL(gsl_matrix * LU)
{
  return triangular_mult_L3(CblasUpper, LU);
}

/*
triangular_multsymm_L2()
  Compute L^T L or U U^T

Inputs: Uplo - CblasUpper or CblasLower
        T    - on output the upper (or lower) part of T
               is replaced by L^T L or U U^T

Return: success/error

Notes:
1) Based on LAPACK routine DLAUU2 using Level 2 BLAS
*/

static int
triangular_multsymm_L2(CBLAS_UPLO_t Uplo, gsl_matrix * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      gsl_vector_view v1, v2;
      size_t i;

      if (Uplo == CblasUpper)
        {
        }
      else
        {
          for (i = 0; i < N; ++i)
            {
              double Tii = gsl_matrix_get(T, i, i);

              if (i < N - 1)
                {
                  double tmp;

                  v1 = gsl_matrix_subcolumn(T, i, i, N - i);
                  gsl_blas_ddot(&v1.vector, &v1.vector, &tmp);
                  gsl_matrix_set(T, i, i, tmp);

                  if (i > 0)
                    {
                      gsl_matrix_view m = gsl_matrix_submatrix(T, i + 1, 0, N - i - 1, i);

                      v1 = gsl_matrix_subcolumn(T, i, i + 1, N - i - 1);
                      v2 = gsl_matrix_subrow(T, i, 0, i);

                      gsl_blas_dgemv(CblasTrans, 1.0, &m.matrix, &v1.vector, Tii, &v2.vector);
                    }
                }
              else
                {
                  v1 = gsl_matrix_row(T, N - 1);
                  gsl_blas_dscal(Tii, &v1.vector);
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
triangular_multsymm_L3()
  Compute L^T L or U U^T

Inputs: Uplo - CblasUpper or CblasLower
        T    - on output the upper (or lower) part of T
               is replaced by L^T L or U U^T

Return: success/error

Notes:
1) Based on ReLAPACK routine DLAUUM using Level 3 BLAS
*/

static int
triangular_multsymm_L3(CBLAS_UPLO_t Uplo, gsl_matrix * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (N <= CROSSOVER_TRIMULT)
    {
      return triangular_multsymm_L2(Uplo, T);
    }
  else
    {
      /* partition matrix:
       *
       * T11 T12
       * T21 T22
       *
       * where T11 is N1-by-N1
       */
      int status;
      const size_t N1 = GSL_LINALG_SPLIT(N);
      const size_t N2 = N - N1;
      gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_view T21 = gsl_matrix_submatrix(T, N1, 0, N2, N1);
      gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

      /* recursion on T11 */
      status = triangular_multsymm_L3(Uplo, &T11.matrix);
      if (status)
        return status;

      if (Uplo == CblasLower)
        {
          /* T11 += T21^T T21 */
          gsl_blas_dsyrk(Uplo, CblasTrans, 1.0, &T21.matrix, 1.0, &T11.matrix);

          /* T21 = T22^T * T21 */
          gsl_blas_dtrmm(CblasLeft, Uplo, CblasTrans, CblasNonUnit, 1.0, &T22.matrix, &T21.matrix);
        }
      else
        {
          /* T11 += T12 T12^T */
          gsl_blas_dsyrk(Uplo, CblasNoTrans, 1.0, &T12.matrix, 1.0, &T11.matrix);

          /* T12 = T12 * T22^T */
          gsl_blas_dtrmm(CblasRight, Uplo, CblasTrans, CblasNonUnit, 1.0, &T22.matrix, &T12.matrix);
        }

      /* recursion on T22 */
      status = triangular_multsymm_L3(Uplo, &T22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

/*
triangular_mult_L2()
  Compute U L or L U

Inputs: Uplo - CblasUpper or CblasLower (for the first triangular factor)
        A    - on input, matrix in LU format;
               on output, U L or L U

Return: success/error
*/

static int
triangular_mult_L2(CBLAS_UPLO_t Uplo, gsl_matrix * A)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      size_t i;

      /* quick return */
      if (N == 1)
        return GSL_SUCCESS;

      if (Uplo == CblasUpper)
        {
          /* compute U * L and store in A */

          for (i = 0; i < N; ++i)
            {
              double * Aii = gsl_matrix_ptr(A, i, i);
              double Uii = *Aii;

              if (i < N - 1)
                {
                  gsl_vector_view lb = gsl_matrix_subcolumn(A, i, i + 1, N - i - 1);
                  gsl_vector_view ur = gsl_matrix_subrow(A, i, i + 1, N - i - 1);
                  double tmp;

                  gsl_blas_ddot(&lb.vector, &ur.vector, &tmp);
                  *Aii += tmp;

                  if (i > 0)
                    {
                      gsl_matrix_view U_TR = gsl_matrix_submatrix(A, 0, i + 1, i, N - i - 1);
                      gsl_matrix_view L_BL = gsl_matrix_submatrix(A, i + 1, 0, N - i - 1, i);
                      gsl_vector_view ut = gsl_matrix_subcolumn(A, i, 0, i);
                      gsl_vector_view ll = gsl_matrix_subrow(A, i, 0, i);
                      
                      gsl_blas_dgemv(CblasTrans, 1.0, &L_BL.matrix, &ur.vector, Uii, &ll.vector);
                      gsl_blas_dgemv(CblasNoTrans, 1.0, &U_TR.matrix, &lb.vector, 1.0, &ut.vector);
                    }
                }
              else
                {
                  gsl_vector_view v = gsl_matrix_subrow(A, N - 1, 0, N - 1);
                  gsl_blas_dscal(Uii, &v.vector);
                }
            }
        }
      else
        {
        }

      return GSL_SUCCESS;
    }
}

/*
triangular_mult_L3()
  Compute U L or L U

Inputs: Uplo - CblasUpper or CblasLower (for the first triangular factor)
        A    - on input, matrix in LU format;
               on output, U L or L U

Return: success/error
*/

static int
triangular_mult_L3(CBLAS_UPLO_t Uplo, gsl_matrix * A)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (N <= CROSSOVER_TRIMULT)
    {
      return triangular_mult_L2(Uplo, A);
    }
  else
    {
      /* partition matrix:
       *
       * A11 A12
       * A21 A22
       *
       * where A11 is N1-by-N1
       */
      int status;
      const size_t N1 = GSL_LINALG_SPLIT(N);
      const size_t N2 = N - N1;
      gsl_matrix_view A11 = gsl_matrix_submatrix(A, 0, 0, N1, N1);
      gsl_matrix_view A12 = gsl_matrix_submatrix(A, 0, N1, N1, N2);
      gsl_matrix_view A21 = gsl_matrix_submatrix(A, N1, 0, N2, N1);
      gsl_matrix_view A22 = gsl_matrix_submatrix(A, N1, N1, N2, N2);

      /* recursion on A11 */
      status = triangular_mult_L3(Uplo, &A11.matrix);
      if (status)
        return status;

      if (Uplo == CblasLower)
        {
        }
      else
        {
          /* form U * L */

          /* A11 += A12 A21 */
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &A12.matrix, &A21.matrix, 1.0, &A11.matrix);

          /* A12 = A12 * L22 */
          gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasUnit, 1.0, &A22.matrix, &A12.matrix);

          /* A21 = U22 * A21 */
          gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &A22.matrix, &A21.matrix);
        }

      /* recursion on A22 */
      status = triangular_mult_L3(Uplo, &A22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}
