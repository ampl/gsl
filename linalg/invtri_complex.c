/* linalg/invtri_complex.c
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
 * This module contains code to invert complex triangular matrices
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "recurse.h"

static int complex_tri_invert_L2(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix_complex * T);
static int complex_tri_invert_L3(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix_complex * T);
static int triangular_singular(const gsl_matrix_complex * T);

/*
gsl_linalg_complex_tri_invert()
  Invert a complex triangular matrix T using Level 3 BLAS

Inputs: Uplo - CblasUpper or CblasLower
        Diag - unit triangular?
        T    - on output the upper (or lower) part of T
               is replaced by its inverse

Return: success/error
*/

int
gsl_linalg_complex_tri_invert(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix_complex * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      int status;

      status = triangular_singular(T);
      if (status)
        return status;

      return complex_tri_invert_L3(Uplo, Diag, T);
    }
}

/*
complex_tri_invert_L2()
  Invert a complex triangular matrix T

Inputs: Uplo - CblasUpper or CblasLower
        Diag - unit triangular?
        T    - on output the upper (or lower) part of T
               is replaced by its inverse

Return: success/error

Notes:
1) Based on LAPACK routine ZTRTI2 using Level 2 BLAS
*/

static int
complex_tri_invert_L2(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix_complex * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      size_t i;

      if (Uplo == CblasUpper)
        {
          for (i = 0; i < N; ++i)
            {
              gsl_complex * Tii = gsl_matrix_complex_ptr(T, i, i);
              gsl_complex aii;

              if (Diag == CblasNonUnit)
                {
                  *Tii = gsl_complex_inverse(*Tii);
                  GSL_REAL(aii) = -GSL_REAL(*Tii);
                  GSL_IMAG(aii) = -GSL_IMAG(*Tii);
                }
              else
                aii = GSL_COMPLEX_NEGONE;

              if (i > 0)
                {
                  gsl_matrix_complex_view m = gsl_matrix_complex_submatrix(T, 0, 0, i, i);
                  gsl_vector_complex_view v = gsl_matrix_complex_subcolumn(T, i, 0, i);

                  gsl_blas_ztrmv(CblasUpper, CblasNoTrans, Diag, &m.matrix, &v.vector);
                  gsl_blas_zscal(aii, &v.vector);
                }
            }
        }
      else
        {
          for (i = 0; i < N; ++i)
            {
              size_t j = N - i - 1;
              gsl_complex * Tjj = gsl_matrix_complex_ptr(T, j, j);
              gsl_complex ajj;

              if (Diag == CblasNonUnit)
                {
                  *Tjj = gsl_complex_inverse(*Tjj);
                  GSL_REAL(ajj) = -GSL_REAL(*Tjj);
                  GSL_IMAG(ajj) = -GSL_IMAG(*Tjj);
                }
              else
                ajj = GSL_COMPLEX_NEGONE;

              if (j < N - 1)
                {
                  gsl_matrix_complex_view m = gsl_matrix_complex_submatrix(T, j + 1, j + 1, N - j - 1, N - j - 1);
                  gsl_vector_complex_view v = gsl_matrix_complex_subcolumn(T, j, j + 1, N - j - 1);

                  gsl_blas_ztrmv(CblasLower, CblasNoTrans, Diag, &m.matrix, &v.vector);
                  gsl_blas_zscal(ajj, &v.vector);
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
complex_tri_invert_L3()
  Invert a complex triangular matrix T

Inputs: Uplo - CblasUpper or CblasLower
        Diag - unit triangular?
        T    - on output the upper (or lower) part of T
               is replaced by its inverse

Return: success/error

Notes:
1) Based on ReLAPACK using Level 3 BLAS
*/

static int
complex_tri_invert_L3(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix_complex * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (N <= CROSSOVER_INVTRI)
    {
      /* use Level 2 algorithm */
      return complex_tri_invert_L2(Uplo, Diag, T);
    }
  else
    {
      /*
       * partition matrix:
       *
       * T11 T12
       * T21 T22
       *
       * where T11 is N1-by-N1
       */
      int status;
      const size_t N1 = GSL_LINALG_SPLIT_COMPLEX(N);
      const size_t N2 = N - N1;
      gsl_matrix_complex_view T11 = gsl_matrix_complex_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_complex_view T12 = gsl_matrix_complex_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_complex_view T21 = gsl_matrix_complex_submatrix(T, N1, 0, N2, N1);
      gsl_matrix_complex_view T22 = gsl_matrix_complex_submatrix(T, N1, N1, N2, N2);

      /* recursion on T11 */
      status = complex_tri_invert_L3(Uplo, Diag, &T11.matrix);
      if (status)
        return status;

      if (Uplo == CblasLower)
        {
          /* T21 = - T21 * T11 */
          gsl_blas_ztrmm(CblasRight, Uplo, CblasNoTrans, Diag, GSL_COMPLEX_NEGONE, &T11.matrix, &T21.matrix);

          /* T21 = T22 * T21^{-1} */
          gsl_blas_ztrsm(CblasLeft, Uplo, CblasNoTrans, Diag, GSL_COMPLEX_ONE, &T22.matrix, &T21.matrix);
        }
      else
        {
          /* T12 = - T11 * T12 */
          gsl_blas_ztrmm(CblasLeft, Uplo, CblasNoTrans, Diag, GSL_COMPLEX_NEGONE, &T11.matrix, &T12.matrix);

          /* T12 = T12 * T22^{-1} */
          gsl_blas_ztrsm(CblasRight, Uplo, CblasNoTrans, Diag, GSL_COMPLEX_ONE, &T22.matrix, &T12.matrix);
        }

      /* recursion on T22 */
      status = complex_tri_invert_L3(Uplo, Diag, &T22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

static int
triangular_singular(const gsl_matrix_complex * T)
{
  size_t i;

  for (i = 0; i < T->size1; ++i)
    {
      gsl_complex z = gsl_matrix_complex_get(T, i, i);
      if (GSL_REAL(z) == 0.0 && GSL_IMAG(z) == 0.0)
        return GSL_ESING;
    }

  return GSL_SUCCESS;
}
