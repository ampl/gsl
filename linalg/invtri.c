/* linalg/invtri.c
 *
 * Copyright (C) 2016, 2017, 2018, 2019 Patrick Alken
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
 * This module contains code to invert triangular matrices
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "recurse.h"

static int triangular_inverse_L2(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix * T);
static int triangular_inverse_L3(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix * T);
static int triangular_singular(const gsl_matrix * T);

int
gsl_linalg_tri_invert(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix * T)
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

      return triangular_inverse_L3(Uplo, Diag, T);
    }
}

/*
triangular_inverse_L2()
  Invert a triangular matrix T

Inputs: Uplo - CblasUpper or CblasLower
        Diag - unit triangular?
        T    - on output the upper (or lower) part of T
               is replaced by its inverse

Return: success/error

Notes:
1) Based on LAPACK routine DTRTI2 using Level 2 BLAS
*/

static int
triangular_inverse_L2(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      gsl_matrix_view m;
      gsl_vector_view v;
      size_t i;

      if (Uplo == CblasUpper)
        {
          for (i = 0; i < N; ++i)
            {
              double aii;

              if (Diag == CblasNonUnit)
                {
                  double *Tii = gsl_matrix_ptr(T, i, i);
                  *Tii = 1.0 / *Tii;
                  aii = -(*Tii);
                }
              else
                {
                  aii = -1.0;
                }

              if (i > 0)
                {
                  m = gsl_matrix_submatrix(T, 0, 0, i, i);
                  v = gsl_matrix_subcolumn(T, i, 0, i);

                  gsl_blas_dtrmv(CblasUpper, CblasNoTrans, Diag,
                                 &m.matrix, &v.vector);

                  gsl_blas_dscal(aii, &v.vector);
                }
            } /* for (i = 0; i < N; ++i) */
        }
      else
        {
          for (i = 0; i < N; ++i)
            {
              double ajj;
              size_t j = N - i - 1;

              if (Diag == CblasNonUnit)
                {
                  double *Tjj = gsl_matrix_ptr(T, j, j);
                  *Tjj = 1.0 / *Tjj;
                  ajj = -(*Tjj);
                }
              else
                {
                  ajj = -1.0;
                }

              if (j < N - 1)
                {
                  m = gsl_matrix_submatrix(T, j + 1, j + 1,
                                           N - j - 1, N - j - 1);
                  v = gsl_matrix_subcolumn(T, j, j + 1, N - j - 1);

                  gsl_blas_dtrmv(CblasLower, CblasNoTrans, Diag,
                                 &m.matrix, &v.vector);

                  gsl_blas_dscal(ajj, &v.vector);
                }
            } /* for (i = 0; i < N; ++i) */
        }

      return GSL_SUCCESS;
    }
}

/*
triangular_inverse_L3()
  Invert a triangular matrix T

Inputs: Uplo - CblasUpper or CblasLower
        Diag - unit triangular?
        T    - on output the upper (or lower) part of T
               is replaced by its inverse

Return: success/error

Notes:
1) Based on ReLAPACK routine DTRTRI using Level 3 BLAS
*/

static int
triangular_inverse_L3(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (N <= CROSSOVER_INVTRI)
    {
      /* use Level 2 BLAS code */
      return triangular_inverse_L2(Uplo, Diag, T);
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
      const size_t N1 = GSL_LINALG_SPLIT(N);
      const size_t N2 = N - N1;
      gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_view T21 = gsl_matrix_submatrix(T, N1, 0, N2, N1);
      gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

      /* recursion on T11 */
      status = triangular_inverse_L3(Uplo, Diag, &T11.matrix);
      if (status)
        return status;

      if (Uplo == CblasLower)
        {
          /* T21 = - T21 * T11 */
          gsl_blas_dtrmm(CblasRight, Uplo, CblasNoTrans, Diag, -1.0, &T11.matrix, &T21.matrix);

          /* T21 = T22 * T21^{-1} */
          gsl_blas_dtrsm(CblasLeft, Uplo, CblasNoTrans, Diag, 1.0, &T22.matrix, &T21.matrix);
        }
      else
        {
          /* T12 = - T11 * T12 */
          gsl_blas_dtrmm(CblasLeft, Uplo, CblasNoTrans, Diag, -1.0, &T11.matrix, &T12.matrix);

          /* T12 = T12 * T22^{-1} */
          gsl_blas_dtrsm(CblasRight, Uplo, CblasNoTrans, Diag, 1.0, &T22.matrix, &T12.matrix);
        }

      /* recursion on T22 */
      status = triangular_inverse_L3(Uplo, Diag, &T22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

static int
triangular_singular(const gsl_matrix * T)
{
  size_t i;

  for (i = 0; i < T->size1; ++i)
    {
      double Tii = gsl_matrix_get(T, i, i);
      if (Tii == 0.0)
        return GSL_ESING;
    }

  return GSL_SUCCESS;
}

#ifndef GSL_DISABLE_DEPRECATED

int
gsl_linalg_tri_upper_invert(gsl_matrix * T)
{
  int status = triangular_singular(T);
  if (status)
    return status;
  
  return triangular_inverse_L3(CblasUpper, CblasNonUnit, T);
}

int
gsl_linalg_tri_lower_invert(gsl_matrix * T)
{
  int status = triangular_singular(T);
  if (status)
    return status;

  return triangular_inverse_L3(CblasLower, CblasNonUnit, T);
}

int
gsl_linalg_tri_upper_unit_invert(gsl_matrix * T)
{
  int status = triangular_singular(T);
  if (status)
    return status;

  return triangular_inverse_L3(CblasUpper, CblasUnit, T);
}

int
gsl_linalg_tri_lower_unit_invert(gsl_matrix * T)
{
  int status = triangular_singular(T);
  if (status)
    return status;

  return triangular_inverse_L3(CblasLower, CblasUnit, T);
}

#endif
