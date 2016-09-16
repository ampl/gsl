/* linalg/invtri.c
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

static int triangular_inverse(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix * T);

int
gsl_linalg_tri_upper_invert(gsl_matrix * T)
{
  int status = triangular_inverse(CblasUpper, CblasNonUnit, T);
  return status;
}

int
gsl_linalg_tri_lower_invert(gsl_matrix * T)
{
  int status = triangular_inverse(CblasLower, CblasNonUnit, T);
  return status;
}

int
gsl_linalg_tri_upper_unit_invert(gsl_matrix * T)
{
  int status = triangular_inverse(CblasUpper, CblasUnit, T);
  return status;
}

int
gsl_linalg_tri_lower_unit_invert(gsl_matrix * T)
{
  int status = triangular_inverse(CblasLower, CblasUnit, T);
  return status;
}

/*
triangular_inverse()
  Invert a triangular matrix T

Inputs: Uplo - CblasUpper or CblasLower
        Diag - unit triangular?
        T    - on output the upper (or lower) part of T
               is replaced by its inverse

Return: success/error
*/

static int
triangular_inverse(CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag, gsl_matrix * T)
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
