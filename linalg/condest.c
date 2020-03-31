/* linalg/condest.c
 * 
 * Copyright (C) 2016 Patrick Alken
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/*
 * This module contains routines for estimating the condition number
 * of matrices in the 1-norm. The algorithm is based on the paper,
 *
 * [1] N. J. Higham, "FORTRAN codes for estimating the one-norm of
 * a real or complex matrix, with applications to condition estimation",
 * ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
 */

static double condest_tri_norm1(CBLAS_UPLO_t Uplo, const gsl_matrix * A);
static int condest_tri_rcond(CBLAS_UPLO_t Uplo, const gsl_matrix * A,
                             double * rcond, gsl_vector * work);
static int condest_same_sign(const gsl_vector * x, const gsl_vector * y);
static int condest_invtriu(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params);
static int condest_invtril(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params);

int
gsl_linalg_tri_rcond(CBLAS_UPLO_t Uplo, const gsl_matrix * A, double * rcond, gsl_vector * work)
{
  return condest_tri_rcond(Uplo, A, rcond, work);
}

/*
gsl_linalg_invnorm1()
  Estimate the 1-norm of ||A^{-1}||, where A is a square
N-by-N matrix

Inputs: N        - size of matrix
        Ainvx    - pointer to function which calculates:
                   x := A^{-1} x or x := A^{-t} x
        params   - parameters to pass to Ainvx
        Ainvnorm - (output) estimate of ||A^{-1}||_1
        work     - workspace, length 3*N
*/

int
gsl_linalg_invnorm1(const size_t N,
                    int (* Ainvx)(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params),
                    void * params, double * Ainvnorm, gsl_vector * work)
{
  if (work->size != 3 * N)
    {
      GSL_ERROR ("work vector must have length 3*N", GSL_EBADLEN);
    }
  else
    {
      const size_t maxit = 5;
      gsl_vector_view x = gsl_vector_subvector(work, 0, N);
      gsl_vector_view v = gsl_vector_subvector(work, N, N);
      gsl_vector_view xi = gsl_vector_subvector(work, 2*N, N);
      double gamma, gamma_old, temp;
      size_t i, k;

      for (i = 0; i < N; ++i)
        gsl_vector_set(&x.vector, i, 1.0 / (double) N);

      /* compute v = A^{-1} x */
      gsl_vector_memcpy(&v.vector, &x.vector);
      (*Ainvx)(CblasNoTrans, &v.vector, params);

      /* gamma = ||v||_1 */
      gamma = gsl_blas_dasum(&v.vector);

      /* xi = sign(v) */
      for (i = 0; i < N; ++i)
        {
          double vi = gsl_vector_get(&v.vector, i);
          gsl_vector_set(&xi.vector, i, GSL_SIGN(vi));
        }

      /* x = A^{-t} xi */
      gsl_vector_memcpy(&x.vector, &xi.vector);
      (*Ainvx)(CblasTrans, &x.vector, params);

      for (k = 0; k < maxit; ++k)
        {
          size_t j = (size_t) gsl_blas_idamax(&x.vector);

          /* v := A^{-1} e_j */
          gsl_vector_set_zero(&v.vector);
          gsl_vector_set(&v.vector, j, 1.0);
          (*Ainvx)(CblasNoTrans, &v.vector, params);

          gamma_old = gamma;
          gamma = gsl_blas_dasum(&v.vector);

          /* check for repeated sign vector (algorithm has converged) */
          if (condest_same_sign(&v.vector, &xi.vector) || (gamma < gamma_old))
            break;

          /* xi = sign(v) */
          for (i = 0; i < N; ++i)
            {
              double vi = gsl_vector_get(&v.vector, i);
              gsl_vector_set(&xi.vector, i, GSL_SIGN(vi));
            }

          /* x = A^{-t} sign(v) */
          gsl_vector_memcpy(&x.vector, &xi.vector);
          (*Ainvx)(CblasTrans, &x.vector, params);
        }

      temp = 1.0; /* (-1)^i */
      for (i = 0; i < N; ++i)
        {
          double term = 1.0 + (double) i / (N - 1.0);
          gsl_vector_set(&x.vector, i, temp * term);
          temp = -temp;
        }

      /* x := A^{-1} x */
      (*Ainvx)(CblasNoTrans, &x.vector, params);

      temp = 2.0 * gsl_blas_dasum(&x.vector) / (3.0 * N);
      if (temp > gamma)
        {
          gsl_vector_memcpy(&v.vector, &x.vector);
          gamma = temp;
        }

      *Ainvnorm = gamma;

      return GSL_SUCCESS;
    }
}

static int
condest_tri_rcond(CBLAS_UPLO_t Uplo, const gsl_matrix * A, double * rcond, gsl_vector * work)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (work->size != 3 * N)
    {
      GSL_ERROR ("work vector must have length 3*N", GSL_EBADLEN);
    }
  else
    {
      int status;
      double Anorm = condest_tri_norm1(Uplo, A); /* ||A||_1 */
      double Ainvnorm;                           /* ||A^{-1}||_1 */

      *rcond = 0.0;

      /* don't continue if matrix is singular */
      if (Anorm == 0.0)
        return GSL_SUCCESS;

      /* estimate ||A^{-1}||_1 */
      if (Uplo == CblasUpper)
        status = gsl_linalg_invnorm1(N, condest_invtriu, (void *) A, &Ainvnorm, work);
      else
        status = gsl_linalg_invnorm1(N, condest_invtril, (void *) A, &Ainvnorm, work);

      if (status)
        return status;

      if (Ainvnorm != 0.0)
        *rcond = (1.0 / Anorm) / Ainvnorm;

      return GSL_SUCCESS;
    }
}

/* calculate 1 norm of triangular matrix */
static double
condest_tri_norm1(CBLAS_UPLO_t Uplo, const gsl_matrix * A)
{
  const size_t N = A->size2;
  double max = 0.0;
  size_t i, j;

  if (Uplo == CblasUpper)
    {
      for (j = 0; j < N; ++j)
        {
          double sum = 0.0;
          for (i = 0; i <= j; ++i)
            {
              double Aij = gsl_matrix_get(A, i, j);
              sum += fabs(Aij);
            }

          max = GSL_MAX(max, sum);
        }
    }
  else
    {
      for (j = 0; j < N; ++j)
        {
          double sum = 0.0;
          for (i = j; i < N; ++i)
            {
              double Aij = gsl_matrix_get(A, i, j);
              sum += fabs(Aij);
            }

          max = GSL_MAX(max, sum);
        }
    }

  return max;
}

/* return 1 if sign(x) = sign(y), 0 otherwise */
static int
condest_same_sign(const gsl_vector * x, const gsl_vector * y)
{
  const size_t n = x->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(y, i);
      if (GSL_SIGN(xi) != GSL_SIGN(yi))
        return 0;
    }

  return 1;
}

/* x := A^{-1} x, A upper triangular */
static int
condest_invtriu(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params)
{
  gsl_matrix * A = (gsl_matrix *) params;
  return gsl_blas_dtrsv(CblasUpper, TransA, CblasNonUnit, A, x);
}

/* x := A^{-1} x, A lower triangular */
static int
condest_invtril(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params)
{
  gsl_matrix * A = (gsl_matrix *) params;
  return gsl_blas_dtrsv(CblasLower, TransA, CblasNonUnit, A, x);
}

#ifndef GSL_DISABLE_DEPRECATED

int
gsl_linalg_tri_upper_rcond(const gsl_matrix * A, double * rcond, gsl_vector * work)
{
  int status = condest_tri_rcond(CblasUpper, A, rcond, work);
  return status;
}

int
gsl_linalg_tri_lower_rcond(const gsl_matrix * A, double * rcond, gsl_vector * work)
{
  int status = condest_tri_rcond(CblasLower, A, rcond, work);
  return status;
}

#endif
