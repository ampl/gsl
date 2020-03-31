/* eigen/genherm.c
 * 
 * Copyright (C) 2007, 2019 Patrick Alken
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

#include <stdlib.h>

#include <config.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "recurse.h"

/*
 * This module computes the eigenvalues of a complex generalized
 * hermitian-definite eigensystem A x = \lambda B x, where A and
 * B are hermitian, and B is positive-definite.
 */

static int genherm_standardize_L2(gsl_matrix_complex *A, const gsl_matrix_complex *B);
static int genherm_standardize_L3(gsl_matrix_complex *A, const gsl_matrix_complex *B);

/*
gsl_eigen_genherm_alloc()

Allocate a workspace for solving the generalized hermitian-definite
eigenvalue problem. The size of this workspace is O(3n).

Inputs: n - size of matrices

Return: pointer to workspace
*/

gsl_eigen_genherm_workspace *
gsl_eigen_genherm_alloc(const size_t n)
{
  gsl_eigen_genherm_workspace *w;

  if (n == 0)
    {
      GSL_ERROR_NULL ("matrix dimension must be positive integer",
                      GSL_EINVAL);
    }

  w = (gsl_eigen_genherm_workspace *) calloc (1, sizeof (gsl_eigen_genherm_workspace));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->size = n;

  w->herm_workspace_p = gsl_eigen_herm_alloc(n);
  if (!w->herm_workspace_p)
    {
      gsl_eigen_genherm_free(w);
      GSL_ERROR_NULL("failed to allocate space for herm workspace", GSL_ENOMEM);
    }

  return (w);
} /* gsl_eigen_genherm_alloc() */

/*
gsl_eigen_genherm_free()
  Free workspace w
*/

void
gsl_eigen_genherm_free (gsl_eigen_genherm_workspace * w)
{
  RETURN_IF_NULL (w);

  if (w->herm_workspace_p)
    gsl_eigen_herm_free(w->herm_workspace_p);

  free(w);
} /* gsl_eigen_genherm_free() */

/*
gsl_eigen_genherm()

Solve the generalized hermitian-definite eigenvalue problem

A x = \lambda B x

for the eigenvalues \lambda.

Inputs: A    - complex hermitian matrix
        B    - complex hermitian and positive definite matrix
        eval - where to store eigenvalues
        w    - workspace

Return: success or error
*/

int
gsl_eigen_genherm (gsl_matrix_complex * A, gsl_matrix_complex * B,
                   gsl_vector * eval, gsl_eigen_genherm_workspace * w)
{
  const size_t N = A->size1;

  /* check matrix and vector sizes */

  if (N != A->size2)
    {
      GSL_ERROR ("matrix must be square to compute eigenvalues", GSL_ENOTSQR);
    }
  else if ((N != B->size1) || (N != B->size2))
    {
      GSL_ERROR ("B matrix dimensions must match A", GSL_EBADLEN);
    }
  else if (eval->size != N)
    {
      GSL_ERROR ("eigenvalue vector must match matrix size", GSL_EBADLEN);
    }
  else if (w->size != N)
    {
      GSL_ERROR ("matrix size does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int s;

      /* compute Cholesky factorization of B */
      s = gsl_linalg_complex_cholesky_decomp(B);
      if (s != GSL_SUCCESS)
        return s; /* B is not positive definite */

      /* transform to standard hermitian eigenvalue problem */
      gsl_eigen_genherm_standardize(A, B);

      s = gsl_eigen_herm(A, eval, w->herm_workspace_p);

      return s;
    }
} /* gsl_eigen_genherm() */

/*
gsl_eigen_genherm_standardize()
  Reduce the generalized hermitian-definite eigenproblem to
the standard hermitian eigenproblem by computing

C = L^{-1} A L^{-H}

where L L^H is the Cholesky decomposition of B

Inputs: A - (input/output) complex hermitian matrix
        B - complex hermitian, positive definite matrix in Cholesky form

Return: success

Notes: A is overwritten by L^{-1} A L^{-H}
*/

int
gsl_eigen_genherm_standardize(gsl_matrix_complex *A, const gsl_matrix_complex *B)
{
  return genherm_standardize_L3(A, B);
}

/*
genherm_standardize_L2()
  Reduce the generalized hermitian-definite eigenproblem to
the standard hermitian eigenproblem by computing

C = L^{-1} A L^{-H}

where L L^H is the Cholesky decomposition of B

Inputs: A - (input/output) complex hermitian matrix
        B - complex hermitian, positive definite matrix in Cholesky form

Return: success

Notes:
1) A is overwritten by L^{-1} A L^{-H}

2) Based on LAPACK ZHEGS2 using Level 2 BLAS
*/

static int
genherm_standardize_L2(gsl_matrix_complex *A, const gsl_matrix_complex *B)
{
  const size_t N = A->size1;
  size_t i;
  double a, b;
  gsl_complex y, z;

  GSL_SET_IMAG(&z, 0.0);

  for (i = 0; i < N; ++i)
    {
      /* update lower triangle of A(i:n, i:n) */

      y = gsl_matrix_complex_get(A, i, i);
      a = GSL_REAL(y);
      y = gsl_matrix_complex_get(B, i, i);
      b = GSL_REAL(y);
      a /= b * b;
      GSL_SET_REAL(&z, a);
      gsl_matrix_complex_set(A, i, i, z);

      if (i < N - 1)
        {
          gsl_vector_complex_view ai =
            gsl_matrix_complex_subcolumn(A, i, i + 1, N - i - 1);
          gsl_matrix_complex_view ma =
            gsl_matrix_complex_submatrix(A, i + 1, i + 1, N - i - 1, N - i - 1);
          gsl_vector_complex_const_view bi =
            gsl_matrix_complex_const_subcolumn(B, i, i + 1, N - i - 1);
          gsl_matrix_complex_const_view mb =
            gsl_matrix_complex_const_submatrix(B, i + 1, i + 1, N - i - 1, N - i - 1);

          gsl_blas_zdscal(1.0 / b, &ai.vector);

          GSL_SET_REAL(&z, -0.5 * a);
          gsl_blas_zaxpy(z, &bi.vector, &ai.vector);

          gsl_blas_zher2(CblasLower,
                         GSL_COMPLEX_NEGONE,
                         &ai.vector,
                         &bi.vector,
                         &ma.matrix);

          gsl_blas_zaxpy(z, &bi.vector, &ai.vector);

          gsl_blas_ztrsv(CblasLower,
                         CblasNoTrans,
                         CblasNonUnit,
                         &mb.matrix,
                         &ai.vector);
        }
    }

  return GSL_SUCCESS;
}

/*
genherm_standardize_L3()
  Reduce the generalized hermitian-definite eigenproblem to
the standard hermitian eigenproblem by computing

C = L^{-1} A L^{-H}

where L L^H is the Cholesky decomposition of B

Inputs: A - (input/output) complex hermitian matrix
        B - complex hermitian, positive definite matrix in Cholesky form

Return: success

Notes:
1) A is overwritten by L^{-1} A L^{-H}

2) Based on ReLAPACK using Level 3 BLAS
*/

static int
genherm_standardize_L3(gsl_matrix_complex *A, const gsl_matrix_complex *B)
{
  const size_t N = A->size1;

  if (N <= CROSSOVER_GENHERM)
    {
      /* use Level 2 algorithm */
      return genherm_standardize_L2(A, B);
    }
  else
    {
      /*
       * partition matrices:
       *
       * A11 A12  and  B11 B12
       * A21 A22       B21 B22
       *
       * where A11 and B11 are N1-by-N1
       */
      int status;
      const size_t N1 = GSL_EIGEN_SPLIT_COMPLEX(N);
      const size_t N2 = N - N1;

      gsl_matrix_complex_view A11 = gsl_matrix_complex_submatrix(A, 0, 0, N1, N1);
      gsl_matrix_complex_view A21 = gsl_matrix_complex_submatrix(A, N1, 0, N2, N1);
      gsl_matrix_complex_view A22 = gsl_matrix_complex_submatrix(A, N1, N1, N2, N2);

      gsl_matrix_complex_const_view B11 = gsl_matrix_complex_const_submatrix(B, 0, 0, N1, N1);
      gsl_matrix_complex_const_view B21 = gsl_matrix_complex_const_submatrix(B, N1, 0, N2, N1);
      gsl_matrix_complex_const_view B22 = gsl_matrix_complex_const_submatrix(B, N1, N1, N2, N2);

      const gsl_complex MHALF = gsl_complex_rect(-0.5, 0.0);

      /* recursion on (A11, B11) */
      status = genherm_standardize_L3(&A11.matrix, &B11.matrix);
      if (status)
        return status;

      /* A21 = A21 * B11^{-1} */
      gsl_blas_ztrsm(CblasRight, CblasLower, CblasConjTrans, CblasNonUnit, GSL_COMPLEX_ONE, &B11.matrix, &A21.matrix);

      /* A21 = A21 - 1/2 B21 A11 */
      gsl_blas_zhemm(CblasRight, CblasLower, MHALF, &A11.matrix, &B21.matrix, GSL_COMPLEX_ONE, &A21.matrix);

      /* A22 = A22 - A21 * B21' - B21 * A21' */
      gsl_blas_zher2k(CblasLower, CblasNoTrans, GSL_COMPLEX_NEGONE, &A21.matrix, &B21.matrix, 1.0, &A22.matrix);

      /* A21 = A21 - 1/2 B21 A11 */
      gsl_blas_zhemm(CblasRight, CblasLower, MHALF, &A11.matrix, &B21.matrix, GSL_COMPLEX_ONE, &A21.matrix);

      /* A21 = B22 * A21^{-1} */
      gsl_blas_ztrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, GSL_COMPLEX_ONE, &B22.matrix, &A21.matrix);

      /* recursion on (A22, B22) */
      status = genherm_standardize_L3(&A22.matrix, &B22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}
