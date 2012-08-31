/* linalg/choleskyc.c
 * 
 * Copyright (C) 2007 Patrick Alken
 * Copyright (C) 2010 Huan Wu (gsl_linalg_complex_cholesky_invert)
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

/*
 * This module contains routines related to the Cholesky decomposition
 * of a complex Hermitian positive definite matrix.
 */

static void cholesky_complex_conj_vector(gsl_vector_complex *v);

/*
gsl_linalg_complex_cholesky_decomp()
  Perform the Cholesky decomposition on a Hermitian positive definite
matrix. See Golub & Van Loan, "Matrix Computations" (3rd ed),
algorithm 4.2.2.

Inputs: A - (input/output) complex postive definite matrix

Return: success or error

The lower triangle of A is overwritten with the Cholesky decomposition
*/

int
gsl_linalg_complex_cholesky_decomp(gsl_matrix_complex *A)
{
  const size_t N = A->size1;
  
  if (N != A->size2)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      size_t i, j;
      gsl_complex z;
      double ajj;

      for (j = 0; j < N; ++j)
        {
          z = gsl_matrix_complex_get(A, j, j);
          ajj = GSL_REAL(z);

          if (j > 0)
            {
              gsl_vector_complex_const_view aj =
                gsl_matrix_complex_const_subrow(A, j, 0, j);

              gsl_blas_zdotc(&aj.vector, &aj.vector, &z);
              ajj -= GSL_REAL(z);
            }

          if (ajj <= 0.0)
            {
              GSL_ERROR("matrix is not positive definite", GSL_EDOM);
            }

          ajj = sqrt(ajj);
          GSL_SET_COMPLEX(&z, ajj, 0.0);
          gsl_matrix_complex_set(A, j, j, z);

          if (j < N - 1)
            {
              gsl_vector_complex_view av =
                gsl_matrix_complex_subcolumn(A, j, j + 1, N - j - 1);

              if (j > 0)
                {
                  gsl_vector_complex_view aj =
                    gsl_matrix_complex_subrow(A, j, 0, j);
                  gsl_matrix_complex_view am =
                    gsl_matrix_complex_submatrix(A, j + 1, 0, N - j - 1, j);

                  cholesky_complex_conj_vector(&aj.vector);

                  gsl_blas_zgemv(CblasNoTrans,
                                 GSL_COMPLEX_NEGONE,
                                 &am.matrix,
                                 &aj.vector,
                                 GSL_COMPLEX_ONE,
                                 &av.vector);

                  cholesky_complex_conj_vector(&aj.vector);
                }

              gsl_blas_zdscal(1.0 / ajj, &av.vector);
            }
        }

      /* Now store L^H in upper triangle */
      for (i = 1; i < N; ++i)
        {
          for (j = 0; j < i; ++j)
            {
              z = gsl_matrix_complex_get(A, i, j);
              gsl_matrix_complex_set(A, j, i, gsl_complex_conjugate(z));
            }
        }

      return GSL_SUCCESS;
    }
} /* gsl_linalg_complex_cholesky_decomp() */

/*
gsl_linalg_complex_cholesky_solve()
  Solve A x = b where A is in cholesky form
*/

int
gsl_linalg_complex_cholesky_solve (const gsl_matrix_complex * cholesky,
                                   const gsl_vector_complex * b,
                                   gsl_vector_complex * x)
{
  if (cholesky->size1 != cholesky->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (cholesky->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (cholesky->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_complex_memcpy (x, b);

      /* solve for y using forward-substitution, L y = b */

      gsl_blas_ztrsv (CblasLower, CblasNoTrans, CblasNonUnit, cholesky, x);

      /* perform back-substitution, L^H x = y */

      gsl_blas_ztrsv (CblasLower, CblasConjTrans, CblasNonUnit, cholesky, x);

      return GSL_SUCCESS;
    }
} /* gsl_linalg_complex_cholesky_solve() */

/*
gsl_linalg_complex_cholesky_svx()
  Solve A x = b in place where A is in cholesky form
*/

int
gsl_linalg_complex_cholesky_svx (const gsl_matrix_complex * cholesky,
                                 gsl_vector_complex * x)
{
  if (cholesky->size1 != cholesky->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (cholesky->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* solve for y using forward-substitution, L y = b */

      gsl_blas_ztrsv (CblasLower, CblasNoTrans, CblasNonUnit, cholesky, x);

      /* perform back-substitution, L^H x = y */

      gsl_blas_ztrsv (CblasLower, CblasConjTrans, CblasNonUnit, cholesky, x);

      return GSL_SUCCESS;
    }
} /* gsl_linalg_complex_cholesky_svx() */


/******************************************************************************

gsl_linalg_complex_cholesky_invert()
  Compute the inverse of an Hermitian positive definite matrix in
  Cholesky form.

Inputs: LLT - matrix in cholesky form on input
              A^{-1} = L^{-H} L^{-1} on output

Return: success or error
******************************************************************************/

int
gsl_linalg_complex_cholesky_invert(gsl_matrix_complex * LLT)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      size_t N = LLT->size1;
      size_t i, j;
      gsl_vector_complex_view v1;

      /* invert the lower triangle of LLT */
      for (i = 0; i < N; ++i)
        {
          double ajj;
          gsl_complex z;

          j = N - i - 1;

          { 
            gsl_complex z0 = gsl_matrix_complex_get(LLT, j, j);
            ajj = 1.0 / GSL_REAL(z0); 
          }

          GSL_SET_COMPLEX(&z, ajj, 0.0);
          gsl_matrix_complex_set(LLT, j, j, z);

          {
            gsl_complex z1 = gsl_matrix_complex_get(LLT, j, j);
            ajj = -GSL_REAL(z1);
          }

          if (j < N - 1)
            {
              gsl_matrix_complex_view m;
              
              m = gsl_matrix_complex_submatrix(LLT, j + 1, j + 1,
                                       N - j - 1, N - j - 1);
              v1 = gsl_matrix_complex_subcolumn(LLT, j, j + 1, N - j - 1);

              gsl_blas_ztrmv(CblasLower, CblasNoTrans, CblasNonUnit,
                             &m.matrix, &v1.vector);

              gsl_blas_zdscal(ajj, &v1.vector);
            }
        } /* for (i = 0; i < N; ++i) */

      /*
       * The lower triangle of LLT now contains L^{-1}. Now compute
       * A^{-1} = L^{-H} L^{-1}
       *
       * The (ij) element of A^{-1} is column i of conj(L^{-1}) dotted into
       * column j of L^{-1}
       */

      for (i = 0; i < N; ++i)
        {
          gsl_complex sum;
          for (j = i + 1; j < N; ++j)
            {
              gsl_vector_complex_view v2;
              v1 = gsl_matrix_complex_subcolumn(LLT, i, j, N - j);
              v2 = gsl_matrix_complex_subcolumn(LLT, j, j, N - j);

              /* compute Ainv[i,j] = sum_k{conj(Linv[k,i]) * Linv[k,j]} */
              gsl_blas_zdotc(&v1.vector, &v2.vector, &sum);

              /* store in upper triangle */
              gsl_matrix_complex_set(LLT, i, j, sum);
            }

          /* now compute the diagonal element */
          v1 = gsl_matrix_complex_subcolumn(LLT, i, i, N - i);
          gsl_blas_zdotc(&v1.vector, &v1.vector, &sum);
          gsl_matrix_complex_set(LLT, i, i, sum);
        }

      /* copy the Hermitian upper triangle to the lower triangle */

      for (j = 1; j < N; j++)
        {
          for (i = 0; i < j; i++)
            {
              gsl_complex z = gsl_matrix_complex_get(LLT, i, j);
              gsl_matrix_complex_set(LLT, j, i, gsl_complex_conjugate(z));
            }
        } 

      return GSL_SUCCESS;
    }
} /* gsl_linalg_complex_cholesky_invert() */


/********************************************
 *           INTERNAL ROUTINES              *
 ********************************************/

static void
cholesky_complex_conj_vector(gsl_vector_complex *v)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      gsl_complex z = gsl_vector_complex_get(v, i);
      gsl_vector_complex_set(v, i, gsl_complex_conjugate(z));
    }
} /* cholesky_complex_conj_vector() */
