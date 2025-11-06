/* linalg/lu.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Gerard Jungman, Brian Gough
 * Copyright (C) 2019, 2021 Patrick Alken
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

#include "recurse.h"

static int LU_decomp_L2 (gsl_matrix * A, gsl_vector_uint * ipiv);
static int LU_decomp_L3 (gsl_matrix * A, gsl_vector_uint * ipiv);
static int singular (const gsl_matrix * LU);
static int apply_pivots(gsl_matrix * A, const gsl_vector_uint * ipiv);

/* Factorise a general N x N matrix A into,
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
 * signum gives the sign of the permutation, (-1)^n, where n is the
 * number of interchanges in the permutation. 
 *
 * See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss
 * Elimination with Partial Pivoting).
 */

int
gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int *signum)
{
  const size_t M = A->size1;

  if (p->size != M)
    {
      GSL_ERROR ("permutation length must match matrix size1", GSL_EBADLEN);
    }
  else
    {
      int status;
      const size_t N = A->size2;
      const size_t minMN = GSL_MIN(M, N);
      gsl_vector_uint * ipiv = gsl_vector_uint_alloc(minMN);
      gsl_matrix_view AL = gsl_matrix_submatrix(A, 0, 0, M, minMN);
      size_t i;

      status = LU_decomp_L3 (&AL.matrix, ipiv);

      /* process remaining right matrix */
      if (M < N)
        {
          gsl_matrix_view AR = gsl_matrix_submatrix(A, 0, M, M, N - M);

          /* apply pivots to AR */
          apply_pivots(&AR.matrix, ipiv);

          /* AR = AL^{-1} AR */
          gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 1.0, &AL.matrix, &AR.matrix);
        }

      /* convert ipiv array to permutation */

      gsl_permutation_init(p);
      *signum = 1;

      for (i = 0; i < minMN; ++i)
        {
          unsigned int pivi = gsl_vector_uint_get(ipiv, i);

          if (p->data[pivi] != p->data[i])
            {
              size_t tmp = p->data[pivi];
              p->data[pivi] = p->data[i];
              p->data[i] = tmp;
              *signum = -(*signum);
            }
        }

      gsl_vector_uint_free(ipiv);

      return status;
    }
}

/*
LU_decomp_L2
  LU decomposition with partial pivoting using Level 2 BLAS

Inputs: A    - on input, matrix to be factored; on output, L and U factors
        ipiv - (output) array containing row swaps

Notes:
1) Based on LAPACK DGETF2

Return: GSL_SUCCESS if successful; otherwise a positive integer k
in [1,GSL_MIN(M,N)] indicating that U(k,k) is zero and therefore singular;
in this case, the LU factorization is still completed
*/

static int
LU_decomp_L2 (gsl_matrix * A, gsl_vector_uint * ipiv)
{
  int status = GSL_SUCCESS;
  const size_t M = A->size1;
  const size_t N = A->size2;
  const size_t minMN = GSL_MIN(M, N);

  if (ipiv->size != minMN)
    {
      GSL_ERROR ("ipiv length must equal MIN(M,N)", GSL_EBADLEN);
    }
  else
    {
      size_t i, j;

      for (j = 0; j < minMN; ++j)
        {
          /* find maximum in the j-th column */
          gsl_vector_view v = gsl_matrix_subcolumn(A, j, j, M - j);
          size_t j_pivot = j + gsl_blas_idamax(&v.vector);
          double Ajpj = gsl_matrix_get(A, j_pivot, j);
          gsl_vector_view v1, v2;

          gsl_vector_uint_set(ipiv, j, j_pivot);

          if (Ajpj != 0.0)
            {
              if (j_pivot != j)
                {
                  /* swap rows j and j_pivot */
                  v1 = gsl_matrix_row(A, j);
                  v2 = gsl_matrix_row(A, j_pivot);
                  gsl_blas_dswap(&v1.vector, &v2.vector);
                }

              if (j < M - 1)
                {
                  double Ajj = gsl_matrix_get(A, j, j);

                  if (fabs(Ajj) >= GSL_DBL_MIN)
                    {
                      v1 = gsl_matrix_subcolumn(A, j, j + 1, M - j - 1);
                      gsl_blas_dscal(1.0 / Ajj, &v1.vector);
                    }
                  else
                    {
                      for (i = 1; i < M - j; ++i)
                        {
                          double * ptr = gsl_matrix_ptr(A, j + i, j);
                          *ptr /= Ajj;
                        }
                    }
                }
            }
          else
            {
              status = (int) j + 1;
            }

          if (j < minMN - 1)
            {
              gsl_matrix_view A22 = gsl_matrix_submatrix(A, j + 1, j + 1, M - j - 1, N - j - 1);
              v1 = gsl_matrix_subcolumn(A, j, j + 1, M - j - 1);
              v2 = gsl_matrix_subrow(A, j, j + 1, N - j - 1);

              gsl_blas_dger(-1.0, &v1.vector, &v2.vector, &A22.matrix);
            }
        }

      return status;
    }
}

/*
LU_decomp_L3
  LU decomposition with partial pivoting using Level 3 BLAS

Inputs: A    - on input, matrix to be factored; on output, L and U factors
        ipiv - (output) array containing row swaps

Notes:
1) Based on ReLAPACK DGETRF
*/

static int
LU_decomp_L3 (gsl_matrix * A, gsl_vector_uint * ipiv)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M < N)
    {
      GSL_ERROR ("matrix must have M >= N", GSL_EBADLEN);
    }
  else if (ipiv->size != GSL_MIN(M, N))
    {
      GSL_ERROR ("ipiv length must equal MIN(M,N)", GSL_EBADLEN);
    }
  else if (N <= CROSSOVER_LU)
    {
      /* use Level 2 algorithm */
      return LU_decomp_L2(A, ipiv);
    }
  else
    {
      /*
       * partition matrix:
       *
       *       N1  N2
       * N1  [ A11 A12 ]
       * M2  [ A21 A22 ]
       *
       * and
       *      N1  N2
       * M  [ AL  AR  ]
       */
      int status;
      const size_t N1 = GSL_LINALG_SPLIT(N);
      const size_t N2 = N - N1;
      const size_t M2 = M - N1;
      gsl_matrix_view A11 = gsl_matrix_submatrix(A, 0, 0, N1, N1);
      gsl_matrix_view A12 = gsl_matrix_submatrix(A, 0, N1, N1, N2);
      gsl_matrix_view A21 = gsl_matrix_submatrix(A, N1, 0, M2, N1);
      gsl_matrix_view A22 = gsl_matrix_submatrix(A, N1, N1, M2, N2);

      gsl_matrix_view AL = gsl_matrix_submatrix(A, 0, 0, M, N1);
      gsl_matrix_view AR = gsl_matrix_submatrix(A, 0, N1, M, N2);

      /*
       * partition ipiv = [ ipiv1 ] N1
       *                  [ ipiv2 ] N2
       */
      gsl_vector_uint_view ipiv1 = gsl_vector_uint_subvector(ipiv, 0, N1);
      gsl_vector_uint_view ipiv2 = gsl_vector_uint_subvector(ipiv, N1, N2);

      size_t i;

      /* recursion on (AL, ipiv1) */
      status = LU_decomp_L3(&AL.matrix, &ipiv1.vector);
      if (status)
        return status;

      /* apply ipiv1 to AR */
      apply_pivots(&AR.matrix, &ipiv1.vector);

      /* A12 = A11^{-1} A12 */
      gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 1.0, &A11.matrix, &A12.matrix);

      /* A22 = A22 - A21 * A12 */
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &A21.matrix, &A12.matrix, 1.0, &A22.matrix);

      /* recursion on (A22, ipiv2) */
      status = LU_decomp_L3(&A22.matrix, &ipiv2.vector);
      if (status)
        return status;

      /* apply pivots to A21 */
      apply_pivots(&A21.matrix, &ipiv2.vector);

      /* shift pivots */
      for (i = 0; i < N2; ++i)
        {
          unsigned int * ptr = gsl_vector_uint_ptr(&ipiv2.vector, i);
          *ptr += N1;
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_LU_solve (const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LU->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (singular (LU)) 
    {
      GSL_ERROR ("matrix is singular", GSL_EDOM);
    }
  else
    {
      int status;

      /* copy x <- b */
      gsl_vector_memcpy (x, b);

      /* solve for x */
      status = gsl_linalg_LU_svx (LU, p, x);

      return status;
    }
}


int
gsl_linalg_LU_svx (const gsl_matrix * LU, const gsl_permutation * p, gsl_vector * x)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution/rhs size", GSL_EBADLEN);
    }
  else if (singular (LU)) 
    {
      GSL_ERROR ("matrix is singular", GSL_EDOM);
    }
  else
    {
      /* apply permutation to RHS */
      gsl_permute_vector (p, x);

      /* solve for c using forward-substitution, L c = P b */
      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasUnit, LU, x);

      /* perform back-substitution, U x = c */
      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LU, x);

      return GSL_SUCCESS;
    }
}


int
gsl_linalg_LU_refine (const gsl_matrix * A, const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x, gsl_vector * work)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("matrix a must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (A->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be decomposition of a", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LU->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (LU->size1 != work->size)
    {
      GSL_ERROR ("matrix size must match workspace size", GSL_EBADLEN);
    }
  else if (singular (LU)) 
    {
      GSL_ERROR ("matrix is singular", GSL_EDOM);
    }
  else
    {
      int status;

      /* compute residual = (A * x  - b) */
      gsl_vector_memcpy (work, b);
      gsl_blas_dgemv (CblasNoTrans, 1.0, A, x, -1.0, work);

      /* find correction, delta = - (A^-1) * residual, and apply it */
      status = gsl_linalg_LU_svx (LU, p, work);
      gsl_blas_daxpy (-1.0, work, x);

      return status;
    }
}

int
gsl_linalg_LU_invert (const gsl_matrix * LU, const gsl_permutation * p, gsl_matrix * inverse)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (inverse->size1 != LU->size1 || inverse->size2 != LU->size2)
    {
      GSL_ERROR ("inverse matrix must match LU matrix dimensions", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_memcpy(inverse, LU);
      return gsl_linalg_LU_invx (inverse, p);
    }
}

int
gsl_linalg_LU_invx (gsl_matrix * LU, const gsl_permutation * p)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (singular (LU)) 
    {
      GSL_ERROR ("matrix is singular", GSL_EDOM);
    }
  else
    {
      int status;
      const size_t N = LU->size1;
      size_t i;

      /* compute U^{-1} */
      status = gsl_linalg_tri_invert(CblasUpper, CblasNonUnit, LU);
      if (status)
        return status;

      /* compute L^{-1} */
      status = gsl_linalg_tri_invert(CblasLower, CblasUnit, LU);
      if (status)
        return status;

      /* compute U^{-1} L^{-1} */
      status = gsl_linalg_tri_UL(LU);
      if (status)
        return status;

      /* apply permutation to columns of A^{-1} */
      for (i = 0; i < N; ++i)
        {
          gsl_vector_view v = gsl_matrix_row(LU, i);
          gsl_permute_vector_inverse(p, &v.vector);
        }

      return GSL_SUCCESS;
    }
}

double
gsl_linalg_LU_det (gsl_matrix * LU, int signum)
{
  size_t i, n = LU->size1;

  double det = (double) signum;

  for (i = 0; i < n; i++)
    {
      det *= gsl_matrix_get (LU, i, i);
    }

  return det;
}


double
gsl_linalg_LU_lndet (gsl_matrix * LU)
{
  size_t i, n = LU->size1;

  double lndet = 0.0;

  for (i = 0; i < n; i++)
    {
      lndet += log (fabs (gsl_matrix_get (LU, i, i)));
    }

  return lndet;
}

int
gsl_linalg_LU_sgndet (gsl_matrix * LU, int signum)
{
  size_t i, n = LU->size1;

  int s = signum;

  for (i = 0; i < n; i++)
    {
      double u = gsl_matrix_get (LU, i, i);

      if (u < 0)
        {
          s *= -1;
        }
      else if (u == 0)
        {
          s = 0;
          break;
        }
    }

  return s;
}

static int
singular (const gsl_matrix * LU)
{
  size_t i, n = LU->size1;

  for (i = 0; i < n; i++)
    {
      double u = gsl_matrix_get (LU, i, i);
      if (u == 0) return 1;
    }
 
  return 0;
}

static int
apply_pivots(gsl_matrix * A, const gsl_vector_uint * ipiv)
{
  if (A->size1 < ipiv->size)
    {
      GSL_ERROR("matrix does not match pivot vector", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      for (i = 0; i < ipiv->size; ++i)
        {
          size_t pi = gsl_vector_uint_get(ipiv, i);

          if (i != pi)
            {
              /* swap rows i and pi */
              gsl_vector_view v1 = gsl_matrix_row(A, i);
              gsl_vector_view v2 = gsl_matrix_row(A, pi);
              gsl_blas_dswap(&v1.vector, &v2.vector);
            }
        }

      return GSL_SUCCESS;
    }
}
