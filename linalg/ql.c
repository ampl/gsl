/* linalg/ql.c
 * 
 * Copyright (C) 2019 Patrick Alken
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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/* Factorise a general M x N matrix A into
 *  
 *   A = Q L
 *
 * where Q is orthogonal (M x M) and L is lower triangular (M x N).
 *
 * Q is stored as a packed set of Householder transformations in the
 * strict upper triangular part of the input matrix.
 *
 * L is stored in the diagonal and lower triangle of the input matrix.
 *
 * The full matrix for Q can be obtained as the product
 *
 *       Q = Q_k .. Q_2 Q_1
 *
 * where k = MIN(M,N) and
 *
 *       Q_i = (I - tau_i * v_i * v_i')
 *
 * and where v_i is a Householder vector
 *
 *       v_i = [A(1,N-k+i), A(2,N-k+i), ... , A(M-k+i,N-k+i), 1, 0, ..., 0]
 *
 * This storage scheme is the same as in LAPACK.  */

/*
gsl_linalg_QL_decomp()
  Perform QL decomposition of a matrix A

Inputs: A   - M-by-N matrix
        tau - (output) Householder coefficients, length N

Notes:
1) The K = MIN(M, N) Householder scalars are stored in tau(N-K+1:N)
on output; the rest of tau is used as temporary workspace
*/

int
gsl_linalg_QL_decomp (gsl_matrix * A, gsl_vector * tau)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau->size != N)
    {
      GSL_ERROR ("size of tau must be N", GSL_EBADLEN);
    }
  else
    {
      const size_t K = GSL_MIN(M, N);
      size_t i;

      for (i = 0; i < K; i++)
        {
          /* compute the Householder transformation to annihilate the (N-K+i)-th
             column of the matrix */

          gsl_vector_view c = gsl_matrix_subcolumn (A, N - i - 1, 0, M - i);
          double * alpha = gsl_matrix_ptr(A, M - i - 1, N - i - 1);
          double tau_j = gsl_linalg_householder_transform2 (alpha, &(c.vector));

          /* apply the transformation to A(1:M-i,1:N-i-2) from the left */
          if (i + 1 < N)
            {
              gsl_vector_view work = gsl_vector_subvector(tau, 0, N - i - 1);
              gsl_matrix_view m = gsl_matrix_submatrix (A, 0, 0, M - i, N - i - 1);
              double tmp = *alpha;

              *alpha = 1.0;
              gsl_linalg_householder_left (tau_j, &(c.vector), &(m.matrix), &work.vector);
              *alpha = tmp;
            }

          gsl_vector_set (tau, N - i - 1, tau_j);
        }

      return GSL_SUCCESS;
    }
}

/* form the orthogonal matrix Q from the packed QL matrix */
int
gsl_linalg_QL_unpack (const gsl_matrix * QL, const gsl_vector * tau, gsl_matrix * Q, gsl_matrix * L)
{
  const size_t M = QL->size1;
  const size_t N = QL->size2;

  if (Q->size1 != M || Q->size2 != M)
    {
      GSL_ERROR ("Q matrix must be M x M", GSL_ENOTSQR);
    }
  else if (L->size1 != M || L->size2 != N)
    {
      GSL_ERROR ("L matrix must be M x N", GSL_ENOTSQR);
    }
  else if (tau->size != N)
    {
      GSL_ERROR ("size of tau must be N", GSL_EBADLEN);
    }
  else
    {
      const size_t K = GSL_MIN(M, N);
      size_t i;

      /* initialize Q to the identity */
      gsl_matrix_set_identity (Q);

      for (i = 0; i < K; ++i)
        {
          gsl_vector_const_view h = gsl_matrix_const_subcolumn (QL, N - K + i, 0, M - K + i + 1);
          gsl_matrix_view m = gsl_matrix_submatrix (Q, 0, 0, M - K + i + 1, M - K + i + 1);
          gsl_vector_view work = gsl_matrix_subcolumn(L, 0, 0, M - K + i + 1);
          double ti = gsl_vector_get (tau, N - K + i);
          double * ptr = gsl_matrix_ptr((gsl_matrix *) QL, M - K + i, N - K + i);
          double tmp = *ptr;

          *ptr = 1.0;
          gsl_linalg_householder_left (ti, &h.vector, &m.matrix, &work.vector);
          *ptr = tmp;
        }

      /* form the left triangular matrix L from a packed QL matrix */
      gsl_matrix_set_zero(L);

      if (M >= N)
        {
          gsl_matrix_const_view src = gsl_matrix_const_submatrix(QL, M - N, 0, N, N);
          gsl_matrix_view dest = gsl_matrix_submatrix(L, M - N, 0, N, N);
          gsl_matrix_tricpy(CblasLower, CblasNonUnit, &dest.matrix, &src.matrix);
        }
      else
        {
          gsl_matrix_const_view src1 = gsl_matrix_const_submatrix(QL, 0, 0, M, N - M);
          gsl_matrix_view dest1 = gsl_matrix_submatrix(L, 0, 0, M, N - M);

          gsl_matrix_const_view src2 = gsl_matrix_const_submatrix(QL, 0, N - M, M, M);
          gsl_matrix_view dest2 = gsl_matrix_submatrix(L, 0, N - M, M, M);

          gsl_matrix_memcpy(&dest1.matrix, &src1.matrix);
          gsl_matrix_tricpy(CblasLower, CblasNonUnit, &dest2.matrix, &src2.matrix);
        }

      return GSL_SUCCESS;
    }
}
