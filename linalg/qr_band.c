/* linalg/qr_band.c
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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/* Factorise a (p,q) banded M x N matrix A into
 *  
 *   A = Q R
 *
 * where Q is orthogonal (M x M) and R is upper triangular (M x N).
 *
 * Example with M = 7, N = 6, (p,q) = (2,1)
 *
 * A = [ A11 A12  0   0   0   0  ]
 *     [ A21 A22 A23  0   0   0  ]
 *     [ A31 A32 A33 A34  0   0  ]
 *     [  0  A42 A43 A44 A45  0  ]
 *     [  0   0  A53 A54 A55 A56 ]
 *     [  0   0   0  A64 A65 A66 ]
 *     [  0   0   0   0  A75 A76 ]
 *
 * AB has dimensions N-by-(2p + q + 1)
 *
 * INPUT:                          OUTPUT:
 *
 * AB = [ *  *  *  A11 A21 A31 ]   AB = [  *   *   *  R11 V21 V31 ]
 *      [ *  * A12 A22 A32 A42 ]        [  *   *  R12 R22 V32 V42 ]
 *      [ *  0 A23 A33 A43 A53 ]        [  *  R13 R23 R33 V43 V53 ]
 *      [ 0  0 A34 A44 A54 A64 ]        [ R14 R24 R34 R44 V54 V64 ]
 *      [ 0  0 A45 A55 A65 A75 ]        [ R25 R35 R45 R55 V65 V75 ]
 *      [ 0  0 A56 A66 A76  *  ]        [ R36 R46 R56 R66 V76  *  ]
 *        -p-- -q- -1- ---p---            ---p--- -q- -1- ---p---
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
 */

int
gsl_linalg_QR_band_decomp_L2 (const size_t M, const size_t p, const size_t q, gsl_matrix * AB, gsl_vector * tau)
{
  const size_t N = AB->size1;

  if (tau->size != N)
    {
      GSL_ERROR ("tau must have length N", GSL_EBADLEN);
    }
  else if (AB->size2 != 2*p + q + 1)
    {
      GSL_ERROR ("dimensions of AB are inconsistent with (p,q)", GSL_EBADLEN);
    }
  else
    {
      const size_t minMN = GSL_MIN(M, N);
      size_t j;

      /* set AB(:,1:p) to zero */
      if (p > 0)
        {
          gsl_matrix_view m = gsl_matrix_submatrix(AB, 0, 0, N, p);
          gsl_matrix_set_zero(&m.matrix);
        }

      for (j = 0; j < minMN; ++j)
        {
          /* Compute the Householder transformation to reduce the j-th
             column of the matrix to a multiple of the j-th unit vector */

          size_t k1 = GSL_MIN(p + 1, M - j);     /* number of non-zero elements of this column, including diagonal element */
          size_t k2 = GSL_MIN(p + q, N - j - 1); /* number of columns to update */
          gsl_vector_view c = gsl_matrix_subrow(AB, j, p + q, k1);
          double tau_j = gsl_linalg_householder_transform (&(c.vector));
          double * ptr = gsl_vector_ptr(&(c.vector), 0);

          gsl_vector_set (tau, j, tau_j);

          /* apply the transformation to the remaining columns */
          if (k2 > 0)
            {
              gsl_matrix_view m = gsl_matrix_submatrix (AB, j + 1, p + q - 1, k2, k1);
              gsl_vector_view work = gsl_vector_subvector(tau, j + 1, k2);
              double tmp = *ptr;

              m.matrix.tda -= 1; /* unskew matrix */

              /* we want to compute H*A(j:j+k1-1,j+1:j+k2), but due to our using row-major order, the
               * matrix m contains A(j:j+k1-1,j+1:j+k2)^T. So therefore we apply H from the right,
               *
               * [H*A]^T = A^T H^T = A^T H
               */

              *ptr = 1.0;
              gsl_linalg_householder_right(tau_j, &(c.vector), &(m.matrix), &(work.vector));
              *ptr = tmp;
            }
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_QR_band_unpack_L2 (const size_t p, const size_t q, const gsl_matrix * QRB, const gsl_vector * tau,
                              gsl_matrix * Q, gsl_matrix * R)
{
  const size_t M = Q->size1;
  const size_t N = QRB->size1;

  if (Q->size2 != M)
    {
      GSL_ERROR ("Q matrix must be square", GSL_ENOTSQR);
    }
  else if (R->size1 != M || R->size2 != N)
    {
      GSL_ERROR ("R matrix must be M x N", GSL_ENOTSQR);
    }
  else if (tau->size < GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be at least MIN(M,N)", GSL_EBADLEN);
    }
  else if (QRB->size2 != 2*p + q + 1)
    {
      GSL_ERROR ("dimensions of QRB are inconsistent with (p,q)", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* form matrix Q */
      gsl_matrix_set_identity (Q);

      for (i = GSL_MIN (M, N); i-- > 0;)
        {
          size_t k1 = GSL_MIN(p + 1, M - i);     /* number of non-zero elements of this column, including diagonal element */
          gsl_vector_const_view h = gsl_matrix_const_subrow(QRB, i, p + q, k1);
          gsl_matrix_view m = gsl_matrix_submatrix (Q, i, i, k1, M - i);
          double ti = gsl_vector_get (tau, i);
          gsl_vector_view work = gsl_matrix_subcolumn(R, 0, 0, M - i);
          double * ptr = gsl_vector_ptr((gsl_vector *) &h.vector, 0);
          double tmp = *ptr;

          *ptr = 1.0;
          gsl_linalg_householder_left (ti, &h.vector, &m.matrix, &work.vector);
          *ptr = tmp;
        }

      /* form matrix R */
      gsl_matrix_set_zero(R);

      for (i = 0; i <= GSL_MIN(p + q, N - 1); ++i)
        {
          gsl_vector_const_view src = gsl_matrix_const_subcolumn(QRB, p + q - i, i, GSL_MIN(M, N - i));
          gsl_vector_view dest = gsl_matrix_superdiagonal(R, i);
          gsl_vector_memcpy(&dest.vector, &src.vector);
        }

      return GSL_SUCCESS;
    }
}
