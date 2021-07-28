/* linalg/qrc.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Copyright (C) 2017 Christian Krueger
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

/* Author:  G. Jungman, modified by C. Krueger */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/* Factorise a general complex-valued M x N matrix A into
 *
 *   A = Q R
 *
 * where Q is unitary (M x M) and R is upper triangular (M x N).
 *
 * Q is stored as a packed set of Householder transformations in the
 * strict lower triangular part of the input matrix.
 *
 * R is stored in the diagonal and upper triangle of the input matrix.
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
 *       v_i = [1, m(i+1,i), m(i+2,i), ... , m(M,i)]
 *
 * This storage scheme is the same as in LAPACK.  */

int
gsl_linalg_complex_QR_decomp (gsl_matrix_complex * A, gsl_vector_complex * tau)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau->size != N)
    {
      GSL_ERROR ("size of tau must be N", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      for (i = 0; i < GSL_MIN (M, N); i++)
        {
          /* Compute the Householder transformation to reduce the j-th
             column of the matrix to a multiple of the j-th unit vector */

          gsl_vector_complex_view c = gsl_matrix_complex_subcolumn (A, i, i, M - i);
          gsl_complex tau_i = gsl_linalg_complex_householder_transform (&(c.vector));

          gsl_vector_complex_set (tau, i, tau_i);

          /* apply the transformation to the remaining columns and update the norms */

          if (i + 1 < N)
            {
              gsl_matrix_complex_view m = gsl_matrix_complex_submatrix (A, i, i + 1, M - i, N - i - 1);
              gsl_complex tau_i_conj = gsl_complex_conjugate(tau_i);
              gsl_vector_complex_view work = gsl_vector_complex_subvector(tau, i + 1, N - i - 1);

              gsl_linalg_complex_householder_left(tau_i_conj, &(c.vector), &(m.matrix), &(work.vector));
            }
        }

      return GSL_SUCCESS;
    }
}

/* Solves the system A x = b using the QR factorisation,

 *  R x = Q^H b
 *
 * to obtain x.
 */

int
gsl_linalg_complex_QR_solve (const gsl_matrix_complex * QR, const gsl_vector_complex * tau,
                             const gsl_vector_complex * b, gsl_vector_complex * x)
{
  if (QR->size1 != QR->size2)
    {
      GSL_ERROR ("QR matrix must be square", GSL_ENOTSQR);
    }
  else if (QR->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (QR->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* copy x <- b */
      gsl_vector_complex_memcpy (x, b);

      /* solve for x */
      gsl_linalg_complex_QR_svx (QR, tau, x);

      return GSL_SUCCESS;
    }
}

/* Solves the system A x = b in place using the QR factorisation,

 *  R x = Q^H b
 *
 * to obtain x.
 */

int
gsl_linalg_complex_QR_svx (const gsl_matrix_complex * QR, const gsl_vector_complex * tau, gsl_vector_complex * x)
{

  if (QR->size1 != QR->size2)
    {
      GSL_ERROR ("QR matrix must be square", GSL_ENOTSQR);
    }
  else if (QR->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match x/rhs size", GSL_EBADLEN);
    }
  else
    {
      /* compute rhs = Q^H b */
      gsl_linalg_complex_QR_QHvec (QR, tau, x);

      /* solve R x = rhs, storing x in-place */
      gsl_blas_ztrsv (CblasUpper, CblasNoTrans, CblasNonUnit, QR, x);

      return GSL_SUCCESS;
    }
}


/* Find the least squares solution to the overdetermined system
 *
 *   A x = b
 *
 * for M >= N using the QR factorization A = Q R.
 */

int
gsl_linalg_complex_QR_lssolve (const gsl_matrix_complex * QR, const gsl_vector_complex * tau,
                               const gsl_vector_complex * b, gsl_vector_complex * x,
                               gsl_vector_complex * residual)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (M < N)
    {
      GSL_ERROR ("QR matrix must have M>=N", GSL_EBADLEN);
    }
  else if (M != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (N != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (M != residual->size)
    {
      GSL_ERROR ("matrix size must match residual size", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_complex_const_view R = gsl_matrix_complex_const_submatrix (QR, 0, 0, N, N);
      gsl_vector_complex_view c = gsl_vector_complex_subvector(residual, 0, N);

      gsl_vector_complex_memcpy(residual, b);

      /* compute rhs = Q^H b */
      gsl_linalg_complex_QR_QHvec (QR, tau, residual);

      /* solve R x = rhs */
      gsl_vector_complex_memcpy(x, &(c.vector));
      gsl_blas_ztrsv (CblasUpper, CblasNoTrans, CblasNonUnit, &(R.matrix), x);

      /* compute residual = b - A x = Q (Q^H b - R x) */
      gsl_vector_complex_set_zero(&(c.vector));
      gsl_linalg_complex_QR_Qvec(QR, tau, residual);

      return GSL_SUCCESS;
    }
}

/* form the product v := Q^H v from a QR factorized matrix */

int
gsl_linalg_complex_QR_QHvec (const gsl_matrix_complex * QR, const gsl_vector_complex * tau, gsl_vector_complex * v)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (tau->size != N)
    {
      GSL_ERROR ("size of tau must be N", GSL_EBADLEN);
    }
  else if (v->size != M)
    {
      GSL_ERROR ("vector size must be M", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* compute Q^H v */

      for (i = 0; i < GSL_MIN (M, N); i++)
        {
          gsl_vector_complex_const_view h = gsl_matrix_complex_const_subcolumn (QR, i, i, M - i);
          gsl_vector_complex_view w = gsl_vector_complex_subvector (v, i, M - i);
          gsl_complex ti = gsl_vector_complex_get (tau, i);
          gsl_complex ti_conj = gsl_complex_conjugate(ti);
          gsl_linalg_complex_householder_hv (ti_conj, &(h.vector), &(w.vector));
        }

      return GSL_SUCCESS;
    }
}


int
gsl_linalg_complex_QR_Qvec (const gsl_matrix_complex * QR, const gsl_vector_complex * tau, gsl_vector_complex * v)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (v->size != M)
    {
      GSL_ERROR ("vector size must be M", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* compute Q v */

      for (i = GSL_MIN (M, N); i-- > 0;)
        {
          gsl_vector_complex_const_view c = gsl_matrix_complex_const_column (QR, i);
          gsl_vector_complex_const_view h = gsl_vector_complex_const_subvector (&(c.vector),
                                                                i, M - i);
          gsl_vector_complex_view w = gsl_vector_complex_subvector (v, i, M - i);
          gsl_complex ti = gsl_vector_complex_get (tau, i);
          /* we do not need the conjugate of ti here */
          gsl_linalg_complex_householder_hv (ti, &h.vector, &w.vector);
        }
      return GSL_SUCCESS;
    }
}

/* form the unitary matrix Q from the packed QR matrix */
int
gsl_linalg_complex_QR_unpack (const gsl_matrix_complex * QR, const gsl_vector_complex * tau,
                              gsl_matrix_complex * Q, gsl_matrix_complex * R)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (Q->size1 != M || Q->size2 != M)
    {
      GSL_ERROR ("Q matrix must be M x M", GSL_ENOTSQR);
    }
  else if (R->size1 != M || R->size2 != N)
    {
      GSL_ERROR ("R matrix must be M x N", GSL_ENOTSQR);
    }
  else if (tau->size != N)
    {
      GSL_ERROR ("size of tau must be N", GSL_EBADLEN);
    }
  else
    {
      size_t i, j;

      /* initialize Q to the identity */
      gsl_matrix_complex_set_identity (Q);

      for (i = GSL_MIN (M, N); i-- > 0;)
        {
          gsl_vector_complex_const_view c = gsl_matrix_complex_const_column (QR, i);
          gsl_vector_complex_const_view h = gsl_vector_complex_const_subvector (&c.vector, i, M - i);
          gsl_matrix_complex_view m = gsl_matrix_complex_submatrix (Q, i, i, M - i, M - i);
          gsl_complex ti = gsl_vector_complex_get (tau, i);
          gsl_vector_complex_view work = gsl_matrix_complex_subcolumn(R, 0, 0, M - i);

          /* we do not need the conjugate of ti here */
          gsl_linalg_complex_householder_left (ti, &h.vector, &m.matrix, &work.vector);
        }

      /* form the right triangular matrix R from a packed QR matrix */

      for (i = 0; i < M; i++)
        {
          for (j = 0; j < i && j < N; j++)
            gsl_matrix_complex_set (R, i, j, GSL_COMPLEX_ZERO);

          for (j = i; j < N; j++)
            gsl_matrix_complex_set (R, i, j, gsl_matrix_complex_get (QR, i, j));
        }

      return GSL_SUCCESS;
    }
}
