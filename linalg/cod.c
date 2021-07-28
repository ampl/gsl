/* linalg/cod.c
 * 
 * Copyright (C) 2016, 2017 Patrick Alken
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

/*
 * This module contains routines for factoring an M-by-N matrix A as:
 *
 * A P = Q R Z^T
 *
 * known as the Complete Orthogonal Decomposition, where:
 *
 * P is a N-by-N permutation matrix
 * Q is M-by-M orthogonal
 * R has an r-by-r upper triangular block
 * Z is N-by-N orthogonal
 *
 * When A is full rank, Z = I and this becomes the QR decomposition
 * with column pivoting. When A is rank deficient, then
 *
 * R = [ R11 0 ] where R11 is r-by-r and r = rank(A)
 *     [  0  0 ]
 */

static int cod_RZ(gsl_matrix * A, gsl_vector * tau);
static double cod_householder_transform(double *alpha, gsl_vector * v);
static int cod_householder_mh(const double tau, const gsl_vector * v,
                              gsl_matrix * A, gsl_vector * work);
static int cod_householder_hv(const double tau, const gsl_vector * v, gsl_vector * w);
static int cod_householder_Zvec(const gsl_matrix * QRZT, const gsl_vector * tau_Z, const size_t rank,
                                gsl_vector * v);
static int cod_trireg_solve(const gsl_matrix * R, const double lambda, const gsl_vector * b,
                            gsl_matrix * S, gsl_vector * x, gsl_vector * work);

int
gsl_linalg_COD_decomp_e(gsl_matrix * A, gsl_vector * tau_Q, gsl_vector * tau_Z,
                        gsl_permutation * p, double tol, size_t * rank, gsl_vector * work)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau_Q->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau_Q must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (tau_Z->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau_Z must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (p->size != N)
    {
      GSL_ERROR ("permutation size must be N", GSL_EBADLEN);
    }
  else if (work->size != N)
    {
      GSL_ERROR ("work size must be N", GSL_EBADLEN);
    }
  else
    {
      int status, signum;
      size_t r;

      /* decompose: A P = Q R */
      status = gsl_linalg_QRPT_decomp(A, tau_Q, p, &signum, work);
      if (status)
        return status;

      /* estimate rank of A */
      r = gsl_linalg_QRPT_rank(A, tol);

      if (r < N)
        {
          /*
           * matrix is rank-deficient, so that the R factor is
           *
           * R = [ R11 R12 ] =~ [ R11 R12 ]
           *     [  0  R22 ]    [  0   0  ]
           *
           * compute RZ decomposition of upper trapezoidal matrix
           * [ R11 R12 ] = [ R11~ 0 ] Z
           */
          gsl_matrix_view R_upper = gsl_matrix_submatrix(A, 0, 0, r, N);
          gsl_vector_view t = gsl_vector_subvector(tau_Z, 0, r);

          cod_RZ(&R_upper.matrix, &t.vector);
        }

      *rank = r;

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_COD_decomp(gsl_matrix * A, gsl_vector * tau_Q, gsl_vector * tau_Z,
                      gsl_permutation * p, size_t * rank, gsl_vector * work)
{
  return gsl_linalg_COD_decomp_e(A, tau_Q, tau_Z, p, -1.0, rank, work);
}

/*
gsl_linalg_COD_lssolve()
  Find the least squares solution to the overdetermined system 

min ||b - A x||^2

for M >= N using the COD factorization A P = Q R Z

Inputs: QRZT     - matrix A, in COD compressed format, M-by-N
        tau_Q    - Householder scalars for Q, length min(M,N)
        tau_Z    - Householder scalars for Z, length min(M,N)
        perm     - permutation matrix
        rank     - rank of A
        b        - rhs vector, length M
        x        - (output) solution vector, length N
        residual - (output) residual vector, b - A x, length M
*/

int
gsl_linalg_COD_lssolve (const gsl_matrix * QRZT, const gsl_vector * tau_Q, const gsl_vector * tau_Z,
                        const gsl_permutation * perm, const size_t rank, const gsl_vector * b,
                        gsl_vector * x, gsl_vector * residual)
{
  const size_t M = QRZT->size1;
  const size_t N = QRZT->size2;

  if (M < N)
    {
      GSL_ERROR ("QRZT matrix must have M>=N", GSL_EBADLEN);
    }
  else if (M != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (rank > GSL_MIN (M, N))
    {
      GSL_ERROR ("rank must be <= MIN(M,N)", GSL_EBADLEN);
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
      gsl_matrix_const_view R11 = gsl_matrix_const_submatrix (QRZT, 0, 0, rank, rank);
      gsl_vector_view QTb1 = gsl_vector_subvector(residual, 0, rank);
      gsl_vector_view x1 = gsl_vector_subvector(x, 0, rank);

      gsl_vector_set_zero(x);

      /* compute residual = Q^T b = [ c1 ; c2 ] */
      gsl_vector_memcpy(residual, b);
      gsl_linalg_QR_QTvec (QRZT, tau_Q, residual);

      /* solve x1 := R11^{-1} (Q^T b)(1:r) */
      gsl_vector_memcpy(&(x1.vector), &(QTb1.vector));
      gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, &(R11.matrix), &(x1.vector));

      /* compute Z ( R11^{-1} x1; 0 ) */
      cod_householder_Zvec(QRZT, tau_Z, rank, x);

      /* compute x = P Z^T ( R11^{-1} x1; 0 ) */
      gsl_permute_vector_inverse(perm, x);

      /* compute residual = b - A x = Q (Q^T b - R [ R11^{-1} x1; 0 ]) = Q [ 0 ; c2 ] */
      gsl_vector_set_zero(&(QTb1.vector));
      gsl_linalg_QR_Qvec(QRZT, tau_Q, residual);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_COD_lssolve2()
  Find the least squares solution to the Tikhonov regularized
system in standard form:

min ||b - A x||^2 + lambda^2 ||x||^2

for M >= N using the COD factorization A P = Q R Z

Inputs: lambda   - parameter
        QRZT     - matrix A, in COD compressed format, M-by-N
        tau_Q    - Householder scalars for Q, length min(M,N)
        tau_Z    - Householder scalars for Z, length min(M,N)
        perm     - permutation matrix
        rank     - rank of A
        b        - rhs vector, length M
        x        - (output) solution vector, length N
        residual - (output) residual vector, b - A x, length M
        S        - workspace, rank-by-rank
        work     - workspace, length rank
*/

int
gsl_linalg_COD_lssolve2 (const double lambda, const gsl_matrix * QRZT, const gsl_vector * tau_Q, const gsl_vector * tau_Z,
                         const gsl_permutation * perm, const size_t rank, const gsl_vector * b,
                         gsl_vector * x, gsl_vector * residual, gsl_matrix * S, gsl_vector * work)
{
  const size_t M = QRZT->size1;
  const size_t N = QRZT->size2;

  if (M < N)
    {
      GSL_ERROR ("QRZT matrix must have M>=N", GSL_EBADLEN);
    }
  else if (M != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (rank > GSL_MIN (M, N))
    {
      GSL_ERROR ("rank must be <= MIN(M,N)", GSL_EBADLEN);
    }
  else if (N != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (M != residual->size)
    {
      GSL_ERROR ("matrix size must match residual size", GSL_EBADLEN);
    }
  else if (S->size1 != rank || S->size2 != rank)
    {
      GSL_ERROR ("S must be rank-by-rank", GSL_EBADLEN);
    }
  else if (work->size != rank)
    {
      GSL_ERROR ("work must be length rank", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_const_view R11 = gsl_matrix_const_submatrix (QRZT, 0, 0, rank, rank);
      gsl_vector_view c1 = gsl_vector_subvector(residual, 0, rank);
      gsl_vector_view y1 = gsl_vector_subvector(x, 0, rank);

      gsl_vector_set_zero(x);

      /* compute residual = Q^T b = [ c1 ; c2 ]*/
      gsl_vector_memcpy(residual, b);
      gsl_linalg_QR_QTvec (QRZT, tau_Q, residual);

      /* solve [ R11 ; lambda*I ] y1 = [ (Q^T b)(1:r) ; 0 ] */
      cod_trireg_solve(&(R11.matrix), lambda, &(c1.vector), S, &(y1.vector), work);

      /* save y1 for later residual calculation */
      gsl_vector_memcpy(work, &(y1.vector));

      /* compute Z [ y1; 0 ] */
      cod_householder_Zvec(QRZT, tau_Z, rank, x);

      /* compute x = P Z^T ( y1; 0 ) */
      gsl_permute_vector_inverse(perm, x);

      /* compute residual = b - A x = Q (Q^T b - [ R11 y1; 0 ]) = Q [ c1 - R11*y1 ; c2 ] */

      /* work = R11*y1 */
      gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, &(R11.matrix), work);

      gsl_vector_sub(&(c1.vector), work);
      gsl_linalg_QR_Qvec(QRZT, tau_Q, residual);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_COD_unpack()
  Unpack encoded COD decomposition into the matrices Q,R,Z,P

Inputs: QRZT  - encoded COD decomposition
        tau_Q - Householder scalars for Q
        tau_Z - Householder scalars for Z
        rank  - rank of matrix (as determined from gsl_linalg_COD_decomp)
        Q     - (output) M-by-M matrix Q
        R     - (output) M-by-N matrix R
        Z     - (output) N-by-N matrix Z
*/

int
gsl_linalg_COD_unpack(const gsl_matrix * QRZT, const gsl_vector * tau_Q,
                      const gsl_vector * tau_Z, const size_t rank, gsl_matrix * Q,
                      gsl_matrix * R, gsl_matrix * Z)
{
  const size_t M = QRZT->size1;
  const size_t N = QRZT->size2;

  if (tau_Q->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau_Q must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (tau_Z->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau_Z must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (rank > GSL_MIN (M, N))
    {
      GSL_ERROR ("rank must be <= MIN(M,N)", GSL_EBADLEN);
    }
  else if (Q->size1 != M || Q->size2 != M)
    {
      GSL_ERROR ("Q must by M-by-M", GSL_EBADLEN);
    }
  else if (R->size1 != M || R->size2 != N)
    {
      GSL_ERROR ("R must by M-by-N", GSL_EBADLEN);
    }
  else if (Z->size1 != N || Z->size2 != N)
    {
      GSL_ERROR ("Z must by N-by-N", GSL_EBADLEN);
    }
  else
    {
      size_t i;
      gsl_matrix_view R11 = gsl_matrix_submatrix(R, 0, 0, rank, rank);
      gsl_matrix_const_view QRZT11 = gsl_matrix_const_submatrix(QRZT, 0, 0, rank, rank);

      /* form Q matrix */

      gsl_matrix_set_identity(Q);

      for (i = GSL_MIN (M, N); i-- > 0;)
        {
          gsl_vector_const_view h = gsl_matrix_const_subcolumn (QRZT, i, i, M - i);
          gsl_matrix_view m = gsl_matrix_submatrix (Q, i, i, M - i, M - i);
          gsl_vector_view work = gsl_matrix_subcolumn (R, 0, 0, M - i);
          double ti = gsl_vector_get (tau_Q, i);
          double * ptr = gsl_vector_ptr((gsl_vector *) &h.vector, 0);
          double tmp = *ptr;

          *ptr = 1.0;
          gsl_linalg_householder_left (ti, &h.vector, &m.matrix, &work.vector);
          *ptr = tmp;
        }

      /* form Z matrix */
      gsl_matrix_set_identity(Z);

      if (rank < N)
        {
          gsl_vector_view work = gsl_matrix_row(R, 0); /* temporary workspace, size N */

          /* multiply I by Z from the right */
          gsl_linalg_COD_matZ(QRZT, tau_Z, rank, Z, &work.vector);
        }

      /* copy rank-by-rank upper triangle of QRZT into R and zero the rest */
      gsl_matrix_set_zero(R);
      gsl_matrix_tricpy(CblasUpper, CblasNonUnit, &R11.matrix, &QRZT11.matrix);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_COD_matZ
  Multiply an M-by-N matrix A on the right by Z (N-by-N)

Inputs: QRZT  - encoded COD matrix
        tau_Z - Householder scalars for Z
        rank  - matrix rank
        A     - on input, M-by-N matrix
                on output, A * Z
        work  - workspace of length M
*/

int
gsl_linalg_COD_matZ(const gsl_matrix * QRZT, const gsl_vector * tau_Z, const size_t rank,
                    gsl_matrix * A, gsl_vector * work)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau_Z->size != GSL_MIN (QRZT->size1, QRZT->size2))
    {
      GSL_ERROR("tau_Z must be GSL_MIN(M,N)", GSL_EBADLEN);
    }
  else if (QRZT->size2 != N)
    {
      GSL_ERROR("QRZT must have N columns", GSL_EBADLEN);
    }
  else if (work->size != M)
    {
      GSL_ERROR("workspace must be length M", GSL_EBADLEN);
    }
  else
    {
      /* if rank == N, then Z = I and there is nothing to do */
      if (rank < N)
        {
          size_t i;

          for (i = rank; i > 0 && i--; )
            {
              gsl_vector_const_view h = gsl_matrix_const_subrow (QRZT, i, rank, N - rank);
              gsl_matrix_view m = gsl_matrix_submatrix (A, 0, i, M, N - i);
              double ti = gsl_vector_get (tau_Z, i);
              cod_householder_mh (ti, &h.vector, &m.matrix, work);
            }
        }

      return GSL_SUCCESS;
    }
}


/*********************************************
 * INTERNAL ROUTINES                         *
 *********************************************/

/*
cod_RZ()
  Perform RZ decomposition of an upper trapezoidal matrix,

A = [ A11 A12 ] = [ R 0 ] Z

where A is M-by-N with N >= M, A11 is M-by-M upper triangular,
and A12 is M-by-(N-M). On output, Z is stored as Householder
reflectors in the A12 portion of A,

Z = Z(1) Z(2) ... Z(M)

Inputs: A   - M-by-N matrix with N >= M
              On input, upper trapezoidal matrix [ A11 A12 ]
              On output, A11 is overwritten by R (subdiagonal elements
              are not touched), and A12 is overwritten by Z in packed storage
        tau - (output) Householder scalars, size M
*/

static int
cod_RZ(gsl_matrix * A, gsl_vector * tau)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau->size != M)
    {
      GSL_ERROR("tau has wrong size", GSL_EBADLEN);
    }
  else if (N < M)
    {
      GSL_ERROR("N must be >= M", GSL_EINVAL);
    }
  else if (M == N)
    {
      /* quick return */
      gsl_vector_set_all(tau, 0.0);
      return GSL_SUCCESS;
    }
  else
    {
      size_t k;

      for (k = M; k > 0 && k--; )
        {
          double *alpha = gsl_matrix_ptr(A, k, k);
          gsl_vector_view z = gsl_matrix_subrow(A, k, M, N - M);
          double tauk;

          /* compute Householder reflection to zero [ A(k,k) A(k,M+1:N) ] */
          tauk = cod_householder_transform(alpha, &z.vector);
          gsl_vector_set(tau, k, tauk);

          if ((tauk != 0) && (k > 0))
            {
              gsl_vector_view w = gsl_vector_subvector(tau, 0, k);
              gsl_matrix_view B = gsl_matrix_submatrix(A, 0, k, k, N - k);

              cod_householder_mh(tauk, &z.vector, &B.matrix, &w.vector);
            }
        }

      return GSL_SUCCESS;
    }
}

static double
cod_householder_transform(double *alpha, gsl_vector * v)
{
  double beta, tau;
  double xnorm = gsl_blas_dnrm2(v);

  if (xnorm == 0)
    {
      return 0.0; /* tau = 0 */
    }

  beta = - (*alpha >= 0.0 ? +1.0 : -1.0) * gsl_hypot(*alpha, xnorm);
  tau = (beta - *alpha) / beta;

  {
    double s = (*alpha - beta);
    
    if (fabs(s) > GSL_DBL_MIN) 
      {
        gsl_blas_dscal (1.0 / s, v);
      }
    else
      {
        gsl_blas_dscal (GSL_DBL_EPSILON / s, v);
        gsl_blas_dscal (1.0 / GSL_DBL_EPSILON, v);
      }

    *alpha = beta;
  }
  
  return tau;
}

/*
cod_householder_hv
  Apply Householder reflection H = (I - tau*v*v') to vector v from the left,

w' = H * w

Inputs: tau  - Householder scalar
        v    - Householder vector, size M
        w    - on input, w vector, size M
               on output, H * w

Notes:
1) Based on LAPACK routine DLARZ
*/

static int
cod_householder_hv(const double tau, const gsl_vector * v, gsl_vector * w)
{
  if (tau == 0)
    {
      return GSL_SUCCESS; /* H = I */
    }
  else
    {
      const size_t M = w->size;
      const size_t L = v->size;
      double w0 = gsl_vector_get(w, 0);
      gsl_vector_view w1 = gsl_vector_subvector(w, M - L, L);
      double d1, d;

      /* d1 := v . w(M-L:M) */
      gsl_blas_ddot(v, &w1.vector, &d1);

      /* d := w(1) + v . w(M-L:M) */
      d = w0 + d1;

      /* w(1) = w(1) - tau * d */
      gsl_vector_set(w, 0, w0 - tau * d);

      /* w(M-L:M) = w(M-L:M) - tau * d * v */
      gsl_blas_daxpy(-tau * d, v, &w1.vector);

      return GSL_SUCCESS;
    }
}

/*
cod_householder_mh
  Apply Householder reflection H = (I - tau*v*v') to matrix A from the right

Inputs: tau  - Householder scalar
        v    - Householder vector, size N-M
        A    - matrix, size M-by-N
        work - workspace, size M

Notes:
1) Based on LAPACK routine DLARZ
*/

static int
cod_householder_mh(const double tau, const gsl_vector * v, gsl_matrix * A,
                   gsl_vector * work)
{
  if (tau == 0)
    {
      return GSL_SUCCESS; /* H = I */
    }
  else
    {
      const size_t M = A->size1;
      const size_t N = A->size2;
      const size_t L = v->size;
      gsl_vector_view A1 = gsl_matrix_subcolumn(A, 0, 0, M);
      gsl_matrix_view C = gsl_matrix_submatrix(A, 0, N - L, M, L);

      /* work(1:M) = A(1:M,1) */
      gsl_vector_memcpy(work, &A1.vector);

      /* work(1:M) = work(1:M) + A(1:M,M+1:N) * v(1:N-M) */
      gsl_blas_dgemv(CblasNoTrans, 1.0, &C.matrix, v, 1.0, work);

      /* A(1:M,1) = A(1:M,1) - tau * work(1:M) */
      gsl_blas_daxpy(-tau, work, &A1.vector);

      /* A(1:M,M+1:N) = A(1:M,M+1:N) - tau * work(1:M) * v(1:N-M)' */
      gsl_blas_dger(-tau, work, v, &C.matrix);

      return GSL_SUCCESS;
    }
}

/*
cod_householder_Zvec
  Multiply a vector by Z

Inputs: QRZT  - encoded COD matrix
        tau_Z - Householder scalars for Z
        rank  - matrix rank
        v     - on input, vector of length N
                on output, Z * v
*/

static int
cod_householder_Zvec(const gsl_matrix * QRZT, const gsl_vector * tau_Z, const size_t rank,
                     gsl_vector * v)
{
  const size_t M = QRZT->size1;
  const size_t N = QRZT->size2;

  if (tau_Z->size != GSL_MIN (M, N))
    {
      GSL_ERROR("tau_Z must be GSL_MIN(M,N)", GSL_EBADLEN);
    }
  else if (v->size != N)
    {
      GSL_ERROR("v must be length N", GSL_EBADLEN);
    }
  else
    {
      if (rank < N)
        {
          size_t i;

          for (i = 0; i < rank; ++i)
            {
              gsl_vector_const_view h = gsl_matrix_const_subrow (QRZT, i, rank, N - rank);
              gsl_vector_view w = gsl_vector_subvector (v, i, N - i);
              double ti = gsl_vector_get (tau_Z, i);
              cod_householder_hv(ti, &h.vector, &w.vector);
            }
        }

      return GSL_SUCCESS;
    }
}

/*
cod_trireg_solve()

  This function computes the solution to the least squares system

  [    R     ] x = [ b ]
  [ lambda*I ]     [ 0 ]

where R is an N-by-N upper triangular matrix, lambda is a scalar parameter,
and b is a vector of length N. This is done by computing the QR factorization

[    R     ] = W S^T
[ lambda*I ]

where S^T is upper triangular, and solving

S^T x = W^T [ b ]
            [ 0 ]

Inputs: R      - full rank upper triangular matrix; the diagonal
                 elements are modified but restored on output
        lambda - scalar parameter lambda
        b      - right hand side vector b
        S      - workspace, N-by-N
        x      - (output) least squares solution of the system
        work   - workspace of length N
*/

static int
cod_trireg_solve (const gsl_matrix * R, const double lambda, const gsl_vector * b,
                  gsl_matrix * S, gsl_vector * x, gsl_vector * work)
{
  const size_t N = R->size2;
  gsl_vector_const_view diag = gsl_matrix_const_diagonal(R);
  size_t i, j, k;

  if (lambda <= 0.0)
    {
      GSL_ERROR("lambda must be positive", GSL_EINVAL);
    }

  /* copy R and b to preserve input and initialise S; store diag(R) in work */
  gsl_matrix_transpose_tricpy(CblasUpper, CblasUnit, S, R);
  gsl_vector_memcpy(work, &diag.vector);
  gsl_vector_memcpy(x, b);

  /* eliminate the diagonal matrix lambda*I using Givens rotations */

  for (j = 0; j < N; j++)
    {
      double bj = 0.0;

      gsl_matrix_set (S, j, j, lambda);

      for (k = j + 1; k < N; k++)
        {
          gsl_matrix_set (S, k, k, 0.0);
        }

      /* the transformations to eliminate the row of lambda*I modify only a
         single element of b beyond the first n, which is initially
         zero */

      for (k = j; k < N; k++)
        {
          /* determine a Givens rotation which eliminates the
             appropriate element in the current row of lambda*I */

          double sine, cosine;

          double xk = gsl_vector_get (x, k);
          double rkk = gsl_vector_get (work, k);
          double skk = gsl_matrix_get (S, k, k);

          if (skk == 0)
            {
              continue;
            }

          if (fabs (rkk) < fabs (skk))
            {
              double cotangent = rkk / skk;
              sine = 0.5 / sqrt (0.25 + 0.25 * cotangent * cotangent);
              cosine = sine * cotangent;
            }
          else
            {
              double tangent = skk / rkk;
              cosine = 0.5 / sqrt (0.25 + 0.25 * tangent * tangent);
              sine = cosine * tangent;
            }

          /* Compute the modified diagonal element of r and the
             modified element of [b,0] */

          {
            double new_rkk = cosine * rkk + sine * skk;
            double new_xk = cosine * xk + sine * bj;
            
            bj = -sine * xk + cosine * bj;

            gsl_vector_set(work, k, new_rkk);
            gsl_matrix_set(S, k, k, new_rkk);
            gsl_vector_set(x, k, new_xk);
          }

          /* Accumulate the transformation in the row of s */

          for (i = k + 1; i < N; i++)
            {
              double sik = gsl_matrix_get (S, i, k);
              double sii = gsl_matrix_get (S, i, i);
              
              double new_sik = cosine * sik + sine * sii;
              double new_sii = -sine * sik + cosine * sii;

              gsl_matrix_set(S, i, k, new_sik);
              gsl_matrix_set(S, i, i, new_sii);
            }
        }
    }

  /* solve: S^T x = rhs in place */
  gsl_blas_dtrsv(CblasLower, CblasTrans, CblasNonUnit, S, x);

  return GSL_SUCCESS;
}
