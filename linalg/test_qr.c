/* linalg/test_qr.c
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>

static int
test_QR_decomp_dim(const gsl_matrix * m, double eps)
{
  int s = 0;
  unsigned long i,j, M = m->size1, N = m->size2;

  gsl_matrix * qr = gsl_matrix_alloc(M,N);
  gsl_matrix * a  = gsl_matrix_alloc(M,N);
  gsl_matrix * q  = gsl_matrix_alloc(M,M);
  gsl_matrix * r  = gsl_matrix_alloc(M,N);
  gsl_vector * d = gsl_vector_alloc(GSL_MIN(M,N));

  gsl_matrix_memcpy(qr,m);

  s += gsl_linalg_QR_decomp(qr, d);
  s += gsl_linalg_QR_unpack(qr, d, q, r);
  
  /* compute a = q r */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, q, r, 0.0, a);

  for(i=0; i<M; i++) {
    for(j=0; j<N; j++) {
      double aij = gsl_matrix_get(a, i, j);
      double mij = gsl_matrix_get(m, i, j);
      int foo = check(aij, mij, eps);
      if(foo) {
        printf("(%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n", M, N, i,j, aij, mij);
      }
      s += foo;
    }
  }

  gsl_vector_free(d);
  gsl_matrix_free(qr);
  gsl_matrix_free(a);
  gsl_matrix_free(q);
  gsl_matrix_free(r);

  return s;
}

static int
test_QR_decomp(void)
{
  int f;
  int s = 0;

  f = test_QR_decomp_dim(m35, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp m(3,5)");
  s += f;

  f = test_QR_decomp_dim(m53, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp m(5,3)");
  s += f;

  f = test_QR_decomp_dim(hilb2, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp hilbert(2)");
  s += f;

  f = test_QR_decomp_dim(hilb3, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp hilbert(3)");
  s += f;

  f = test_QR_decomp_dim(hilb4, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp hilbert(4)");
  s += f;

  f = test_QR_decomp_dim(hilb12, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp hilbert(12)");
  s += f;

  f = test_QR_decomp_dim(vander2, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp vander(2)");
  s += f;

  f = test_QR_decomp_dim(vander3, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp vander(3)");
  s += f;

  f = test_QR_decomp_dim(vander4, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  QR_decomp vander(4)");
  s += f;

  f = test_QR_decomp_dim(vander12, 0.0005); /* FIXME: bad accuracy */
  gsl_test(f, "  QR_decomp vander(12)");
  s += f;

  return s;
}

static int
test_QR_decomp_r_eps(const gsl_matrix * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  gsl_matrix * QR = gsl_matrix_alloc(M, N);
  gsl_matrix * T = gsl_matrix_alloc(N, N);
  gsl_matrix * A  = gsl_matrix_alloc(M, N);
  gsl_matrix * R  = gsl_matrix_alloc(N, N);
  gsl_matrix * Q  = gsl_matrix_alloc(M, M);
  gsl_matrix_view Q1 = gsl_matrix_submatrix(Q, 0, 0, M, N);

  gsl_matrix_memcpy(QR, m);

  s += gsl_linalg_QR_decomp_r(QR, T);
  s += gsl_linalg_QR_unpack_r(QR, T, Q, R);
  
  /* compute A = Q R */
  gsl_matrix_memcpy(A, &Q1.matrix);
  gsl_blas_dtrmm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, R, A);

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          double aij = gsl_matrix_get(A, i, j);
          double mij = gsl_matrix_get(m, i, j);

          gsl_test_rel(aij, mij, eps, "%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i,j, aij, mij);
        }
    }

  if (M > N)
    {
      gsl_matrix * R_alt  = gsl_matrix_alloc(M, N);
      gsl_matrix * Q_alt  = gsl_matrix_alloc(M, M);
      gsl_vector_view tau = gsl_matrix_diagonal(T);

      /* test that Q2 was computed correctly by comparing with Level 2 algorithm */
      gsl_linalg_QR_unpack(QR, &tau.vector, Q_alt, R_alt);

      for (i = 0; i < M; i++)
        {
          for (j = 0; j < M; j++)
            {
              double aij = gsl_matrix_get(Q, i, j);
              double bij = gsl_matrix_get(Q_alt, i, j);

              gsl_test_rel(aij, bij, eps, "%s Q (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, aij, bij);
            }
        }

      gsl_matrix_free(R_alt);
      gsl_matrix_free(Q_alt);
    }

  gsl_matrix_free(QR);
  gsl_matrix_free(T);
  gsl_matrix_free(A);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);

  return s;
}

static int
test_QR_decomp_r(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 50; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix * A = gsl_matrix_alloc(M, N);

          create_random_matrix(A, r);
          s += test_QR_decomp_r_eps(A, 1.0e5 * M * GSL_DBL_EPSILON, "QR_decomp_r random");

          gsl_matrix_free(A);
        }
    }

  s += test_QR_decomp_r_eps(m53,   1.0e2 * GSL_DBL_EPSILON, "QR_decomp_r m(5,3)");
  s += test_QR_decomp_r_eps(hilb2, 1.0e2 * GSL_DBL_EPSILON, "QR_decomp_r hilbert(2)");
  s += test_QR_decomp_r_eps(hilb3, 1.0e2 * GSL_DBL_EPSILON, "QR_decomp_r hilbert(3)");
  s += test_QR_decomp_r_eps(hilb4, 1.0e2 * GSL_DBL_EPSILON, "QR_decomp_r hilbert(4)");
  s += test_QR_decomp_r_eps(hilb12, 1.0e2 * GSL_DBL_EPSILON, "QR_decomp_r hilbert(12)");
  s += test_QR_decomp_r_eps(vander2, 1.0e1 * GSL_DBL_EPSILON, "QR_decomp_r vander(2)");
  s += test_QR_decomp_r_eps(vander3, 1.0e1 * GSL_DBL_EPSILON, "QR_decomp_r vander(3)");
  s += test_QR_decomp_r_eps(vander4, 1.0e1 * GSL_DBL_EPSILON, "QR_decomp_r vander(4)");

  return s;
}

static int
test_QR_QTmat_r_eps(const gsl_matrix * A, const gsl_matrix * B, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  const size_t K = B->size2;
  size_t i, j;

  gsl_matrix * QR = gsl_matrix_alloc(M, N);
  gsl_matrix * T = gsl_matrix_alloc(N, N);
  gsl_matrix * work  = gsl_matrix_alloc(N, K);
  gsl_matrix * B1 = gsl_matrix_alloc(M, K);
  gsl_matrix * B2 = gsl_matrix_alloc(M, K);
  gsl_vector_view tau = gsl_matrix_diagonal(T);

  gsl_matrix_memcpy(QR, A);
  gsl_matrix_memcpy(B1, B);
  gsl_matrix_memcpy(B2, B);

  s += gsl_linalg_QR_decomp_r(QR, T);

  /* compute Q^T B with both recursive and non-recursive methods and compare */
  s += gsl_linalg_QR_QTmat_r(QR, T, B1, work);
  s += gsl_linalg_QR_QTmat(QR, &tau.vector, B2);
  
  for (i = 0; i < M; i++)
    {
      for (j = 0; j < K; j++)
        {
          double aij = gsl_matrix_get(B1, i, j);
          double bij = gsl_matrix_get(B2, i, j);

          gsl_test_rel(aij, bij, eps, "%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, K, i,j, aij, bij);
        }
    }

  gsl_matrix_free(QR);
  gsl_matrix_free(T);
  gsl_matrix_free(B1);
  gsl_matrix_free(B2);
  gsl_matrix_free(work);

  return s;
}

static int
test_QR_QTmat_r(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 50; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          const size_t K = GSL_MAX(N / 2, 1);
          gsl_matrix * A = gsl_matrix_alloc(M, N);
          gsl_matrix * B = gsl_matrix_alloc(M, K);

          create_random_matrix(A, r);
          create_random_matrix(B, r);
          s += test_QR_QTmat_r_eps(A, B, 1.0e6 * M * GSL_DBL_EPSILON, "QR_QTmat_r random");

          gsl_matrix_free(A);
          gsl_matrix_free(B);
        }
    }

  return s;
}

static int
test_QR_solve_r_eps(const gsl_matrix * A, const gsl_vector * rhs, const gsl_vector * sol,
                    const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i;

  gsl_matrix * QR = gsl_matrix_alloc(M, N);
  gsl_matrix * T = gsl_matrix_alloc(N, N);
  gsl_vector * x  = gsl_vector_alloc(N);

  gsl_matrix_memcpy(QR, A);

  s += gsl_linalg_QR_decomp_r(QR, T);
  s += gsl_linalg_QR_solve_r(QR, T, rhs, x);

  for (i = 0; i < M; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps, "%s (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_matrix_free(QR);
  gsl_matrix_free(T);
  gsl_vector_free(x);

  return s;
}

static int
test_QR_solve_r(gsl_rng * r)
{
  int s = 0;
  size_t N;

  for (N = 1; N <= 50; ++N)
    {
      gsl_matrix * A = gsl_matrix_alloc(N, N);
      gsl_vector * sol = gsl_vector_alloc(N);
      gsl_vector * rhs = gsl_vector_alloc(N);

      create_random_matrix(A, r);
      create_random_vector(sol, r);
      gsl_blas_dgemv(CblasNoTrans, 1.0, A, sol, 0.0, rhs);

      s += test_QR_solve_r_eps(A, rhs, sol, 1.0e5 * N * GSL_DBL_EPSILON, "QR_solve_r random");

      gsl_matrix_free(A);
      gsl_vector_free(sol);
      gsl_vector_free(rhs);
    }

  return s;
}

static int
test_QR_lssolve_r_eps(const gsl_matrix * A, const gsl_vector * b, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i;

  gsl_matrix * QR = gsl_matrix_alloc(M, N);
  gsl_matrix * T = gsl_matrix_alloc(N, N);
  gsl_vector * x = gsl_vector_alloc(M);
  gsl_vector * work  = gsl_vector_alloc(N);

  gsl_matrix * U = gsl_matrix_alloc(M, N);
  gsl_matrix * V = gsl_matrix_alloc(N, N);
  gsl_vector * S = gsl_vector_alloc(N);
  gsl_vector * x_svd = gsl_vector_alloc(N);

  gsl_vector * residual = gsl_vector_alloc(M);

  gsl_matrix_memcpy(QR, A);
  s += gsl_linalg_QR_decomp_r(QR, T);
  s += gsl_linalg_QR_lssolve_r(QR, T, b, x, work);

  gsl_matrix_memcpy(U, A);
  gsl_linalg_SV_decomp(U, V, S, work);
  gsl_linalg_SV_solve(U, V, S, b, x_svd);

  /* compare QR with SVD solution */
  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(x_svd, i);

      gsl_test_rel(xi, yi, eps, "%s (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, M, N, i, xi, yi);
    }

  if (M > N)
    {
      gsl_vector_view x1 = gsl_vector_subvector(x, 0, N);
      gsl_vector_view x2 = gsl_vector_subvector(x, N, M - N);
      double norm = gsl_blas_dnrm2(&x2.vector);
      double norm_expected;

      /* compute residual and check || x(N+1:end) || = || b - A x || */
      gsl_vector_memcpy(residual, b);
      gsl_blas_dgemv(CblasNoTrans, -1.0, A, &x1.vector, 1.0, residual);
      norm_expected = gsl_blas_dnrm2(residual);

      gsl_test_rel(norm, norm_expected, eps, "%s rnorm (%3lu,%3lu): %22.18g   %22.18g\n",
                   desc, M, N, norm, norm_expected);
    }

  gsl_matrix_free(QR);
  gsl_matrix_free(T);
  gsl_vector_free(x);
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(x_svd);
  gsl_vector_free(work);
  gsl_vector_free(S);
  gsl_vector_free(residual);

  return s;
}

static int
test_QR_lssolve_r(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 30; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix * A = gsl_matrix_alloc(M, N);
          gsl_vector * b = gsl_vector_alloc(M);

          create_random_matrix(A, r);
          create_random_vector(b, r);

          s += test_QR_lssolve_r_eps(A, b, 1.0e4 * M * GSL_DBL_EPSILON, "QR_lssolve_r random");

          gsl_matrix_free(A);
          gsl_vector_free(b);
        }
    }

  return s;
}

static int
test_QR_TR_decomp_eps(const gsl_matrix * S, const gsl_matrix * A, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i, j;

  gsl_matrix * R = gsl_matrix_calloc(N, N);
  gsl_matrix * V = gsl_matrix_alloc(M, N);
  gsl_matrix * T = gsl_matrix_alloc(N, N);
  gsl_matrix * B1  = gsl_matrix_calloc(N, N);
  gsl_matrix * B2  = gsl_matrix_alloc(M, N);

  gsl_matrix_tricpy(CblasUpper, CblasNonUnit, R, S);
  gsl_matrix_memcpy(V, A);

  s += gsl_linalg_QR_TR_decomp(R, V, T);
  
  /*
   * compute B = Q R = [ R - T R ]
   *                   [ -V~ T R ]
   */

  gsl_matrix_tricpy(CblasUpper, CblasNonUnit, B1, T);
  gsl_matrix_memcpy(B2, V);

  gsl_blas_dtrmm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, R, B1);   /* B1 = T R */
  gsl_blas_dtrmm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, B1, B2); /* B2 = -V~ T R */
  gsl_matrix_sub (B1, R);                                                            /* B1 = T R - R */
  gsl_matrix_scale (B1, -1.0);                                                       /* B1 = R - T R */

  /* test S = B1 */
  for (i = 0; i < N; i++)
    {
      for (j = i; j < N; j++)
        {
          double aij = gsl_matrix_get(B1, i, j);
          double mij = gsl_matrix_get(S, i, j);

          gsl_test_rel(aij, mij, eps, "%s B1 (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i,j, aij, mij);
        }
    }

  /* test A = B2 */
  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          double aij = gsl_matrix_get(B2, i, j);
          double mij = gsl_matrix_get(A, i, j);

          gsl_test_rel(aij, mij, eps, "%s B2 (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i,j, aij, mij);
        }
    }

  gsl_matrix_free(R);
  gsl_matrix_free(V);
  gsl_matrix_free(T);
  gsl_matrix_free(B1);
  gsl_matrix_free(B2);

  return s;
}

static int
test_QR_TR_decomp(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 50; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix * S = gsl_matrix_alloc(N, N);
          gsl_matrix * A = gsl_matrix_alloc(M, N);

          create_random_matrix(A, r);
          create_random_matrix(S, r);
          s += test_QR_TR_decomp_eps(S, A, 1.0e5 * M * GSL_DBL_EPSILON, "QR_TR_decomp random");

          gsl_matrix_free(S);
          gsl_matrix_free(A);
        }
    }

  return s;
}
