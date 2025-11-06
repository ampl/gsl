/* linalg/test_cod.c
 *
 * Copyright (C) 2017 Patrick Alken
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

static int test_COD_decomp_eps(const gsl_matrix * m, const double eps, const char *desc);
static int test_COD_decomp(gsl_rng *r);
static int test_COD_lssolve_eps(const gsl_matrix * m, const double * actual, const double eps, const char *desc);
static int test_COD_lssolve(void);
static int test_COD_lssolve2_eps(const double lambda, const gsl_matrix * A, const gsl_vector * b, const double eps, const char *desc);
static int test_COD_lssolve2(gsl_rng * r);

static int
test_COD_decomp_eps(const gsl_matrix * m, const double eps, const char *desc)
{
  int s = 0;
  size_t i, j, M = m->size1, N = m->size2;
  size_t rank;

  gsl_matrix * QRZT = gsl_matrix_alloc(M, N);
  gsl_matrix * Q = gsl_matrix_alloc(M, M);
  gsl_matrix * R = gsl_matrix_alloc(M, N);
  gsl_matrix * QR = gsl_matrix_alloc(M, N);
  gsl_matrix * Z = gsl_matrix_alloc(N, N);
  gsl_vector * tau_Q = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * tau_Z = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * work = gsl_vector_alloc(N);
  gsl_matrix * lhs = gsl_matrix_alloc(M, N);
  gsl_matrix * rhs = gsl_matrix_alloc(M, N);

  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_matrix_memcpy(QRZT, m);

  s += gsl_linalg_COD_decomp(QRZT, tau_Q, tau_Z, perm, &rank, work);
  s += gsl_linalg_COD_unpack(QRZT, tau_Q, tau_Z, rank, Q, R, Z);

  /* compute lhs = m P */
  gsl_matrix_memcpy(lhs, m);
  gsl_permute_matrix(perm, lhs);

  /* compute rhs = Q R Z^T */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, R, 0.0, QR);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, QR, Z, 0.0, rhs);

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          double aij = gsl_matrix_get(rhs, i, j);
          double bij = gsl_matrix_get(lhs, i, j);

          gsl_test_rel(aij, bij, eps, "%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i, j, aij, bij);
        }
    }

  gsl_permutation_free (perm);
  gsl_vector_free(work);
  gsl_vector_free(tau_Q);
  gsl_vector_free(tau_Z);
  gsl_matrix_free(QRZT);
  gsl_matrix_free(lhs);
  gsl_matrix_free(rhs);
  gsl_matrix_free(QR);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);
  gsl_matrix_free(Z);

  return s;
}

static int
test_COD_decomp(gsl_rng *r)
{
  int s = 0;
  size_t N;
  gsl_matrix *m;

  /* test COD decomposition on Hilbert matrices */
  for (N = 1; N <= 12; ++N)
    {
      m = gsl_matrix_alloc(N, N);
      create_hilbert_matrix2(m);

      test_COD_decomp_eps(m, 256.0 * N * GSL_DBL_EPSILON, "COD_decomp hilbert");

      gsl_matrix_free(m);
    }

  /* build some matrices of a given rank and test */

  m = gsl_matrix_alloc(100, 50);
  create_rank_matrix(26, m, r);
  test_COD_decomp_eps(m, 1.0e2 * GSL_DBL_EPSILON, "COD_decomp rank 26");
  gsl_matrix_free(m);

  m = gsl_matrix_alloc(550, 200);
  create_rank_matrix(176, m, r);
  test_COD_decomp_eps(m, 1.0e3 * GSL_DBL_EPSILON, "COD_decomp rank 176");
  gsl_matrix_free(m);

  test_COD_decomp_eps(m35, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp m(3,5)");
  test_COD_decomp_eps(m53, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp m(5,3)");

  test_COD_decomp_eps(s35, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp s(3,5)");
  test_COD_decomp_eps(s53, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp s(5,3)");

  test_COD_decomp_eps(vander2, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp vander(2)");
  test_COD_decomp_eps(vander3, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp vander(3)");
  test_COD_decomp_eps(vander4, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp vander(4)");
  test_COD_decomp_eps(vander12, 1e-3, "COD_decomp vander(12)"); /* FIXME: large tolerance needed */

  return s;
}

static int
test_COD_lssolve_eps(const gsl_matrix * m, const double * actual, const double eps, const char *desc)
{
  int s = 0;
  size_t i, M = m->size1, N = m->size2;

  gsl_vector * lhs = gsl_vector_alloc(M);
  gsl_vector * rhs = gsl_vector_alloc(M);
  gsl_matrix * QRZT  = gsl_matrix_alloc(M, N);
  gsl_vector * tau_Q = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * tau_Z = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * work = gsl_vector_alloc(N);
  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * r = gsl_vector_alloc(M);
  gsl_vector * res = gsl_vector_alloc(M);
  gsl_permutation * perm = gsl_permutation_alloc(N);
  size_t rank;

  gsl_matrix_memcpy(QRZT, m);

  for (i = 0; i < M; i++)
    gsl_vector_set(rhs, i, i + 1.0);

  s += gsl_linalg_COD_decomp(QRZT, tau_Q, tau_Z, perm, &rank, work);
  s += gsl_linalg_COD_lssolve(QRZT, tau_Q, tau_Z, perm, rank, rhs, x, res);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      gsl_test_rel(xi, actual[i], eps,
                   "%s (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, M, N, i, xi, actual[i]);
    }

  /* compute residual r = b - m x */
  if (M == N)
    {
      gsl_vector_set_zero(r);
    }
  else
    {
      gsl_vector_memcpy(r, rhs);
      gsl_blas_dgemv(CblasNoTrans, -1.0, m, x, 1.0, r);
    }

  for (i = 0; i < N; i++)
    {
      double r1 = gsl_vector_get(res, i);
      double r2 = gsl_vector_get(r, i);

      if (fabs(r2) < 1.0e3 * GSL_DBL_EPSILON)
        {
          gsl_test_abs(r1, r2, 10.0 * eps,
                       "%s res (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i, r1, r2);
        }
      else
        {
          gsl_test_rel(r1, r2, eps,
                       "%s res (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i, r1, r2);
        }
    }

  gsl_vector_free(r);
  gsl_vector_free(res);
  gsl_vector_free(x);
  gsl_vector_free(tau_Q);
  gsl_vector_free(tau_Z);
  gsl_matrix_free(QRZT);
  gsl_vector_free(rhs);
  gsl_vector_free(lhs);
  gsl_vector_free(work);
  gsl_permutation_free(perm);

  return s;
}

static int
test_COD_lssolve(void)
{
  int s = 0;

  test_COD_lssolve_eps(m53, m53_lssolution, 1.0e4 * GSL_DBL_EPSILON, "COD_lssolve m(5,3)");

  test_COD_lssolve_eps(hilb2, hilb2_solution, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve hilbert(2)");
  test_COD_lssolve_eps(hilb3, hilb3_solution, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve hilbert(3)");
  test_COD_lssolve_eps(hilb4, hilb4_solution, 1.0e4 * GSL_DBL_EPSILON, "COD_lssolve hilbert(4)");

  test_COD_lssolve_eps(vander2, vander2_solution, 1.0e1 * GSL_DBL_EPSILON, "COD_lssolve vander(2)");
  test_COD_lssolve_eps(vander3, vander3_solution, 1.0e1 * GSL_DBL_EPSILON, "COD_lssolve vander(3)");
  test_COD_lssolve_eps(vander4, vander4_solution, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve vander(4)");

  /* rank-1 least squares problem from the 'lin2' test dataset for multifit_nlinear */
  {
    const size_t M = 20;
    const size_t N = 5;

    /* unique minimum norm solution from "Factorize" matlab package */
    const double x_ls[] = { 1.818181818181817e-02, 3.636363636363636e-02, 5.454545454545454e-02,
                            7.272727272727272e-02, 9.090909090909088e-02 };

    gsl_matrix *m = gsl_matrix_alloc(M, N);
    size_t i, j;

    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            gsl_matrix_set(m, i, j, (i + 1.0) * (j + 1.0));
          }
      }

    test_COD_lssolve_eps(m, x_ls, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve lin2");

    gsl_matrix_free(m);
  }

  return s;
}

/* solve: min ||b - A x||^2 + lambda^2 ||x||^2 */
static int
test_COD_lssolve2_eps(const double lambda, const gsl_matrix * A, const gsl_vector * b, const double eps, const char *desc)
{
  int s = 0;
  size_t i, M = A->size1, N = A->size2;

  gsl_vector * lhs = gsl_vector_alloc(M);
  gsl_matrix * QRZT  = gsl_matrix_alloc(M, N);
  gsl_vector * tau_Q = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * tau_Z = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * work = gsl_vector_alloc(N);
  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * x_aug = gsl_vector_alloc(N);
  gsl_vector * r = gsl_vector_alloc(M);
  gsl_vector * res = gsl_vector_alloc(M);
  gsl_permutation * perm = gsl_permutation_alloc(N);
  size_t rank;

  /* form full rank augmented system B = [ A ; lambda*I_N ], f = [ rhs ; 0 ] and solve with QRPT */
  {
    gsl_vector_view v;
    gsl_matrix_view m;
    gsl_permutation *p = gsl_permutation_alloc(N);
    gsl_matrix * B = gsl_matrix_calloc(M + N, N);
    gsl_vector * f = gsl_vector_calloc(M + N);
    gsl_vector * tau = gsl_vector_alloc(N);
    gsl_vector * residual = gsl_vector_alloc(M + N);
    int signum;
    
    m = gsl_matrix_submatrix(B, 0, 0, M, N);
    gsl_matrix_memcpy(&m.matrix, A);

    m = gsl_matrix_submatrix(B, M, 0, N, N);
    v = gsl_matrix_diagonal(&m.matrix);
    gsl_vector_set_all(&v.vector, lambda);

    v = gsl_vector_subvector(f, 0, M);
    gsl_vector_memcpy(&v.vector, b);

    /* solve: [ A ; lambda*I ] x_aug = [ b ; 0 ] */
    gsl_linalg_QRPT_decomp(B, tau, p, &signum, work);
    gsl_linalg_QRPT_lssolve(B, tau, p, f, x_aug, residual);

    gsl_permutation_free(p);
    gsl_matrix_free(B);
    gsl_vector_free(f);
    gsl_vector_free(tau);
    gsl_vector_free(residual);
  }

  gsl_matrix_memcpy(QRZT, A);

  s += gsl_linalg_COD_decomp(QRZT, tau_Q, tau_Z, perm, &rank, work);

  {
    gsl_matrix *S = gsl_matrix_alloc(rank, rank);
    gsl_vector *workr = gsl_vector_alloc(rank);

    s += gsl_linalg_COD_lssolve2(lambda, QRZT, tau_Q, tau_Z, perm, rank, b, x, res, S, workr);

    gsl_matrix_free(S);
    gsl_vector_free(workr);
  }

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(x_aug, i);
      gsl_test_rel(xi, yi, eps,
                   "%s (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, M, N, i, xi, yi);
    }

  /* compute residual r = b - A x */
  if (M == N)
    {
      gsl_vector_set_zero(r);
    }
  else
    {
      gsl_vector_memcpy(r, b);
      gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);
    }

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(res, i);
      double yi = gsl_vector_get(r, i);

      gsl_test_rel(xi, yi, sqrt(eps),
                   "%s res (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, M, N, i, xi, yi);
    }

  gsl_vector_free(r);
  gsl_vector_free(res);
  gsl_vector_free(x);
  gsl_vector_free(x_aug);
  gsl_vector_free(tau_Q);
  gsl_vector_free(tau_Z);
  gsl_matrix_free(QRZT);
  gsl_vector_free(lhs);
  gsl_vector_free(work);
  gsl_permutation_free(perm);

  return s;
}

static int
test_COD_lssolve2(gsl_rng * r)
{
  int s = 0;
  gsl_matrix *m;
  gsl_vector *b;
  size_t M, N;
  const double lambda = 2.3;

  M = 100;
  N = 50;
  m = gsl_matrix_alloc(M, N);
  b = gsl_vector_alloc(M);
  create_rank_matrix(26, m, r);
  create_random_vector(b, r);
  test_COD_lssolve2_eps(lambda, m, b, 1.0e4 * M * GSL_DBL_EPSILON, "COD_lssolve2 rank 26");
  gsl_matrix_free(m);
  gsl_vector_free(b);

  M = 500;
  N = 450;
  m = gsl_matrix_alloc(M, N);
  b = gsl_vector_alloc(M);
  create_rank_matrix(278, m, r);
  create_random_vector(b, r);
  test_COD_lssolve2_eps(lambda, m, b, 1.0e6 * M * GSL_DBL_EPSILON, "COD_lssolve2 rank 278");
  gsl_matrix_free(m);
  gsl_vector_free(b);

  return s;
}
