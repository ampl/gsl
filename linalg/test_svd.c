/* linalg/test_svd.c
 *
 * Copyright (C) 2024 Patrick Alken
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

static int test_SV_lssolve_eps(const double lambda, const gsl_matrix * A, const gsl_vector * b, const double eps, const char *desc);
static int test_SV_lssolve(gsl_rng * r);

/* solve: min ||b - A x||^2 + lambda^2 ||x||^2 */
static int
test_SV_lssolve_eps(const double lambda, const gsl_matrix * A, const gsl_vector * b, const double eps, const char *desc)
{
  int s = 0;
  size_t i, M = A->size1, N = A->size2;

  gsl_matrix * U = gsl_matrix_alloc(M, N);
  gsl_matrix * V = gsl_matrix_alloc(N, N);
  gsl_vector * S = gsl_vector_alloc(N);
  gsl_vector * work = gsl_vector_alloc(N);
  gsl_vector * work2 = gsl_vector_alloc(M+N);
  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * x_aug = gsl_vector_alloc(N);
  gsl_vector * r = gsl_vector_alloc(M);
  double rnorm, rnorm_expected;

  {
    /* form full rank augmented system B = [ A ; lambda*I_N ], f = [ rhs ; 0 ] and solve with QRPT */
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

  gsl_matrix_memcpy(U, A);
  s += gsl_linalg_SV_decomp(U, V, S, work);
  s += gsl_linalg_SV_lssolve(lambda, U, V, S, b, x, &rnorm, work2);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(x_aug, i);
      gsl_test_rel(xi, yi, eps,
                   "%s lambda=%g (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, lambda, M, N, i, xi, yi);
    }

  /* compute residual r = b - A x */
  gsl_vector_memcpy(r, b);
  gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);
  rnorm_expected = gsl_blas_dnrm2(r);

  gsl_test_rel(rnorm, rnorm_expected, sqrt(eps),
               "%s rnorm lambda=%g (%3lu,%3lu): %22.18g   %22.18g\n",
               desc, lambda, M, N, rnorm, rnorm_expected);

  gsl_vector_free(r);
  gsl_vector_free(x);
  gsl_vector_free(x_aug);
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);
  gsl_vector_free(work2);

  return s;
}

static int
test_SV_lssolve(gsl_rng * r)
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
  create_random_vector(b, r);

  create_rank_matrix(26, m, r);
  test_SV_lssolve_eps(lambda, m, b, 1.0e4 * M * GSL_DBL_EPSILON, "SV_lssolve rank 26");

  create_rank_matrix(N, m, r);
  test_SV_lssolve_eps(lambda, m, b, 1.0e4 * M * GSL_DBL_EPSILON, "SV_lssolve rank N");
  test_SV_lssolve_eps(0.0, m, b, 1.0e5 * M * GSL_DBL_EPSILON, "SV_lssolve rank N unreg");

  gsl_matrix_free(m);
  gsl_vector_free(b);

  M = 500;
  N = 450;
  m = gsl_matrix_alloc(M, N);
  b = gsl_vector_alloc(M);
  create_rank_matrix(278, m, r);
  create_random_vector(b, r);
  test_SV_lssolve_eps(lambda, m, b, 1.0e6 * M * GSL_DBL_EPSILON, "SV_lssolve rank 278");
  gsl_matrix_free(m);
  gsl_vector_free(b);

  return s;
}
