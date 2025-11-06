/* linalg/test_lu.c
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
test_LU_decomp_eps(const gsl_matrix * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;

  gsl_matrix * A = gsl_matrix_alloc(M, N);
  gsl_matrix * PLU = gsl_matrix_calloc(M, N);
  gsl_permutation * p = gsl_permutation_alloc(M);
  size_t i, j;
  int signum;

  gsl_matrix_memcpy(A, m);
  gsl_linalg_LU_decomp(A, p, &signum);

  if (M >= N)
    {
      gsl_matrix_view U = gsl_matrix_submatrix(A, 0, 0, N, N);

      /* copy L factor into PLU */
      for (j = 0; j < GSL_MIN(M, N); ++j)
        {
          gsl_vector_view v1 = gsl_matrix_subcolumn(A, j, j, M - j);
          gsl_vector_view v2 = gsl_matrix_subcolumn(PLU, j, j, M - j);

          gsl_vector_memcpy(&v2.vector, &v1.vector);
          gsl_matrix_set(PLU, j, j, 1.0);
        }

      /* PLU := L * U */
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &U.matrix, PLU);
    }
  else
    {
      gsl_matrix_view L = gsl_matrix_submatrix(A, 0, 0, M, M);

      /* copy U factor into PLU */
      for (i = 0; i < GSL_MIN(M, N); ++i)
        {
          gsl_vector_view v1 = gsl_matrix_subrow(A, i, i, N - i);
          gsl_vector_view v2 = gsl_matrix_subrow(PLU, i, i, N - i);
          gsl_vector_memcpy(&v2.vector, &v1.vector);
        }

      /* PLU := L * U */
      gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 1.0, &L.matrix, PLU);
    }

  /* PLU := P * L * U */
  for (j = 0; j < N; ++j)
    {
      gsl_vector_view v = gsl_matrix_column(PLU, j);
      gsl_permute_vector_inverse(p, &v.vector);
    }

  /* now test m = PLU */
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double aij = gsl_matrix_get(PLU, i, j);
          double mij = gsl_matrix_get(m, i, j);

          if (fabs(mij) > GSL_DBL_MIN)
            {
              gsl_test_rel(aij, mij, eps,
                           "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, aij, mij);
            }
          else
            {
              gsl_test_rel(aij / GSL_DBL_MIN, mij / GSL_DBL_MIN, eps,
                           "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, aij, mij);
            }
        }
    }

  gsl_matrix_free(A);
  gsl_matrix_free(PLU);
  gsl_permutation_free(p);

  return s;
}

static int
test_LU_decomp(gsl_rng * r)
{
  int s = 0;
  size_t n;

  /* test square matrices */
  for (n = 1; n <= 50; ++n)
    {
      gsl_matrix * m = gsl_matrix_alloc(n, n);

      create_random_matrix(m, r);
      test_LU_decomp_eps(m, 4096.0 * n * GSL_DBL_EPSILON, "LU_decomp random");

      create_hilbert_matrix2(m);
      test_LU_decomp_eps(m, 256.0 * n * GSL_DBL_EPSILON, "LU_decomp hilbert");

      gsl_matrix_free(m);
    }

  /* bug #61094 */
  {
    double m_data[] = { 0.0, 0.0, 1.0,
                        0.0, 1.0, 1.0,
                        0.0, 1.0, 2.0 };
    gsl_matrix_view m = gsl_matrix_view_array(m_data, 3, 3);

    test_LU_decomp_eps(&m.matrix, 10.0 * GSL_DBL_EPSILON, "LU_decomp bug #61094");
  }

  {
    double m_data[] = { 0.98*GSL_DBL_MIN, 0.2, 1.0,
                        0.97*GSL_DBL_MIN, 1.0, 1.0,
                        0.96*GSL_DBL_MIN, 1.0, 2.0 };
    gsl_matrix_view m = gsl_matrix_view_array(m_data, 3, 3);

    test_LU_decomp_eps(&m.matrix, 10.0 * GSL_DBL_EPSILON, "LU_decomp GSL_DBL_MIN");
  }

  {
    gsl_matrix * m = gsl_matrix_alloc(100, 50);
    create_random_matrix(m, r);
    test_LU_decomp_eps(m, 256.0 * n * GSL_DBL_EPSILON, "LU_decomp rect1");
    gsl_matrix_free(m);
  }

  {
    gsl_matrix * m = gsl_matrix_alloc(50, 100);
    create_random_matrix(m, r);
    test_LU_decomp_eps(m, 1.0e3 * n * GSL_DBL_EPSILON, "LU_decomp rect2");
    gsl_matrix_free(m);
  }

  {
    gsl_matrix * m = gsl_matrix_alloc(80, 100);
    create_random_matrix(m, r);
    test_LU_decomp_eps(m, 1.0e4 * n * GSL_DBL_EPSILON, "LU_decomp rect3");
    gsl_matrix_free(m);
  }

  return s;
}

static int
test_LU_solve_eps(const gsl_matrix * m, const gsl_vector * rhs, const gsl_vector * sol, const double eps, const char * desc)
{
  int s = 0;
  const size_t N = m->size1;
  int signum;
  size_t i;

  gsl_permutation * perm = gsl_permutation_alloc(N);
  gsl_matrix * lu  = gsl_matrix_alloc(N, N);
  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * residual = gsl_vector_alloc(N);

  gsl_matrix_memcpy(lu, m);

  s += gsl_linalg_LU_decomp(lu, perm, &signum);
  s += gsl_linalg_LU_solve(lu, perm, rhs, x);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps, "%s: %3lu[%lu]: %22.18g    %22.18g\n",
                   desc, N, i, xi, yi);
    }

  s += gsl_linalg_LU_refine(m, lu, perm, rhs, x, residual);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps, "%s: improved %3lu[%lu]: %22.18g    %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_vector_free(residual);
  gsl_vector_free(x);
  gsl_matrix_free(lu);
  gsl_permutation_free(perm);

  return s;
}


static int
test_LU_solve(gsl_rng * r)
{
  int s = 0;
  size_t n;

  for (n = 1; n <= 50; ++n)
    {
      gsl_matrix * m = gsl_matrix_alloc(n, n);
      gsl_vector * rhs = gsl_vector_alloc(n);
      gsl_vector * sol = gsl_vector_alloc(n);

      create_random_matrix(m, r);
      create_random_vector(sol, r);
      gsl_blas_dgemv(CblasNoTrans, 1.0, m, sol, 0.0, rhs);
      test_LU_solve_eps(m, rhs, sol, 1.0e5 * n * GSL_DBL_EPSILON, "LU_solve random");

      if (n <= 4)
        {
          create_hilbert_matrix2(m);
          gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);
          test_LU_solve_eps(m, rhs, sol, 4096.0 * n * GSL_DBL_EPSILON, "LU_solve hilbert");
        }

      gsl_matrix_free(m);
      gsl_vector_free(rhs);
      gsl_vector_free(sol);
    }

  return s;
}

int
test_LU_invert_eps(const gsl_matrix * m, const double eps, const char *desc)
{
  int s = 0;
  size_t i, j, N = m->size1;

  gsl_matrix * lu  = gsl_matrix_alloc(N, N);
  gsl_matrix * c  = gsl_matrix_alloc(N, N);
  gsl_permutation * perm = gsl_permutation_alloc(N);
  int signum;

  gsl_matrix_memcpy(lu, m);

  s += gsl_linalg_LU_decomp(lu, perm, &signum);
  s += gsl_linalg_LU_invx(lu, perm);

  /* c = m m^{-1} */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m, lu, 0.0, c);

  /* c should be the identity matrix */

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double cij = gsl_matrix_get(c, i, j);
          double expected = (i == j) ? 1.0 : 0.0;

          gsl_test_rel(cij, expected, eps, "%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, cij, expected);
        }
    }

  gsl_matrix_free(lu);
  gsl_matrix_free(c);
  gsl_permutation_free(perm);

  return s;
}

static int
test_LU_invert(gsl_rng * r)
{
  int s = 0;
  size_t n;

  for (n = 1; n <= 50; ++n)
    {
      gsl_matrix * m = gsl_matrix_alloc(n, n);

      create_random_matrix(m, r);
      test_LU_invert_eps(m, 1.0e5 * n * GSL_DBL_EPSILON, "LU_invert random");

      gsl_matrix_free(m);
    }

  return s;
}
