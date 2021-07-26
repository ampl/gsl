/* linalg/test_lu_band.c
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>

static int
test_LU_band_decomp_eps(const size_t p, const size_t q, const gsl_matrix * A, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  const size_t minMN = GSL_MIN(M, N);
  size_t i, j;

  gsl_matrix * AB = gsl_matrix_alloc(N, 2*p + q + 1);
  gsl_vector_uint * piv = gsl_vector_uint_alloc(minMN);
  gsl_matrix * L = gsl_matrix_alloc(M, minMN);
  gsl_matrix * U = gsl_matrix_alloc(minMN, N);
  gsl_matrix * B = gsl_matrix_alloc(M, N);

  /* convert A to packed banded format */
  gen2band_matrix(p, q, A, AB);

  s += gsl_linalg_LU_band_decomp(M, p, q, AB, piv);
  s += gsl_linalg_LU_band_unpack(M, p, q, AB, piv, L, U);

  /* compute B = L U */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, U, 0.0, B);

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double aij = gsl_matrix_get(A, i, j);
          double bij = gsl_matrix_get(B, i, j);

          gsl_test_rel(bij, aij, eps, "%s (M=%lu,N=%lu)(p=%lu,q=%lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, p, q, i, j, aij, bij);
        }
    }

  gsl_matrix_free(AB);
  gsl_vector_uint_free(piv);
  gsl_matrix_free(L);
  gsl_matrix_free(U);
  gsl_matrix_free(B);

  return s;
}

static int
test_LU_band_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 20;
  gsl_matrix * A = gsl_matrix_alloc(N_max, N_max);
  size_t M, N, p, q;

  for (M = 1; M <= N_max; ++M)
    {
      for (N = 1; N <= N_max; ++N)
        {
          gsl_matrix_view m = gsl_matrix_submatrix(A, 0, 0, M, N);

          for (p = 0; p < GSL_MIN(M, 10); ++p)
            {
              for (q = 0; q < GSL_MIN(N, 10); ++q)
                {
                  create_band_matrix(p, q, &m.matrix, r);
                  s += test_LU_band_decomp_eps(p, q, &m.matrix, 1.0e5 * GSL_MAX(M,N) * GSL_DBL_EPSILON, "LU_band_decomp random");
                }
            }
        }
    }

  gsl_matrix_free(A);

  return s;
}

static int
test_LU_band_solve_eps(const size_t p, const size_t q, const gsl_matrix * A,
                       const gsl_vector * rhs, const gsl_vector * sol,
                       const double eps, const char * desc)
{
  int s = 0;
  const size_t N = A->size2;
  size_t i;

  gsl_matrix * AB = gsl_matrix_alloc(N, 2*p + q + 1);
  gsl_vector_uint * piv = gsl_vector_uint_alloc(N);
  gsl_matrix * L = gsl_matrix_alloc(N, N);
  gsl_matrix * U = gsl_matrix_alloc(N, N);
  gsl_vector * x = gsl_vector_alloc(N);

  /* convert A to packed banded format */
  gen2band_matrix(p, q, A, AB);

  s += gsl_linalg_LU_band_decomp(N, p, q, AB, piv);
  s += gsl_linalg_LU_band_solve(p, q, AB, piv, rhs, x);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps, "%s: %3lu[%lu]: %22.18g    %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_matrix_free(AB);
  gsl_vector_uint_free(piv);
  gsl_matrix_free(L);
  gsl_matrix_free(U);
  gsl_vector_free(x);

  return s;
}

static int
test_LU_band_solve(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 20;
  gsl_matrix * A = gsl_matrix_alloc(N_max, N_max);
  gsl_vector * b = gsl_vector_alloc(N_max);
  gsl_vector * x = gsl_vector_alloc(N_max);
  size_t N, p, q;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix_view m = gsl_matrix_submatrix(A, 0, 0, N, N);
      gsl_vector_view sol = gsl_vector_subvector(x, 0, N);
      gsl_vector_view rhs = gsl_vector_subvector(b, 0, N);

      for (p = 0; p < GSL_MIN(N, 10); ++p)
        {
          for (q = 0; q < GSL_MIN(N, 10); ++q)
            {
              create_band_matrix(p, q, &m.matrix, r);
              create_random_vector(&sol.vector, r);
              gsl_blas_dgemv(CblasNoTrans, 1.0, &m.matrix, &sol.vector, 0.0, &rhs.vector);

              s += test_LU_band_solve_eps(p, q, &m.matrix, &rhs.vector, &sol.vector,
                                          1.0e8 * N * GSL_DBL_EPSILON, "LU_band_solve random");
            }
        }
    }

  gsl_matrix_free(A);
  gsl_vector_free(b);
  gsl_vector_free(x);

  return s;
}
