/* linalg/test_luc.c
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
test_LUc_decomp_eps(const gsl_matrix_complex * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;

  gsl_matrix_complex * A = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * PLU = gsl_matrix_complex_calloc(M, N);
  gsl_permutation * p = gsl_permutation_alloc(M);
  size_t i, j;
  int signum;

  gsl_matrix_complex_memcpy(A, m);
  gsl_linalg_complex_LU_decomp(A, p, &signum);

  if (M >= N)
    {
      gsl_matrix_complex_view U = gsl_matrix_complex_submatrix(A, 0, 0, N, N);

      /* copy L factor into PLU */
      for (j = 0; j < GSL_MIN(M, N); ++j)
        {
          gsl_vector_complex_view v1 = gsl_matrix_complex_subcolumn(A, j, j, M - j);
          gsl_vector_complex_view v2 = gsl_matrix_complex_subcolumn(PLU, j, j, M - j);

          gsl_vector_complex_memcpy(&v2.vector, &v1.vector);
          gsl_matrix_complex_set(PLU, j, j, GSL_COMPLEX_ONE);
        }

      /* PLU := L * U */
      gsl_blas_ztrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, GSL_COMPLEX_ONE, &U.matrix, PLU);
    }
  else
    {
      gsl_matrix_complex_view L = gsl_matrix_complex_submatrix(A, 0, 0, M, M);

      /* copy U factor into PLU */
      for (i = 0; i < GSL_MIN(M, N); ++i)
        {
          gsl_vector_complex_view v1 = gsl_matrix_complex_subrow(A, i, i, N - i);
          gsl_vector_complex_view v2 = gsl_matrix_complex_subrow(PLU, i, i, N - i);
          gsl_vector_complex_memcpy(&v2.vector, &v1.vector);
        }

      /* PLU := L * U */
      gsl_blas_ztrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, GSL_COMPLEX_ONE, &L.matrix, PLU);
    }

  /* PLU := P * L * U */
  for (j = 0; j < N; ++j)
    {
      gsl_vector_complex_view v = gsl_matrix_complex_column(PLU, j);
      gsl_permute_vector_complex_inverse(p, &v.vector);
    }

  /* now test m = PLU */
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex aij = gsl_matrix_complex_get(PLU, i, j);
          gsl_complex mij = gsl_matrix_complex_get(m, i, j);

          if (GSL_REAL(mij) == 0.0 || fabs(GSL_REAL(mij)) > GSL_DBL_MIN)
            {
              gsl_test_rel(GSL_REAL(aij), GSL_REAL(mij), eps,
                           "%s real: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, GSL_REAL(aij), GSL_REAL(mij));
            }
          else
            {
              gsl_test_rel(GSL_REAL(aij) / GSL_DBL_MIN, GSL_REAL(mij) / GSL_DBL_MIN, eps,
                           "%s real: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, GSL_REAL(aij), GSL_REAL(mij));
            }

          if (GSL_IMAG(mij) == 0.0 || fabs(GSL_IMAG(mij)) > GSL_DBL_MIN)
            {
              gsl_test_rel(GSL_IMAG(aij), GSL_IMAG(mij), eps,
                           "%s imag: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, GSL_IMAG(aij), GSL_IMAG(mij));
            }
          else
            {
              gsl_test_rel(GSL_IMAG(aij) / GSL_DBL_MIN, GSL_IMAG(mij) / GSL_DBL_MIN, eps,
                           "%s imag: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, GSL_IMAG(aij), GSL_IMAG(mij));
            }
        }
    }

  gsl_matrix_complex_free(A);
  gsl_matrix_complex_free(PLU);
  gsl_permutation_free(p);

  return s;
}

static int
test_LUc_decomp(gsl_rng * r)
{
  int s = 0;
  size_t n;

  /* test square matrices */
  for (n = 1; n <= 50; ++n)
    {
      gsl_matrix_complex * m = gsl_matrix_complex_alloc(n, n);

      create_random_complex_matrix(m, r);
      test_LUc_decomp_eps(m, 1.0e4 * n * GSL_DBL_EPSILON, "complex_LU_decomp random");

      gsl_matrix_complex_free(m);
    }

  /* bug #61094 */
  {
    double m_data[] = { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 
                        0.0, 0.0, 1.0, 0.0, 2.0, 0.0 };
    gsl_matrix_complex_view m = gsl_matrix_complex_view_array(m_data, 3, 3);

    test_LUc_decomp_eps(&m.matrix, 10.0 * GSL_DBL_EPSILON, "complex_LU_decomp bug #61094");
  }

  {
    double m_data[] = { 0.98*GSL_DBL_MIN, 0.0, 0.2, 0.0, 1.0, 0.0,
                        0.97*GSL_DBL_MIN, 0.0, 1.0, 0.0, 1.0, 0.0,
                        0.96*GSL_DBL_MIN, 0.0, 1.0, 0.0, 2.0, 0.0 };
    gsl_matrix_complex_view m = gsl_matrix_complex_view_array(m_data, 3, 3);

    test_LUc_decomp_eps(&m.matrix, 10.0 * GSL_DBL_EPSILON, "complex_LU_decomp GSL_DBL_MIN");
  }

  {
    gsl_matrix_complex * m = gsl_matrix_complex_alloc(100, 50);
    create_random_complex_matrix(m, r);
    test_LUc_decomp_eps(m, 256.0 * n * GSL_DBL_EPSILON, "complex_LU_decomp rect1");
    gsl_matrix_complex_free(m);
  }

  {
    gsl_matrix_complex * m = gsl_matrix_complex_alloc(50, 100);
    create_random_complex_matrix(m, r);
    test_LUc_decomp_eps(m, 1.0e3 * n * GSL_DBL_EPSILON, "complex_LU_decomp rect2");
    gsl_matrix_complex_free(m);
  }

  {
    gsl_matrix_complex * m = gsl_matrix_complex_alloc(80, 100);
    create_random_complex_matrix(m, r);
    test_LUc_decomp_eps(m, 1.0e3 * n * GSL_DBL_EPSILON, "complex_LU_decomp rect3");
    gsl_matrix_complex_free(m);
  }

  return s;
}

static int
test_LUc_solve_eps(const gsl_matrix_complex * m, const gsl_vector_complex * rhs, const gsl_vector_complex * sol, const double eps, const char * desc)
{
  int s = 0;
  const size_t N = m->size1;
  int signum;
  size_t i;

  gsl_permutation * perm = gsl_permutation_alloc(N);
  gsl_matrix_complex * lu  = gsl_matrix_complex_alloc(N, N);
  gsl_vector_complex * x = gsl_vector_complex_alloc(N);
  gsl_vector_complex * residual = gsl_vector_complex_alloc(N);

  gsl_matrix_complex_memcpy(lu, m);

  s += gsl_linalg_complex_LU_decomp(lu, perm, &signum);
  s += gsl_linalg_complex_LU_solve(lu, perm, rhs, x);

  for (i = 0; i < N; i++)
    {
      gsl_complex xi = gsl_vector_complex_get(x, i);
      gsl_complex yi = gsl_vector_complex_get(sol, i);

      gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps, "%s: real %3lu[%lu]: %22.18g    %22.18g\n",
                   desc, N, i, GSL_REAL(xi), GSL_REAL(yi));

      gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps, "%s: imag %3lu[%lu]: %22.18g    %22.18g\n",
                   desc, N, i, GSL_IMAG(xi), GSL_IMAG(yi));
    }

  s += gsl_linalg_complex_LU_refine(m, lu, perm, rhs, x, residual);

  for (i = 0; i < N; i++)
    {
      gsl_complex xi = gsl_vector_complex_get(x, i);
      gsl_complex yi = gsl_vector_complex_get(sol, i);

      gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps, "%s: improved real %3lu[%lu]: %22.18g    %22.18g\n",
                   desc, N, i, GSL_REAL(xi), GSL_REAL(yi));

      gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps, "%s: improved imag %3lu[%lu]: %22.18g    %22.18g\n",
                   desc, N, i, GSL_IMAG(xi), GSL_IMAG(yi));
    }

  gsl_vector_complex_free(residual);
  gsl_vector_complex_free(x);
  gsl_matrix_complex_free(lu);
  gsl_permutation_free(perm);

  return s;
}

static int
test_LUc_solve(gsl_rng * r)
{
  int s = 0;
  size_t n;

  for (n = 1; n <= 50; ++n)
    {
      gsl_matrix_complex * m = gsl_matrix_complex_alloc(n, n);
      gsl_vector_complex * rhs = gsl_vector_complex_alloc(n);
      gsl_vector_complex * sol = gsl_vector_complex_alloc(n);

      create_random_complex_matrix(m, r);
      create_random_complex_vector(sol, r);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, m, sol, GSL_COMPLEX_ZERO, rhs);
      test_LUc_solve_eps(m, rhs, sol, 1.0e5 * n * GSL_DBL_EPSILON, "complex_LU_solve random");

      gsl_matrix_complex_free(m);
      gsl_vector_complex_free(rhs);
      gsl_vector_complex_free(sol);
    }

  return s;
}

int
test_LUc_invert_eps(const gsl_matrix_complex * m, const double eps, const char *desc)
{
  int s = 0;
  size_t i, j, N = m->size1;

  gsl_matrix_complex * lu  = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * c  = gsl_matrix_complex_alloc(N, N);
  gsl_permutation * perm = gsl_permutation_alloc(N);
  int signum;

  gsl_matrix_complex_memcpy(lu, m);

  s += gsl_linalg_complex_LU_decomp(lu, perm, &signum);
  s += gsl_linalg_complex_LU_invx(lu, perm);

  /* c = m m^{-1} */
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, m, lu, GSL_COMPLEX_ZERO, c);

  /* c should be the identity matrix */

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex cij = gsl_matrix_complex_get(c, i, j);
          double expected = (i == j) ? 1.0 : 0.0;

          gsl_test_rel(GSL_REAL(cij), expected, eps, "%s real (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, GSL_REAL(cij), expected);

          gsl_test_rel(GSL_IMAG(cij), 0.0, eps, "%s imag (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, GSL_IMAG(cij), 0.0);
        }
    }

  gsl_matrix_complex_free(lu);
  gsl_matrix_complex_free(c);
  gsl_permutation_free(perm);

  return s;
}

static int
test_LUc_invert(gsl_rng * r)
{
  int s = 0;
  size_t n;

  for (n = 1; n <= 50; ++n)
    {
      gsl_matrix_complex * m = gsl_matrix_complex_alloc(n, n);

      create_random_complex_matrix(m, r);
      test_LUc_invert_eps(m, 1.0e5 * n * GSL_DBL_EPSILON, "complex_LU_invert random");

      gsl_matrix_complex_free(m);
    }

  return s;
}
