/* linalg/test_cholesky.c
 *
 * Copyright (C) 2016, 2022 Patrick Alken
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

static int test_choleskyc_decomp_eps(const int scale, const gsl_matrix_complex * m, const double eps, const char * desc);
static int test_choleskyc_decomp(gsl_rng * r);
static int test_choleskyc_invert(gsl_rng * r);

static int
test_choleskyc_decomp_eps(const int scale, const gsl_matrix_complex * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  gsl_matrix_complex * v  = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * a  = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * l  = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * lh  = gsl_matrix_complex_alloc(N, N);
  gsl_vector * S = gsl_vector_alloc(N);

  gsl_matrix_complex_memcpy(v, m);
  gsl_matrix_complex_set_zero(l);
  gsl_matrix_complex_set_zero(lh);

  if (scale)
    s += gsl_linalg_complex_cholesky_decomp2(v, S);
  else
    s += gsl_linalg_complex_cholesky_decomp(v);

  /* Compute L L^H */
  
  for (i = 0; i < N ; i++)
    {
      for (j = 0; j <= i; j++)
        {
          gsl_complex vij = gsl_matrix_complex_get(v, i, j);
          gsl_matrix_complex_set (l, i, j, vij);
          gsl_matrix_complex_set (lh, j, i, gsl_complex_conjugate(vij));
        }
    }

  if (scale)
    {
      /* L <- S^{-1} L, LH <- LH S^{-1} */
      for (i = 0; i < N; ++i)
        {
          double Si = gsl_vector_get(S, i);
          gsl_vector_complex_view v = gsl_matrix_complex_row(l, i);
          gsl_vector_complex_view w = gsl_matrix_complex_column(lh, i);
          gsl_complex z;

          GSL_SET_COMPLEX(&z, 1.0 / Si, 0.0);

          gsl_vector_complex_scale(&v.vector, z);
          gsl_vector_complex_scale(&w.vector, z);
        }
    }
            
  /* compute a = l lh */
  gsl_blas_zgemm (CblasNoTrans,
                  CblasNoTrans,
                  GSL_COMPLEX_ONE,
                  l,
                  lh,
                  GSL_COMPLEX_ZERO,
                  a);

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          gsl_complex aij = gsl_matrix_complex_get(a, i, j);
          gsl_complex mij = gsl_matrix_complex_get(m, i, j);

          gsl_test_rel(GSL_REAL(aij), GSL_REAL(mij), eps,
                       "%s: real (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, GSL_REAL(aij), GSL_REAL(mij));

          gsl_test_rel(GSL_IMAG(aij), GSL_IMAG(mij), eps,
                       "%s: imag (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, GSL_IMAG(aij), GSL_IMAG(mij));
        }
    }

  gsl_matrix_complex_free(v);
  gsl_matrix_complex_free(a);
  gsl_matrix_complex_free(l);
  gsl_matrix_complex_free(lh);
  gsl_vector_free(S);

  return s;
}

static int
test_choleskyc_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix_complex * m = gsl_matrix_complex_alloc(N, N);

      create_posdef_complex_matrix(m, r);
      test_choleskyc_decomp_eps(0, m, 1.0e6 * N * GSL_DBL_EPSILON, "cholesky_complex_decomp unscaled random");
      test_choleskyc_decomp_eps(1, m, 1.0e6 * N * GSL_DBL_EPSILON, "cholesky_complex_decomp scaled random");

      gsl_matrix_complex_free(m);
    }

  return s;
}

static int
test_choleskyc_solve_eps(const int scale, const gsl_matrix_complex * m, const gsl_vector_complex * rhs,
                         const gsl_vector_complex * sol, const double eps, const char * desc)
{
  int s = 0;
  size_t i, N = m->size1;
  gsl_matrix_complex * u  = gsl_matrix_complex_alloc(N, N);
  gsl_vector_complex * x = gsl_vector_complex_calloc(N);
  gsl_vector * S = gsl_vector_alloc(N);

  gsl_matrix_complex_memcpy(u, m);

  if (scale)
    {
      s += gsl_linalg_complex_cholesky_decomp2(u, S);
      s += gsl_linalg_complex_cholesky_solve2(u, S, rhs, x);
    }
  else
    {
      s += gsl_linalg_complex_cholesky_decomp(u);
      s += gsl_linalg_complex_cholesky_solve(u, rhs, x);
    }

  for (i = 0; i < N; i++)
    {
      gsl_complex xi = gsl_vector_complex_get(x, i);
      gsl_complex yi = gsl_vector_complex_get(sol, i);

      gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps,
                   "%s: real %3lu[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
      gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps,
                   "%s: imag %3lu[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_vector_complex_free(x);
  gsl_vector_free(S);
  gsl_matrix_complex_free(u);

  return s;
}

static int
test_choleskyc_solve(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix_complex * m = gsl_matrix_complex_alloc(N, N);
      gsl_vector_complex * rhs = gsl_vector_complex_alloc(N);
      gsl_vector_complex * sol = gsl_vector_complex_alloc(N);

      create_posdef_complex_matrix(m, r);
      create_random_complex_vector(sol, r);
      gsl_blas_zhemv(CblasLower, GSL_COMPLEX_ONE, m, sol, GSL_COMPLEX_ZERO, rhs);

      test_choleskyc_solve_eps(0, m, rhs, sol, 64.0 * N * GSL_DBL_EPSILON, "cholesky_complex_solve unscaled random");
      test_choleskyc_solve_eps(1, m, rhs, sol, 64.0 * N * GSL_DBL_EPSILON, "cholesky_complex_solve scaled random");

      gsl_matrix_complex_free(m);
      gsl_vector_complex_free(rhs);
      gsl_vector_complex_free(sol);
    }

  return s;
}

static int
test_choleskyc_invert_eps(const gsl_matrix_complex * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t N = m->size1;
  size_t i, j;
  gsl_matrix_complex * v  = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * c  = gsl_matrix_complex_alloc(N, N);

  gsl_matrix_complex_memcpy(v, m);

  s += gsl_linalg_complex_cholesky_decomp(v);
  s += gsl_linalg_complex_cholesky_invert(v);

  gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, m, v, GSL_COMPLEX_ZERO, c);

  /* c should be the identity matrix */
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex cij = gsl_matrix_complex_get(c, i, j);
          double expected = (i == j) ? 1.0 : 0.0;

          /* check real part */
          gsl_test_rel(GSL_REAL(cij), expected, eps,
                       "%s: real (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, GSL_REAL(cij), expected);

          /* check imaginary part */
          gsl_test_rel(GSL_IMAG(cij), 0.0, eps,
                       "%s: imag (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, GSL_IMAG(cij), 0.0);
        }
    }

  gsl_matrix_complex_free(v);
  gsl_matrix_complex_free(c);

  return s;
}

static int
test_choleskyc_invert(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix_complex * m = gsl_matrix_complex_alloc(N, N);

      create_posdef_complex_matrix(m, r);
      test_choleskyc_invert_eps(m, 32.0 * N * GSL_DBL_EPSILON, "cholesky_complex_invert random");

      gsl_matrix_complex_free(m);
    }

  return s;
}
