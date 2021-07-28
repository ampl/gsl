/* linalg/test_cholesky.c
 *
 * Copyright (C) 2016 Patrick Alken
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

static int test_cholesky_decomp_eps(const int scale, const gsl_matrix * m,
                                    const double expected_rcond, const double eps,
                                    const char * desc);
static int test_cholesky_decomp(gsl_rng * r);
int test_cholesky_invert_eps(const gsl_matrix * m, const double eps, const char *desc);
int test_cholesky_invert(gsl_rng * r);
static int test_pcholesky_decomp_eps(const int scale, const gsl_matrix * m,
                                     const double expected_rcond, const double eps,
                                     const char * desc);
static int test_pcholesky_decomp(gsl_rng * r);
int test_pcholesky_solve_eps(const int scale, const gsl_matrix * m, const gsl_vector * rhs,
                             const gsl_vector * sol, const double eps,
                             const char * desc);
static int test_pcholesky_solve(gsl_rng * r);
int test_pcholesky_invert_eps(const gsl_matrix * m, const double eps, const char *desc);
int test_pcholesky_invert(gsl_rng * r);

static int test_mcholesky_decomp_eps(const int posdef, const int scale, const gsl_matrix * m,
                                     const double expected_rcond, const double eps, const char * desc);

/* Hilbert matrix condition numbers, as calculated by LAPACK DPOSVX */
double hilb_rcond[] = { 1.000000000000e+00, 3.703703703704e-02, 1.336898395722e-03,
                        3.524229074890e-05, 1.059708198754e-06, 3.439939465186e-08,
                        1.015027593823e-09, 2.952221630602e-11, 9.093751565191e-13,
                        2.828277420229e-14, 8.110242564869e-16, 2.409320075800e-17 };

static double
test_cholesky_norm1(const gsl_matrix * m)
{
  const size_t N = m->size2;
  double value = 0.0;
  size_t j;

  for (j = 0; j < N; ++j)
    {
      gsl_vector_const_view v = gsl_matrix_const_column(m, j);
      double sum = gsl_blas_dasum(&v.vector);
      value = GSL_MAX(value, sum);
    }

  return value;
}

static int
test_cholesky_decomp_eps(const int scale, const gsl_matrix * m,
                         const double expected_rcond, const double eps,
                         const char * desc)
{
  int s = 0;
  size_t i, j, N = m->size2;

  gsl_matrix * V  = gsl_matrix_alloc(N, N);
  gsl_matrix * A  = gsl_matrix_alloc(N, N);
  gsl_matrix * L  = gsl_matrix_calloc(N, N);
  gsl_matrix * LT = gsl_matrix_calloc(N, N);
  gsl_vector * S = gsl_vector_alloc(N);

  gsl_matrix_memcpy(V, m);

  if (scale)
    s += gsl_linalg_cholesky_decomp2(V, S);
  else
    s += gsl_linalg_cholesky_decomp1(V);

  /* compute L and LT */
  gsl_matrix_tricpy(CblasLower, CblasNonUnit, L, V);
  gsl_matrix_transpose_tricpy(CblasLower, CblasNonUnit, LT, L);
  
  if (scale)
    {
      /* L <- S^{-1} L, LT <- LT S^{-1} */
      for (i = 0; i < N; ++i)
        {
          double Si = gsl_vector_get(S, i);
          gsl_vector_view v = gsl_matrix_row(L, i);
          gsl_vector_view w = gsl_matrix_column(LT, i);

          gsl_vector_scale(&v.vector, 1.0 / Si);
          gsl_vector_scale(&w.vector, 1.0 / Si);
        }
    }
            
  /* compute A = L LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j);
          double mij = gsl_matrix_get(m, i, j);

          gsl_test_rel(Aij, mij, eps,
                       "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, Aij, mij);
        }
    }

  if (expected_rcond > 0 && !scale)
    {
      gsl_vector *work = gsl_vector_alloc(3 * N);
      double rcond;

      gsl_linalg_cholesky_rcond(V, &rcond, work);

      gsl_test_rel(rcond, expected_rcond, 1.0e-6,
                   "%s rcond: (%3lu,%3lu): %22.18g   %22.18g\n",
                   desc, N, N, rcond, expected_rcond);

      gsl_vector_free(work);
    }

  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);
  gsl_vector_free(S);

  return s;
}

static int
test_cholesky_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 200;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);
      test_cholesky_decomp_eps(0, m, -1.0, 1.0e3 * N * GSL_DBL_EPSILON, "cholesky_decomp unscaled random");
      test_cholesky_decomp_eps(1, m, -1.0, 1.0e3 * N * GSL_DBL_EPSILON, "cholesky_decomp scaled random");

      if (N <= 12)
        {
          double expected_rcond = -1.0;
          
          if (hilb_rcond[N - 1] > 1.0e-12)
            expected_rcond = hilb_rcond[N - 1];

          create_hilbert_matrix2(m);

          test_cholesky_decomp_eps(0, m, expected_rcond, N * GSL_DBL_EPSILON, "cholesky_decomp unscaled hilbert");
          test_cholesky_decomp_eps(1, m, expected_rcond, N * GSL_DBL_EPSILON, "cholesky_decomp scaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

int
test_cholesky_invert_eps(const gsl_matrix * m, const double eps, const char *desc)
{
  int s = 0;
  size_t i, j, N = m->size1;

  gsl_matrix * v  = gsl_matrix_alloc(N, N);
  gsl_matrix * c  = gsl_matrix_alloc(N, N);

  gsl_matrix_memcpy(v, m);

  s += gsl_linalg_cholesky_decomp1(v);
  s += gsl_linalg_cholesky_invert(v);

  /* c = m m^{-1} */
  gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, m, v, 0.0, c);

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

  gsl_matrix_free(v);
  gsl_matrix_free(c);

  return s;
}

int
test_cholesky_invert(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 200;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);

      test_cholesky_invert_eps(m, 1.0e2 * N * GSL_DBL_EPSILON, "cholesky_invert unscaled random");

      if (N <= 4)
        {
          create_hilbert_matrix2(m);
          test_cholesky_invert_eps(m, N * 256.0 * GSL_DBL_EPSILON, "cholesky_invert unscaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

static int
test_mcholesky_decomp_eps(const int posdef, const int scale, const gsl_matrix * m,
                          const double expected_rcond, const double eps, const char * desc)
{
  int s = 0;
  size_t i, j, N = m->size2;

  gsl_matrix * LDLT = gsl_matrix_alloc(N, N);
  gsl_matrix * V  = gsl_matrix_alloc(N, N);
  gsl_matrix * A  = gsl_matrix_alloc(N, N);
  gsl_matrix * L  = gsl_matrix_alloc(N, N);
  gsl_matrix * LT = gsl_matrix_alloc(N, N);
  gsl_vector * S = gsl_vector_alloc(N);
  gsl_vector * E = gsl_vector_alloc(N);
  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_vector_view D = gsl_matrix_diagonal(LDLT);

  gsl_matrix_memcpy(LDLT, m);
  s += gsl_linalg_mcholesky_decomp(LDLT, perm, E);

  /* check that the upper triangle of LDLT equals original matrix */
  for (i = 0; i < N; ++i)
    {
      for (j = i + 1; j < N; ++j)
        {
          double mij = gsl_matrix_get(m, i, j);
          double aij = gsl_matrix_get(LDLT, i, j);

          gsl_test_rel(aij, mij, 1.0e-12,
                       "%s upper triangle: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, aij, mij);
        }
    }

  if (posdef)
    {
      /* ||E|| should be 0 */
      double norm = gsl_blas_dnrm2(E);
      s = norm != 0.0;
      gsl_test(s, "%s: (%zu,%zu): ||E|| = %.12e",
               desc, N, N, norm);

      /* check that D is decreasing */
      s = 0;
      for (i = 1; i < N; ++i)
        {
          double dprev = gsl_vector_get(&D.vector, i - 1);
          double di = gsl_vector_get(&D.vector, i);

          if (di > dprev)
            s = 1;
        }

      gsl_test(s, "%s: (%zu,%zu): D is not decreasing",
               desc, N, N);
    }
  
  /* compute L and LT */
  gsl_matrix_set_identity(L);
  gsl_matrix_set_identity(LT);

  gsl_matrix_tricpy(CblasLower, CblasUnit, L, LDLT);
  gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, LT, L);

  /* compute (L sqrt(D)) and (sqrt(D) LT) */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_column(L, i);
      gsl_vector_view w = gsl_matrix_row(LT, i);
      double di = gsl_vector_get(&D.vector, i);

      gsl_vector_scale(&v.vector, sqrt(di));
      gsl_vector_scale(&w.vector, sqrt(di));
    }

  /* compute A = L D LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  /* compute V = P (S M S + E) P^T */

  gsl_matrix_memcpy(V, m);
  D = gsl_matrix_diagonal(V);

  /* compute S M S */
  if (scale)
    {
      gsl_linalg_cholesky_scale_apply(V, S);
      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, V, V);
    }

  /* compute S M S + E */
  gsl_vector_add(&D.vector, E);

  /* compute M P^T */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(V, i);
      gsl_permute_vector(perm, &v.vector);
    }

  /* compute P M P^T */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_column(V, i);
      gsl_permute_vector(perm, &v.vector);
    }

  for (i = 0; i < N; i++)
    {
      double Ei = gsl_vector_get(E, i);

      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j); /* L D L^T */
          double Bij = gsl_matrix_get(V, i, j); /* P M P^T */
          double Cij;                           /* P M P^T + E */

          if (i == j)
            Cij = Bij + Ei*0;
          else
            Cij = Bij;

          gsl_test_rel(Aij, Cij, eps,
                       "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, Aij, Cij);
        }
    }

  if (expected_rcond > 0 && !scale)
    {
      gsl_vector *work = gsl_vector_alloc(3 * N);
      double rcond;

      gsl_linalg_mcholesky_rcond(LDLT, perm, &rcond, work);

      gsl_test_rel(rcond, expected_rcond, 1.0e-6,
                   "%s rcond: (%3lu,%3lu): %22.18g   %22.18g\n",
                   desc, N, N, rcond, expected_rcond);

      gsl_vector_free(work);
    }

  gsl_matrix_free(LDLT);
  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);
  gsl_vector_free(S);
  gsl_vector_free(E);
  gsl_permutation_free(perm);

  return s;
}

static int
test_mcholesky_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);
      test_mcholesky_decomp_eps(1, 0, m, -1.0, 1.0e3 * N * GSL_DBL_EPSILON, "mcholesky_decomp unscaled random posdef");

      create_symm_matrix(m, r);
      test_mcholesky_decomp_eps(0, 0, m, -1.0, 1.0e5 * N * GSL_DBL_EPSILON, "mcholesky_decomp unscaled random symm");

      if (N <= 8)
        {
          double expected_rcond = -1.0;
          
          if (hilb_rcond[N - 1] > 1.0e-12)
            expected_rcond = hilb_rcond[N - 1];

          create_hilbert_matrix2(m);

          test_mcholesky_decomp_eps(1, 0, m, expected_rcond, 128.0 * N * GSL_DBL_EPSILON, "mcholesky_decomp unscaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

int
test_mcholesky_solve_eps(const gsl_matrix * m, const gsl_vector * rhs,
                         const gsl_vector * sol, const double eps,
                         const char * desc)
{
  int s = 0;
  size_t i, N = m->size1;
  gsl_matrix * u  = gsl_matrix_alloc(N, N);
  gsl_vector * x = gsl_vector_calloc(N);
  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_matrix_memcpy(u, m);

  s += gsl_linalg_mcholesky_decomp(u, perm, NULL);
  s += gsl_linalg_mcholesky_solve(u, perm, rhs, x);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps,
                   "%s: %3lu[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_vector_free(x);
  gsl_matrix_free(u);
  gsl_permutation_free(perm);

  return s;
}

static int
test_mcholesky_solve(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);
      gsl_vector * rhs = gsl_vector_alloc(N);
      gsl_vector * sol = gsl_vector_alloc(N);

      create_posdef_matrix(m, r);
      create_random_vector(sol, r);
      gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);

      test_mcholesky_solve_eps(m, rhs, sol, 64.0 * N * GSL_DBL_EPSILON, "mcholesky_solve random");

      if (N <= 3)
        {
          create_hilbert_matrix2(m);
          gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);
          test_mcholesky_solve_eps(m, rhs, sol, 1.0e3 * N * GSL_DBL_EPSILON, "mcholesky_solve hilbert");
        }

      gsl_matrix_free(m);
      gsl_vector_free(rhs);
      gsl_vector_free(sol);
    }

  return s;
}

int
test_mcholesky_invert_eps(const gsl_matrix * m, const double eps, const char *desc)
{
  int s = 0;
  size_t i, j, N = m->size1;

  gsl_matrix * v  = gsl_matrix_alloc(N, N);
  gsl_matrix * c  = gsl_matrix_alloc(N, N);
  gsl_matrix * minv = gsl_matrix_alloc(N, N);
  gsl_vector * E = gsl_vector_alloc(N);
  gsl_permutation * p = gsl_permutation_alloc(N);

  gsl_matrix_memcpy(v, m);

  s += gsl_linalg_mcholesky_decomp(v, p, E);
  s += gsl_linalg_mcholesky_invert(v, p, minv);

  /* c = m m^{-1} */
  gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, m, minv, 0.0, c);

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

  gsl_matrix_free(v);
  gsl_matrix_free(c);
  gsl_matrix_free(minv);
  gsl_vector_free(E);
  gsl_permutation_free(p);

  return s;
}

int
test_mcholesky_invert(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 30;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);
      test_mcholesky_invert_eps(m, N * GSL_DBL_EPSILON, "mcholesky_invert unscaled random");

      if (N <= 4)
        {
          create_hilbert_matrix2(m);
          test_mcholesky_invert_eps(m, 512.0 * N * GSL_DBL_EPSILON, "mcholesky_invert unscaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

static int
test_pcholesky_decomp_eps(const int scale, const gsl_matrix * m,
                          const double expected_rcond, const double eps,
                          const char * desc)
{
  int s = 0;
  size_t i, j, N = m->size2;

  gsl_matrix * LDLT = gsl_matrix_alloc(N, N);
  gsl_matrix * V  = gsl_matrix_alloc(N, N);
  gsl_matrix * A  = gsl_matrix_alloc(N, N);
  gsl_matrix * L  = gsl_matrix_alloc(N, N);
  gsl_matrix * LT = gsl_matrix_alloc(N, N);
  gsl_vector * S = gsl_vector_alloc(N);
  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_vector_view D = gsl_matrix_diagonal(LDLT);

  gsl_matrix_memcpy(LDLT, m);

  if (scale)
    s += gsl_linalg_pcholesky_decomp2(LDLT, perm, S);
  else
    s += gsl_linalg_pcholesky_decomp(LDLT, perm);

  /* check that the upper triangle of LDLT equals original matrix */
  for (i = 0; i < N; ++i)
    {
      for (j = i + 1; j < N; ++j)
        {
          double mij = gsl_matrix_get(m, i, j);
          double aij = gsl_matrix_get(LDLT, i, j);

          gsl_test_rel(aij, mij, 1.0e-12,
                       "%s upper triangle: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, aij, mij);
        }
    }

  /* check that D is decreasing */
  s = 0;
  for (i = 1; i < N; ++i)
    {
      double dprev = gsl_vector_get(&D.vector, i - 1);
      double di = gsl_vector_get(&D.vector, i);

      if (di > dprev)
        s = 1;
    }

  gsl_test(s, "%s: (%zu,%zu): D is not decreasing",
           desc, N, N);
  
  /* compute L and LT */
  gsl_matrix_set_identity(L);
  gsl_matrix_set_identity(LT);

  gsl_matrix_tricpy(CblasLower, CblasUnit, L, LDLT);
  gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, LT, L);

  /* compute (L sqrt(D)) and (sqrt(D) LT) */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_column(L, i);
      gsl_vector_view w = gsl_matrix_row(LT, i);
      double di = gsl_vector_get(&D.vector, i);

      gsl_vector_scale(&v.vector, sqrt(di));
      gsl_vector_scale(&w.vector, sqrt(di));
    }

  /* compute A = L D LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  /* compute V = P S M S P^T */

  gsl_matrix_memcpy(V, m);

  /* compute S M S */
  if (scale)
    {
      gsl_linalg_cholesky_scale_apply(V, S);
      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, V, V);
    }

  /* compute M P^T */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(V, i);
      gsl_permute_vector(perm, &v.vector);
    }

  /* compute P M P^T */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_column(V, i);
      gsl_permute_vector(perm, &v.vector);
    }

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j); /* L D L^T */
          double Bij = gsl_matrix_get(V, i, j); /* P M P^T */

          gsl_test_rel(Aij, Bij, eps,
                       "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, Aij, Bij);
        }
    }

  if (expected_rcond > 0 && !scale)
    {
      gsl_vector *work = gsl_vector_alloc(3 * N);
      double rcond;

      gsl_linalg_pcholesky_rcond(LDLT, perm, &rcond, work);

      gsl_test_rel(rcond, expected_rcond, 1.0e-6,
                   "%s rcond: (%3lu,%3lu): %22.18g   %22.18g\n",
                   desc, N, N, rcond, expected_rcond);

      gsl_vector_free(work);
    }

  gsl_matrix_free(LDLT);
  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);
  gsl_vector_free(S);
  gsl_permutation_free(perm);

  return s;
}

static int
test_pcholesky_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);
      test_pcholesky_decomp_eps(0, m, -1.0, 1024.0 * N * GSL_DBL_EPSILON, "pcholesky_decomp unscaled random");
      test_pcholesky_decomp_eps(1, m, -1.0, 1024.0 * N * GSL_DBL_EPSILON, "pcholesky_decomp scaled random");

      if (N <= 12)
        {
          double expected_rcond = -1.0;
          
          if (hilb_rcond[N - 1] > 1.0e-12)
            expected_rcond = hilb_rcond[N - 1];

          create_hilbert_matrix2(m);

          test_pcholesky_decomp_eps(0, m, expected_rcond, 1024.0 * N * GSL_DBL_EPSILON, "pcholesky_decomp unscaled hilbert");
          test_pcholesky_decomp_eps(1, m, expected_rcond, 1024.0 * N * GSL_DBL_EPSILON, "pcholesky_decomp scaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

int
test_pcholesky_solve_eps(const int scale, const gsl_matrix * m, const gsl_vector * rhs,
                         const gsl_vector * sol, const double eps,
                         const char * desc)
{
  int s = 0;
  size_t i, N = m->size1;
  gsl_matrix * u  = gsl_matrix_alloc(N, N);
  gsl_vector * x = gsl_vector_calloc(N);
  gsl_vector * S = gsl_vector_alloc(N);
  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_matrix_memcpy(u, m);

  if (scale)
    {
      s += gsl_linalg_pcholesky_decomp2(u, perm, S);
      s += gsl_linalg_pcholesky_solve2(u, perm, S, rhs, x);
    }
  else
    {
      s += gsl_linalg_pcholesky_decomp(u, perm);
      s += gsl_linalg_pcholesky_solve(u, perm, rhs, x);
    }

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps,
                   "%s: %3lu[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_vector_free(x);
  gsl_vector_free(S);
  gsl_matrix_free(u);
  gsl_permutation_free(perm);

  return s;
}

static int
test_pcholesky_solve(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);
      gsl_vector * rhs = gsl_vector_alloc(N);
      gsl_vector * sol = gsl_vector_alloc(N);

      create_posdef_matrix(m, r);
      create_random_vector(sol, r);
      gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);

      test_pcholesky_solve_eps(0, m, rhs, sol, 64.0 * N * GSL_DBL_EPSILON, "pcholesky_solve unscaled random");
      test_pcholesky_solve_eps(1, m, rhs, sol, 64.0 * N * GSL_DBL_EPSILON, "pcholesky_solve scaled random");

      if (N <= 3)
        {
          create_hilbert_matrix2(m);
          gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);
          test_pcholesky_solve_eps(0, m, rhs, sol, 1024.0 * N * GSL_DBL_EPSILON, "pcholesky_solve unscaled hilbert");
          test_pcholesky_solve_eps(1, m, rhs, sol, 2048.0 * N * GSL_DBL_EPSILON, "pcholesky_solve scaled hilbert");
        }

      gsl_matrix_free(m);
      gsl_vector_free(rhs);
      gsl_vector_free(sol);
    }

  return s;
}

int
test_pcholesky_invert_eps(const gsl_matrix * m, const double eps, const char *desc)
{
  int s = 0;
  size_t i, j, N = m->size1;

  gsl_matrix * v  = gsl_matrix_alloc(N, N);
  gsl_matrix * c  = gsl_matrix_alloc(N, N);
  gsl_matrix * minv = gsl_matrix_alloc(N, N);
  gsl_permutation * p = gsl_permutation_alloc(N);

  gsl_matrix_memcpy(v, m);

  s += gsl_linalg_pcholesky_decomp(v, p);
  s += gsl_linalg_pcholesky_invert(v, p, minv);

  /* c = m m^{-1} */
  gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, m, minv, 0.0, c);

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

  gsl_matrix_free(v);
  gsl_matrix_free(c);
  gsl_matrix_free(minv);
  gsl_permutation_free(p);

  return s;
}

int
test_pcholesky_invert(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 200;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);
      test_pcholesky_invert_eps(m, 1.0e2 * N * GSL_DBL_EPSILON, "pcholesky_invert unscaled random");

      if (N <= 4)
        {
          create_hilbert_matrix2(m);
          test_pcholesky_invert_eps(m, 1024.0 * N * GSL_DBL_EPSILON, "pcholesky_invert unscaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

static int
test_cholesky_band_decomp_eps(const size_t p, const gsl_matrix * m, const double eps, const char * desc)
{
  int s = 0;
  size_t i, j, N = m->size2;
  double rcond_expected, rcond;

  gsl_matrix * V  = gsl_matrix_alloc(N, p + 1);
  gsl_matrix * A  = gsl_matrix_alloc(N, N);
  gsl_matrix * L  = gsl_matrix_calloc(N, N);
  gsl_matrix * LT = gsl_matrix_calloc(N, N);
  gsl_vector * work = gsl_vector_alloc(3 * N);

  /* convert m to packed banded format */
  symm2band_matrix(p, m, V);

  s += gsl_linalg_cholesky_band_decomp(V);

  /* compute L and LT */
  gsl_linalg_cholesky_band_unpack(V, L);
  gsl_matrix_transpose_tricpy(CblasLower, CblasNonUnit, LT, L);
  
  /* compute A = L LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j);
          double mij = gsl_matrix_get(m, i, j);

          gsl_test_rel(Aij, mij, eps,
                       "%s: (p=%zu,N=%zu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, p, N, i, j, Aij, mij);
        }
    }

  /* test 1-norm calculation */
  if (p > 0)
    {
      double norm1_expected = test_cholesky_norm1(m);
      double norm1 = gsl_matrix_get(V, N - 1, p);

      gsl_test_rel(norm1, norm1_expected, eps,
                   "%s: (p=%zu,N=%zu) 1-norm: %22.18g   %22.18g\n",
                   desc, p, N, norm1, norm1_expected);
    }

  /* test rcond */
  gsl_matrix_memcpy(A, m);
  s += gsl_linalg_cholesky_decomp1(A);
  s += gsl_linalg_cholesky_rcond(A, &rcond_expected, work);
  s += gsl_linalg_cholesky_band_rcond(V, &rcond, work);
  gsl_test_rel(rcond, rcond_expected, eps,
               "%s: (p=%zu,N=%zu) rcond: %22.18g   %22.18g\n",
               desc, p, N, rcond, rcond_expected);

  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);
  gsl_vector_free(work);

  return s;
}

static int
test_cholesky_band_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N, p;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      for (p = 0; p < GSL_MIN(N, 10); ++p)
        {
          create_posdef_band_matrix(p, m, r);
          s += test_cholesky_band_decomp_eps(p, m, 1.0e5 * N * GSL_DBL_EPSILON, "cholesky_band_decomp random");
        }

      gsl_matrix_free(m);
    }

  return s;
}

int
test_cholesky_band_solve_eps(const int scale, const size_t p, const gsl_matrix * m,
                             const gsl_vector * rhs, const gsl_vector * sol,
                             const double eps, const char * desc)
{
  int status = 0;
  size_t i, N = m->size1;
  gsl_matrix * u  = gsl_matrix_alloc(N, p + 1);
  gsl_vector * x = gsl_vector_calloc(N);
  gsl_vector * s = gsl_vector_alloc(N);

  /* convert m to packed banded format */
  symm2band_matrix(p, m, u);

  if (scale)
    {
      gsl_linalg_cholesky_band_scale(u, s);
      gsl_linalg_cholesky_band_scale_apply(u, s);
    }

  status += gsl_linalg_cholesky_band_decomp(u);

  if (scale)
    {
      gsl_vector_memcpy(x, rhs);
      gsl_vector_mul(x, s);
      status += gsl_linalg_cholesky_band_svx(u, x);
      gsl_vector_mul(x, s);
    }
  else
    {
      status += gsl_linalg_cholesky_band_solve(u, rhs, x);
    }

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps,
                   "%s: p=%zu N=%zu [%lu]: %22.18g   %22.18g\n",
                   desc, p, N, i, xi, yi);
    }

  gsl_vector_free(x);
  gsl_vector_free(s);
  gsl_matrix_free(u);

  return status;
}

static int
test_cholesky_band_solve(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N, p;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);
      gsl_vector * rhs = gsl_vector_alloc(N);
      gsl_vector * sol = gsl_vector_alloc(N);

      for (p = 0; p < GSL_MIN(N, 10); ++p)
        {
          create_posdef_band_matrix(p, m, r);
          create_random_vector(sol, r);
          gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);

          test_cholesky_band_solve_eps(0, p, m, rhs, sol, 1.0e3 * N * GSL_DBL_EPSILON,
                                       "cholesky_band_solve random unscaled");

          test_cholesky_band_solve_eps(1, p, m, rhs, sol, 1.0e3 * N * GSL_DBL_EPSILON,
                                       "cholesky_band_solve random scaled");
        }

      gsl_matrix_free(m);
      gsl_vector_free(rhs);
      gsl_vector_free(sol);
    }

  return s;
}

int
test_cholesky_band_invert_eps(const size_t p, const gsl_matrix * m, const double eps, const char *desc)
{
  int s = 0;
  size_t i, j, N = m->size1;

  gsl_matrix * v  = gsl_matrix_alloc(N, p + 1);
  gsl_matrix * minv  = gsl_matrix_alloc(N, N);
  gsl_matrix * c  = gsl_matrix_alloc(N, N);

  /* convert m to packed banded format */
  symm2band_matrix(p, m, v);

  s += gsl_linalg_cholesky_band_decomp(v);
  s += gsl_linalg_cholesky_band_invert(v, minv);

  /* c = m m^{-1} */
  gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, m, minv, 0.0, c);

  /* c should be the identity matrix */

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double cij = gsl_matrix_get(c, i, j);
          double expected = (i == j) ? 1.0 : 0.0;

          gsl_test_rel(cij, expected, eps, "%s (p=%zu,N=%zu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, p, N, i, j, cij, expected);
        }
    }

  gsl_matrix_free(v);
  gsl_matrix_free(minv);
  gsl_matrix_free(c);

  return s;
}

int
test_cholesky_band_invert(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N, p;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      for (p = 0; p < GSL_MIN(N, 10); ++p)
        {
          create_posdef_band_matrix(p, m, r);
          test_cholesky_band_invert_eps(p, m, N * 1.0e1 * GSL_DBL_EPSILON, "cholesky_band_invert random");
        }

      gsl_matrix_free(m);
    }

  return s;
}
