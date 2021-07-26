/* linalg/test_ql.c
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
test_QL_decomp_eps(const gsl_matrix * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  gsl_matrix * QL = gsl_matrix_alloc(M, N);
  gsl_vector * tau = gsl_vector_alloc(N);
  gsl_matrix * A  = gsl_matrix_alloc(M, N);
  gsl_matrix * L  = gsl_matrix_alloc(M, N);
  gsl_matrix * Q  = gsl_matrix_alloc(M, M);

  gsl_matrix_memcpy(QL, m);

  s += gsl_linalg_QL_decomp(QL, tau);
  s += gsl_linalg_QL_unpack(QL, tau, Q, L);

  /* compute A = Q L */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, L, 0.0, A);

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

  gsl_matrix_free(QL);
  gsl_vector_free(tau);
  gsl_matrix_free(A);
  gsl_matrix_free(Q);
  gsl_matrix_free(L);

  return s;
}

static int
test_QL_decomp(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 30; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix * A = gsl_matrix_alloc(M, N);

          create_random_matrix(A, r);
          s += test_QL_decomp_eps(A, 1.0e5 * GSL_MAX(M, N) * GSL_DBL_EPSILON, "QL_decomp random");

          gsl_matrix_free(A);
        }
    }

  s += test_QL_decomp_eps(m53,   1.0e2 * GSL_DBL_EPSILON, "QL_decomp m(5,3)");
  s += test_QL_decomp_eps(hilb2, 1.0e2 * GSL_DBL_EPSILON, "QL_decomp hilbert(2)");
  s += test_QL_decomp_eps(hilb3, 1.0e2 * GSL_DBL_EPSILON, "QL_decomp hilbert(3)");
  s += test_QL_decomp_eps(hilb4, 1.0e2 * GSL_DBL_EPSILON, "QL_decomp hilbert(4)");
  s += test_QL_decomp_eps(hilb12, 1.0e2 * GSL_DBL_EPSILON, "QL_decomp hilbert(12)");
  s += test_QL_decomp_eps(vander2, 1.0e1 * GSL_DBL_EPSILON, "QL_decomp vander(2)");
  s += test_QL_decomp_eps(vander3, 1.0e1 * GSL_DBL_EPSILON, "QL_decomp vander(3)");
  s += test_QL_decomp_eps(vander4, 1.0e2 * GSL_DBL_EPSILON, "QL_decomp vander(4)");

  return s;
}
