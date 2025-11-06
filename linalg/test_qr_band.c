/* linalg/test_qr_band.c
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
test_QR_band_decomp_eps(const size_t p, const size_t q, const gsl_matrix * A, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i, j;

  gsl_matrix * AB = gsl_matrix_alloc(N, 2*p + q + 1);
  gsl_vector * tau = gsl_vector_alloc(N);
  gsl_matrix * Q = gsl_matrix_alloc(M, M);
  gsl_matrix * R = gsl_matrix_alloc(M, N);
  gsl_matrix * B = gsl_matrix_alloc(M, N);

  /* convert A to packed banded format */
  gen2band_matrix(p, q, A, AB);

#if 0
  print_octave(A, "A");
  printqrb_octave(M, p, q, AB, "AB");
#endif

  s += gsl_linalg_QR_band_decomp_L2(M, p, q, AB, tau);
  s += gsl_linalg_QR_band_unpack_L2(p, q, AB, tau, Q, R);

#if 0
  printqrb_octave(M, p, q, AB, "QR");
  printv_octave(tau, "tau");
  print_octave(Q, "Q");
  print_octave(R, "R");
#endif

  /* compute B = Q R */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, R, 0.0, B);

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double aij = gsl_matrix_get(A, i, j);
          double bij = gsl_matrix_get(B, i, j);

          gsl_test_rel(bij, aij, eps, "%s (M=%3lu,N=%3lu)(p=%3lu,q=%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, p, q, i, j, aij, bij);
        }
    }

  gsl_matrix_free(AB);
  gsl_vector_free(tau);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);
  gsl_matrix_free(B);

  return s;
}

static int
test_QR_band_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 20;
  size_t M, N, p, q;

  /*M = 7; N = 6; p = 2; q = 1;*/
  for (M = 1; M <= N_max; ++M)
    {
      for (N = 1; N <= N_max; ++N)
        {
          gsl_matrix * A = gsl_matrix_alloc(M, N);

          for (p = 0; p < GSL_MIN(M, 10); ++p)
            {
              for (q = 0; q < GSL_MIN(N, 10); ++q)
                {
                  create_band_matrix(p, q, A, r);
                  s += test_QR_band_decomp_eps(p, q, A, 1.0e6 * GSL_MAX(M,N) * GSL_DBL_EPSILON, "QR_band_decomp random");
                }
            }

          gsl_matrix_free(A);
        }
    }

  return s;
}

#if 0
static int
test_QR_band_decomp(void)
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
#endif
