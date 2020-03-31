/* linalg/test_tri.c
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

static int
test_symmtd_decomp_eps(const gsl_matrix * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t N = m->size1;
  size_t i, j;

  gsl_matrix * Q = gsl_matrix_alloc(N, N);
  gsl_matrix * T = gsl_matrix_calloc(N, N);
  gsl_matrix * A  = gsl_matrix_alloc(N, N);
  gsl_matrix * B  = gsl_matrix_alloc(N, N);
  gsl_vector * tau = gsl_vector_alloc(N - 1);
  gsl_vector_view diag = gsl_matrix_diagonal(T);
  gsl_vector_view subdiag = gsl_matrix_subdiagonal(T, 1);
  gsl_vector_view superdiag = gsl_matrix_superdiagonal(T, 1);

  gsl_matrix_memcpy(A, m);
  s += gsl_linalg_symmtd_decomp(A, tau);

  s += gsl_linalg_symmtd_unpack(A, tau, Q, &diag.vector, &subdiag.vector);
  gsl_vector_memcpy(&superdiag.vector, &subdiag.vector);
  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Q, T, 0.0, A); /* A := Q T */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, Q, 0.0, B);   /* B := Q T Q^T */

  for (i = 0; i < N; i++)
    {
      for (j = 0; j <= i; j++)
        {
          double bij = gsl_matrix_get(B, i, j);
          double mij = gsl_matrix_get(m, i, j);

          gsl_test_rel(bij, mij, eps, "%s (%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, i, j, bij, mij);
        }
    }

  gsl_matrix_free(T);
  gsl_matrix_free(A);
  gsl_matrix_free(Q);
  gsl_matrix_free(B);
  gsl_vector_free(tau);

  return s;
}

static int
test_symmtd_decomp(gsl_rng * r)
{
  int s = 0;
  size_t N;

  for (N = 2; N <= 50; ++N)
    {
      gsl_matrix * A = gsl_matrix_alloc(N, N);

      create_symm_matrix(A, r);
      s += test_symmtd_decomp_eps(A, 1.0e5 * N * GSL_DBL_EPSILON, "symmtd_decomp random");

      gsl_matrix_free(A);
    }

  return s;
}

static int
test_hermtd_decomp_eps(const gsl_matrix_complex * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t N = m->size1;
  size_t i, j;

  gsl_matrix_complex * Q = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * T = gsl_matrix_complex_calloc(N, N);
  gsl_matrix_complex * A  = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * B  = gsl_matrix_complex_alloc(N, N);
  gsl_vector_complex * tau = gsl_vector_complex_alloc(N - 1);
  gsl_vector_view diag, subdiag, superdiag;
  gsl_vector_complex_view v;
  
  v = gsl_matrix_complex_diagonal(T);
  diag = gsl_vector_complex_real(&v.vector);

  v = gsl_matrix_complex_subdiagonal(T, 1);
  subdiag = gsl_vector_complex_real(&v.vector);

  v = gsl_matrix_complex_superdiagonal(T, 1);
  superdiag = gsl_vector_complex_real(&v.vector);

  gsl_matrix_complex_memcpy(A, m);
  s += gsl_linalg_hermtd_decomp(A, tau);

  s += gsl_linalg_hermtd_unpack(A, tau, Q, &diag.vector, &subdiag.vector);
  gsl_vector_memcpy(&superdiag.vector, &subdiag.vector);
  
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, Q, T, GSL_COMPLEX_ZERO, A);    /* A := Q T */
  gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, A, Q, GSL_COMPLEX_ZERO, B);  /* B := Q T Q^T */

  for (i = 0; i < N; i++)
    {
      for (j = 0; j <= i; j++)
        {
          gsl_complex bij = gsl_matrix_complex_get(B, i, j);
          gsl_complex mij = gsl_matrix_complex_get(m, i, j);

          gsl_test_rel(GSL_REAL(bij), GSL_REAL(mij), eps, "%s real (%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, i, j, GSL_REAL(bij), GSL_REAL(mij));
          gsl_test_rel(GSL_IMAG(bij), GSL_IMAG(mij), eps, "%s imag (%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, i, j, GSL_IMAG(bij), GSL_IMAG(mij));
        }
    }

  gsl_matrix_complex_free(T);
  gsl_matrix_complex_free(A);
  gsl_matrix_complex_free(Q);
  gsl_matrix_complex_free(B);
  gsl_vector_complex_free(tau);

  return s;
}

static int
test_hermtd_decomp(gsl_rng * r)
{
  int s = 0;
  size_t N;

  for (N = 2; N <= 50; ++N)
    {
      gsl_matrix_complex * A = gsl_matrix_complex_alloc(N, N);

      create_herm_matrix(A, r);
      s += test_hermtd_decomp_eps(A, 1.0e5 * N * GSL_DBL_EPSILON, "hermtd_decomp random");

      gsl_matrix_complex_free(A);
    }

  return s;
}
