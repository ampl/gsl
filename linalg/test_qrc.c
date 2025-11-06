/* linalg/test_qrc.c
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
test_QRc_decomp_eps(const gsl_matrix_complex * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  gsl_matrix_complex * qr = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * a  = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * q  = gsl_matrix_complex_alloc(M, M);
  gsl_matrix_complex * r  = gsl_matrix_complex_alloc(M, N);
  gsl_vector_complex * d = gsl_vector_complex_alloc(N);

  gsl_matrix_complex_memcpy(qr, m);

  s += gsl_linalg_complex_QR_decomp(qr, d);
  s += gsl_linalg_complex_QR_unpack(qr, d, q, r);
  
  /* compute a = q r */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, q, r, GSL_COMPLEX_ZERO, a);

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      gsl_complex aij = gsl_matrix_complex_get(a, i, j);
      gsl_complex mij = gsl_matrix_complex_get(m, i, j);
      int foo_r = check(GSL_REAL(aij), GSL_REAL(mij), eps);
      int foo_i = check(GSL_IMAG(aij), GSL_IMAG(mij), eps);
      if (foo_r || foo_i)
        {
          printf("%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n", desc, M, N, i, j, GSL_REAL(aij), GSL_REAL(mij));
          printf("%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n", desc, M, N, i, j, GSL_IMAG(aij), GSL_IMAG(mij));
        }
      s += foo_r + foo_i;
    }
  }

  gsl_vector_complex_free(d);
  gsl_matrix_complex_free(qr);
  gsl_matrix_complex_free(a);
  gsl_matrix_complex_free(q);
  gsl_matrix_complex_free(r);

  return s;
}

static int
test_QRc_decomp(gsl_rng * r)
{
  int s = 0;
  size_t N, M;

  /* test random matrices */
  for (N = 1; N <= 20; ++N)
    {
      for (M = 1; M <= 20; ++M)
        {
          gsl_matrix_complex * A = gsl_matrix_complex_alloc(N, M);

          create_random_complex_matrix(A, r);
          test_QRc_decomp_eps(A, 1.0e6 * GSL_MAX(N, M) * GSL_DBL_EPSILON, "complex_QR_decomp random");

          gsl_matrix_complex_free(A);
        }
    }

  return s;
}

static int
test_QRc_solve_eps(const gsl_matrix_complex * A, const gsl_vector_complex * rhs, const gsl_vector_complex * sol,
                   const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i;

  gsl_matrix_complex * QR = gsl_matrix_complex_alloc(M, N);
  gsl_vector_complex * tau = gsl_vector_complex_alloc(N);
  gsl_vector_complex * x  = gsl_vector_complex_alloc(N);

  gsl_matrix_complex_memcpy(QR, A);

  s += gsl_linalg_complex_QR_decomp(QR, tau);
  s += gsl_linalg_complex_QR_solve(QR, tau, rhs, x);

  for (i = 0; i < M; i++)
    {
      gsl_complex xi = gsl_vector_complex_get(x, i);
      gsl_complex yi = gsl_vector_complex_get(sol, i);

      gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps, "%s real (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
      gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps, "%s imag (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_matrix_complex_free(QR);
  gsl_vector_complex_free(tau);
  gsl_vector_complex_free(x);

  return s;
}

static int
test_QRc_solve(gsl_rng * r)
{
  int s = 0;
  size_t N;

  for (N = 1; N <= 50; ++N)
    {
      gsl_matrix_complex * A = gsl_matrix_complex_alloc(N, N);
      gsl_vector_complex * sol = gsl_vector_complex_alloc(N);
      gsl_vector_complex * rhs = gsl_vector_complex_alloc(N);

      create_random_complex_matrix(A, r);
      create_random_complex_vector(sol, r);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, A, sol, GSL_COMPLEX_ZERO, rhs);

      s += test_QRc_solve_eps(A, rhs, sol, 1.0e5 * N * GSL_DBL_EPSILON, "complex QR_solve_r random");

      gsl_matrix_complex_free(A);
      gsl_vector_complex_free(sol);
      gsl_vector_complex_free(rhs);
    }

  return s;
}

static int
test_QRc_lssolve_eps(const gsl_matrix_complex * A, const gsl_vector_complex * rhs, const gsl_vector_complex * sol,
                     const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i;

  gsl_matrix_complex * QR = gsl_matrix_complex_alloc(M, N);
  gsl_vector_complex * tau = gsl_vector_complex_alloc(N);
  gsl_vector_complex * x  = gsl_vector_complex_alloc(N);
  gsl_vector_complex * r  = gsl_vector_complex_alloc(M);
  gsl_vector_complex * r2  = gsl_vector_complex_alloc(M);

  gsl_matrix_complex_memcpy(QR, A);

  s += gsl_linalg_complex_QR_decomp(QR, tau);
  s += gsl_linalg_complex_QR_lssolve(QR, tau, rhs, x, r);

  for (i = 0; i < N; i++)
    {
      gsl_complex xi = gsl_vector_complex_get(x, i);
      gsl_complex yi = gsl_vector_complex_get(sol, i);

      gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps, "%s sol real (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, GSL_REAL(xi), GSL_REAL(yi));
      gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps, "%s sol imag (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, GSL_IMAG(xi), GSL_IMAG(yi));
    }

  /* compute residual and check */
  gsl_vector_complex_memcpy(r2, rhs);
  gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_NEGONE, A, sol, GSL_COMPLEX_ONE, r2);

  for (i = 0; i < M; i++)
    {
      gsl_complex xi = gsl_vector_complex_get(r, i);
      gsl_complex yi = gsl_vector_complex_get(r2, i);

      gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps, "%s res real (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, GSL_REAL(xi), GSL_REAL(yi));
      gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps, "%s res imag (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, GSL_IMAG(xi), GSL_IMAG(yi));
    }

  gsl_matrix_complex_free(QR);
  gsl_vector_complex_free(tau);
  gsl_vector_complex_free(x);
  gsl_vector_complex_free(r);
  gsl_vector_complex_free(r2);

  return s;
}

static int
test_QRc_lssolve()
{
  int s = 0;
  const double tol = 1.0e-10;

  {
    const double A_data[] = { 6.70178096995592e-01, 2.51765675689489e-01, 7.18387884271362e-01, 3.61106008613644e-01,
                              3.36342428526987e-02, 5.25755806084206e-01, 5.87596968507183e-01, 5.27700895952418e-01,
                              6.85871226549483e-01, 3.71216051930735e-01, 4.40045953339814e-01, 2.08875308141557e-01,
                              8.62676230425306e-01, 5.79995188556082e-01, 4.86983443637474e-02, 2.82915411954634e-01,
                              2.92843706013013e-01, 5.61536446182319e-01, 6.85137614758495e-02, 8.90853425372411e-01,
                              4.38527971314087e-01, 4.78136089625096e-01, 1.57942824868494e-01, 8.38451530279972e-01,
                              6.36273325487226e-01, 7.74039464290391e-02, 5.45646256301364e-01, 7.80075219450629e-01,
                              9.27400956530167e-01, 7.01700239834713e-02, 1.09682812589509e-01, 5.64047584357803e-01,
                              4.46922541620762e-02, 3.03438549353353e-01, 8.09200219159660e-01, 1.44245237525133e-01 };
    const double rhs_data[] = { 5.82246876454474e-01, 1.42259458199622e-02, 6.25588177982770e-01, 3.79195077388159e-01,
                                8.45448918385455e-02, 3.20808711881935e-01, 8.02701461544476e-01, 5.65141425118050e-01,
                                8.34227120735637e-01, 4.69005326248388e-01, 9.73712086117648e-01, 3.47197650692321e-01 };
    const double sol_data[] = { -2.07103028285037e-01, -3.77392587461962e-01,
                                 7.46590544020302e-01, -1.10976592416587e-01,
                                 6.07653916072011e-01,  2.05471935567110e-01 };
    gsl_matrix_complex_const_view A = gsl_matrix_complex_const_view_array(A_data, 6, 3);
    gsl_vector_complex_const_view rhs = gsl_vector_complex_const_view_array(rhs_data, 6);
    gsl_vector_complex_const_view sol = gsl_vector_complex_const_view_array(sol_data, 3);

    test_QRc_lssolve_eps(&A.matrix, &rhs.vector, &sol.vector, tol, "complex QR_lssolve test1");
  }

  return s;
}

static int
test_QRc_decomp_r_eps(const gsl_matrix_complex * m, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  gsl_matrix_complex * QR = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * T = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * A  = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * R  = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * Q  = gsl_matrix_complex_alloc(M, M);
  gsl_matrix_complex_view Q1 = gsl_matrix_complex_submatrix(Q, 0, 0, M, N);
  gsl_vector_complex_const_view tau = gsl_matrix_complex_const_diagonal(T);

  gsl_matrix_complex_memcpy(QR, m);

  s += gsl_linalg_complex_QR_decomp_r(QR, T);
  s += gsl_linalg_complex_QR_unpack_r(QR, T, Q, R);

  /* compute A = Q R */
  gsl_matrix_complex_memcpy(A, &Q1.matrix);
  gsl_blas_ztrmm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, GSL_COMPLEX_ONE, R, A);

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          gsl_complex aij = gsl_matrix_complex_get(A, i, j);
          gsl_complex mij = gsl_matrix_complex_get(m, i, j);

          gsl_test_rel(GSL_REAL(aij), GSL_REAL(mij), eps, "%s real (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i,j, GSL_REAL(aij), GSL_REAL(mij));
          gsl_test_rel(GSL_IMAG(aij), GSL_IMAG(mij), eps, "%s imag (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i,j, GSL_IMAG(aij), GSL_IMAG(mij));
        }
    }

  if (M > N)
    {
      gsl_matrix_complex * R_alt  = gsl_matrix_complex_alloc(M, N);
      gsl_matrix_complex * Q_alt  = gsl_matrix_complex_alloc(M, M);

      /* test that Q2 was computed correctly by comparing with Level 2 algorithm */
      gsl_linalg_complex_QR_unpack(QR, &tau.vector, Q_alt, R_alt);

      for (i = 0; i < M; i++)
        {
          for (j = 0; j < M; j++)
            {
              gsl_complex aij = gsl_matrix_complex_get(Q, i, j);
              gsl_complex bij = gsl_matrix_complex_get(Q_alt, i, j);

              gsl_test_rel(GSL_REAL(aij), GSL_REAL(bij), eps, "%s real Q (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, GSL_REAL(aij), GSL_REAL(bij));
              gsl_test_rel(GSL_IMAG(aij), GSL_IMAG(bij), eps, "%s imag Q (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                           desc, M, N, i, j, GSL_IMAG(aij), GSL_IMAG(bij));
            }
        }

      gsl_matrix_complex_free(R_alt);
      gsl_matrix_complex_free(Q_alt);
    }

  gsl_matrix_complex_free(QR);
  gsl_matrix_complex_free(T);
  gsl_matrix_complex_free(A);
  gsl_matrix_complex_free(Q);
  gsl_matrix_complex_free(R);

  return s;
}

static int
test_QRc_decomp_r(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 50; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix_complex * A = gsl_matrix_complex_alloc(M, N);

          create_random_complex_matrix(A, r);
          s += test_QRc_decomp_r_eps(A, 1.0e6 * M * GSL_DBL_EPSILON, "complex_QR_decomp_r random");

          gsl_matrix_complex_free(A);
        }
    }

  return s;
}

static int
test_QRc_solve_r_eps(const gsl_matrix_complex * A, const gsl_vector_complex * rhs, const gsl_vector_complex * sol,
                     const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i;

  gsl_matrix_complex * QR = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * T = gsl_matrix_complex_alloc(N, N);
  gsl_vector_complex * x  = gsl_vector_complex_alloc(N);

  gsl_matrix_complex_memcpy(QR, A);

  s += gsl_linalg_complex_QR_decomp_r(QR, T);
  s += gsl_linalg_complex_QR_solve_r(QR, T, rhs, x);

  for (i = 0; i < M; i++)
    {
      gsl_complex xi = gsl_vector_complex_get(x, i);
      gsl_complex yi = gsl_vector_complex_get(sol, i);

      gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps, "%s real (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, GSL_REAL(xi), GSL_REAL(yi));
      gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps, "%s imag (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, GSL_IMAG(xi), GSL_IMAG(yi));
    }

  gsl_matrix_complex_free(QR);
  gsl_matrix_complex_free(T);
  gsl_vector_complex_free(x);

  return s;
}

static int
test_QRc_solve_r(gsl_rng * r)
{
  int s = 0;
  size_t N;

  for (N = 1; N <= 50; ++N)
    {
      gsl_matrix_complex * A = gsl_matrix_complex_alloc(N, N);
      gsl_vector_complex * sol = gsl_vector_complex_alloc(N);
      gsl_vector_complex * rhs = gsl_vector_complex_alloc(N);

      create_random_complex_matrix(A, r);
      create_random_complex_vector(sol, r);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, A, sol, GSL_COMPLEX_ZERO, rhs);

      s += test_QRc_solve_r_eps(A, rhs, sol, 1.0e6 * N * GSL_DBL_EPSILON, "complex_QR_solve_r random");

      gsl_matrix_complex_free(A);
      gsl_vector_complex_free(sol);
      gsl_vector_complex_free(rhs);
    }

  return s;
}

static int
test_QRc_lssolve_r_eps(const gsl_matrix_complex * A, const gsl_vector_complex * rhs, const gsl_vector_complex * sol,
                       const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i;

  gsl_matrix_complex * QR = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * T = gsl_matrix_complex_alloc(N, N);
  gsl_vector_complex * x  = gsl_vector_complex_alloc(M);
  gsl_vector_complex * r  = gsl_vector_complex_alloc(M);
  gsl_vector_complex * work = gsl_vector_complex_alloc(N);
  gsl_vector_complex_const_view x1 = gsl_vector_complex_const_subvector(x, N, M - N);
  double rnorm_expected, rnorm;

  gsl_matrix_complex_memcpy(QR, A);

  s += gsl_linalg_complex_QR_decomp_r(QR, T);
  s += gsl_linalg_complex_QR_lssolve_r(QR, T, rhs, x, work);

  for (i = 0; i < N; i++)
    {
      gsl_complex xi = gsl_vector_complex_get(x, i);
      gsl_complex yi = gsl_vector_complex_get(sol, i);

      gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps, "%s sol real (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, GSL_REAL(xi), GSL_REAL(yi));
      gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps, "%s sol imag (%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, GSL_IMAG(xi), GSL_IMAG(yi));
    }

  /* compute residual and check */
  gsl_vector_complex_memcpy(r, rhs);
  gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_NEGONE, A, sol, GSL_COMPLEX_ONE, r);

  rnorm_expected = gsl_blas_dznrm2(r);
  rnorm = gsl_blas_dznrm2(&x1.vector);
  gsl_test_rel(rnorm, rnorm_expected, eps, "%s rnorm (%3lu): %22.18g   %22.18g\n",
               desc, N, rnorm, rnorm_expected);

  gsl_matrix_complex_free(QR);
  gsl_matrix_complex_free(T);
  gsl_vector_complex_free(x);
  gsl_vector_complex_free(r);
  gsl_vector_complex_free(work);

  return s;
}

static int
test_QRc_lssolve_r()
{
  int s = 0;
  const double tol = 1.0e-10;

  {
    const double A_data[] = { 6.70178096995592e-01, 2.51765675689489e-01, 7.18387884271362e-01, 3.61106008613644e-01,
                              3.36342428526987e-02, 5.25755806084206e-01, 5.87596968507183e-01, 5.27700895952418e-01,
                              6.85871226549483e-01, 3.71216051930735e-01, 4.40045953339814e-01, 2.08875308141557e-01,
                              8.62676230425306e-01, 5.79995188556082e-01, 4.86983443637474e-02, 2.82915411954634e-01,
                              2.92843706013013e-01, 5.61536446182319e-01, 6.85137614758495e-02, 8.90853425372411e-01,
                              4.38527971314087e-01, 4.78136089625096e-01, 1.57942824868494e-01, 8.38451530279972e-01,
                              6.36273325487226e-01, 7.74039464290391e-02, 5.45646256301364e-01, 7.80075219450629e-01,
                              9.27400956530167e-01, 7.01700239834713e-02, 1.09682812589509e-01, 5.64047584357803e-01,
                              4.46922541620762e-02, 3.03438549353353e-01, 8.09200219159660e-01, 1.44245237525133e-01 };
    const double rhs_data[] = { 5.82246876454474e-01, 1.42259458199622e-02, 6.25588177982770e-01, 3.79195077388159e-01,
                                8.45448918385455e-02, 3.20808711881935e-01, 8.02701461544476e-01, 5.65141425118050e-01,
                                8.34227120735637e-01, 4.69005326248388e-01, 9.73712086117648e-01, 3.47197650692321e-01 };
    const double sol_data[] = { -2.07103028285037e-01, -3.77392587461962e-01,
                                 7.46590544020302e-01, -1.10976592416587e-01,
                                 6.07653916072011e-01,  2.05471935567110e-01 };
    gsl_matrix_complex_const_view A = gsl_matrix_complex_const_view_array(A_data, 6, 3);
    gsl_vector_complex_const_view rhs = gsl_vector_complex_const_view_array(rhs_data, 6);
    gsl_vector_complex_const_view sol = gsl_vector_complex_const_view_array(sol_data, 3);

    test_QRc_lssolve_r_eps(&A.matrix, &rhs.vector, &sol.vector, tol, "complex_QR_lssolve_r test1");
  }

  return s;
}

static int
test_QRc_lssolvem_r_eps(const gsl_matrix_complex * A, const gsl_matrix_complex * rhs, const gsl_matrix_complex * sol,
                        const double eps, const char * desc)
{
  int s = 0;
  const size_t M = A->size1;
  const size_t N = A->size2;
  const size_t nrhs = rhs->size2;
  size_t k;

  gsl_matrix_complex * QR = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex * T = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * X  = gsl_matrix_complex_alloc(M, nrhs);
  gsl_vector_complex * r  = gsl_vector_complex_alloc(M);
  gsl_matrix_complex * work = gsl_matrix_complex_alloc(N, nrhs);
  double rnorm_expected, rnorm;

  gsl_matrix_complex_memcpy(QR, A);

  s += gsl_linalg_complex_QR_decomp_r(QR, T);
  s += gsl_linalg_complex_QR_lssolvem_r(QR, T, rhs, X, work);

  for (k = 0; k < nrhs; ++k)
    {
      size_t i;

      for (i = 0; i < N; i++)
        {
          gsl_complex xi = gsl_matrix_complex_get(X, i, k);
          gsl_complex yi = gsl_matrix_complex_get(sol, i, k);

          gsl_test_rel(GSL_REAL(xi), GSL_REAL(yi), eps, "%s sol real nrhs=%zu (%3lu)[%lu]: %22.18g   %22.18g\n",
                       desc, k, N, i, GSL_REAL(xi), GSL_REAL(yi));
          gsl_test_rel(GSL_IMAG(xi), GSL_IMAG(yi), eps, "%s sol imag nrhs=%zu (%3lu)[%lu]: %22.18g   %22.18g\n",
                       desc, k, N, i, GSL_IMAG(xi), GSL_IMAG(yi));
        }
    }

  if (M > N)
    {
      for (k = 0; k < nrhs; ++k)
        {
          gsl_vector_complex_const_view bk = gsl_matrix_complex_const_column(rhs, k);
          gsl_vector_complex_const_view xk = gsl_matrix_complex_const_column(sol, k);
          gsl_vector_complex_const_view xk1 = gsl_matrix_complex_const_subcolumn(X, k, N, M - N);

          /* compute residual and check */
          gsl_vector_complex_memcpy(r, &bk.vector);
          gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_NEGONE, A, &xk.vector, GSL_COMPLEX_ONE, r);

          rnorm_expected = gsl_blas_dznrm2(r);
          rnorm = gsl_blas_dznrm2(&xk1.vector);
          gsl_test_rel(rnorm, rnorm_expected, eps, "%s rnorm (%3lu): %22.18g   %22.18g\n",
                       desc, N, rnorm, rnorm_expected);
        }
    }

  gsl_matrix_complex_free(QR);
  gsl_matrix_complex_free(T);
  gsl_matrix_complex_free(X);
  gsl_vector_complex_free(r);
  gsl_matrix_complex_free(work);

  return s;
}

static int
test_QRc_lssolvem_r()
{
  int s = 0;
  const double tol = 1.0e-10;

  {
    const size_t nrhs = 5;
    const double A_data[] = { 9.04790504331649e-01, 9.80435212633848e-01, 1.07618101679494e-01, 7.79400596246746e-01,
                              7.17074458548151e-02, 9.99465052961555e-01, 3.67623687681568e-01, 7.34016608406609e-01,
                              9.99137872201521e-01, 9.99760542199001e-01, 7.67700087028891e-01, 8.29518609544701e-01,
                              4.17051704278540e-01, 2.94198087611007e-01, 6.32908348756727e-01, 1.81200841370848e-01,
                              2.68194867042750e-01, 4.97933586550018e-01, 6.56758624364592e-02, 2.57540516573773e-01,
                              8.53191450378142e-02, 4.99583529639620e-01, 2.94606770783577e-01, 5.78930337518868e-01,
                              5.30036746483061e-01, 8.76677557563410e-01, 1.29849802106692e-01, 6.95550100670606e-01,
                              4.94664396679213e-01, 7.48233978718139e-01, 5.07341136165628e-01, 3.10669190125673e-01,
                              7.96563087236516e-01, 1.84325544696376e-01, 8.39834331245121e-02, 9.90211235755048e-01 };

    const double rhs_data[] = { 4.50630210620736e-01, 4.97267370874787e-01, 2.46437678597798e-01, 3.13264224526596e-01,
                                8.59862919734670e-01, 4.10534319046098e-01, 8.48140889872288e-01, 7.97459010128614e-01,
                                1.20862599967272e-01, 7.43583140495215e-01, 4.39646575695404e-01, 5.02213530490653e-01,
                                4.45202431731400e-01, 5.29837188142220e-01, 3.77373297010637e-01, 2.04301381620252e-01,
                                7.37610234103160e-01, 3.21593171018591e-02, 2.84065205401820e-01, 6.42441908117899e-01,
                                9.21949893100580e-01, 1.24339834202199e-01, 1.39002028873496e-01, 4.93512416890204e-01,
                                2.83420745911325e-01, 8.53962751052260e-01, 6.90271409983557e-01, 5.59806616331154e-01,
                                5.64317015500352e-01, 8.71286341889171e-01, 5.60509858468075e-01, 6.88618812698512e-02,
                                4.12897987432566e-01, 5.31363799358777e-01, 1.57307252759713e-01, 2.87455856467685e-01,
                                6.64390820899717e-01, 1.11089890870597e-01, 8.99127336065928e-01, 9.40791012566710e-01,
                                4.17132391185994e-01, 6.31758887284191e-01, 9.65913520975680e-01, 1.17825148790175e-02,
                                9.16116408920661e-01, 2.10997841265471e-01, 3.97233306972833e-01, 8.77769863685340e-02,
                                3.67879563639024e-01, 5.52511777336222e-01, 6.96342946659169e-01, 7.06523500046434e-01,
                                7.08686214963970e-01, 6.91010747794824e-01, 1.80929707407004e-01, 5.19890424619853e-01,
                                8.06965146055253e-01, 1.01107591119553e-01, 3.71480780095518e-01, 2.14698063161488e-01 };

    const double sol_data[] = {  2.45487279069336e-01,  3.40446092408598e-02, -8.15927948533375e-02, -4.23937029294787e-01,
                                 5.45869724954213e-01, -4.11848713945504e-01,  3.12148846456684e-01,  2.98676719739628e-01,
                                 9.45419209111303e-02,  3.59201187278154e-02,  2.05691027187336e-01,  3.26423315586977e-01,
                                -6.54555685716176e-02,  5.46341004563229e-01, -5.63805129787749e-01,  3.03538012534692e-01,
                                 1.52247863836554e-01,  2.36910372211693e-02,  2.08544695180720e-01,  3.74036652095719e-02,
                                 2.62083323936319e-01, -5.55962448024165e-01,  7.53064079686481e-01, -4.65489181638711e-01,
                                 6.54511832291737e-01, -2.34400539603719e-01,  1.13123678116990e-01, -7.62902891573756e-01,
                                 5.25979044340358e-01, -9.37725922681537e-02 };

    gsl_matrix_complex_const_view A = gsl_matrix_complex_const_view_array(A_data, 6, 3);
    gsl_matrix_complex_const_view rhs = gsl_matrix_complex_const_view_array(rhs_data, 6, nrhs);
    gsl_matrix_complex_const_view sol = gsl_matrix_complex_const_view_array(sol_data, 3, nrhs);

    test_QRc_lssolvem_r_eps(&A.matrix, &rhs.matrix, &sol.matrix, tol, "complex_QR_lssolvem_r test1");
  }

  return s;
}
