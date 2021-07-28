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
