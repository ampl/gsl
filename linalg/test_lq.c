/* linalg/test_lq.c
 *
 * Copyright (C) 2018 Patrick Alken
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

int test_LQ_solve_dim(const gsl_matrix * m, const double * actual, double eps);
int test_LQ_solve(void);
int test_LQ_LQsolve_dim(const gsl_matrix * m, const double * actual, double eps);
int test_LQ_LQsolve(void);
int test_LQ_lssolve_T_dim(const gsl_matrix * m, const double * actual, double eps);
int test_LQ_lssolve_T(void);
int test_LQ_lssolve_dim(const gsl_matrix * m, const double * actual, double eps);
int test_LQ_lssolve(void);
int test_LQ_decomp_dim(const gsl_matrix * m, double eps);
int test_LQ_decomp(void);
int test_PTLQ_solve_dim(const gsl_matrix * m, const double * actual, double eps);
int test_PTLQ_solve(void);
int test_PTLQ_LQsolve_dim(const gsl_matrix * m, const double * actual, double eps);
int test_PTLQ_LQsolve(void);
int test_PTLQ_decomp_dim(const gsl_matrix * m, double eps);
int test_PTLQ_decomp(void);
int test_LQ_update_dim(const gsl_matrix * m, double eps);
int test_LQ_update(void);

int
test_LQ_solve_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  unsigned long i, dim = m->size1;

  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * lq  = gsl_matrix_alloc(dim,dim);
  gsl_vector * d = gsl_vector_alloc(dim);
  gsl_vector * x = gsl_vector_alloc(dim);

  gsl_matrix_transpose_memcpy(lq,m);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_linalg_LQ_decomp(lq, d);
  s += gsl_linalg_LQ_solve_T(lq, d, rhs, x);
  for(i=0; i<dim; i++) {
    int foo = check(gsl_vector_get(x, i), actual[i], eps);
    if(foo) {
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(x, i), actual[i]);
    }
    s += foo;
  }

  gsl_vector_free(x);
  gsl_vector_free(d);
  gsl_matrix_free(lq);
  gsl_vector_free(rhs);

  return s;
}

int
test_LQ_solve(void)
{
  int f;
  int s = 0;

  f = test_LQ_solve_dim(hilb2, hilb2_solution, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_solve hilbert(2)");
  s += f;

  f = test_LQ_solve_dim(hilb3, hilb3_solution, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_solve hilbert(3)");
  s += f;

  f = test_LQ_solve_dim(hilb4, hilb4_solution, 4 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_solve hilbert(4)");
  s += f;

  f = test_LQ_solve_dim(hilb12, hilb12_solution, 0.5);
  gsl_test(f, "  LQ_solve hilbert(12)");
  s += f;

  f = test_LQ_solve_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_solve vander(2)");
  s += f;

  f = test_LQ_solve_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_solve vander(3)");
  s += f;

  f = test_LQ_solve_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_solve vander(4)");
  s += f;

  f = test_LQ_solve_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  LQ_solve vander(12)");
  s += f;

  return s;
}

int
test_LQ_LQsolve_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  unsigned long i, dim = m->size1;

  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * lq  = gsl_matrix_alloc(dim,dim);
  gsl_matrix * q  = gsl_matrix_alloc(dim,dim);
  gsl_matrix * l  = gsl_matrix_alloc(dim,dim);
  gsl_vector * d = gsl_vector_alloc(dim);
  gsl_vector * x = gsl_vector_alloc(dim);

  gsl_matrix_transpose_memcpy(lq,m);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_linalg_LQ_decomp(lq, d);
  s += gsl_linalg_LQ_unpack(lq, d, q, l);
  s += gsl_linalg_LQ_LQsolve(q, l, rhs, x);
  for(i=0; i<dim; i++) {
    int foo = check(gsl_vector_get(x, i), actual[i], eps);
    if(foo) {
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(x, i), actual[i]);
    }
    s += foo;
  }

  gsl_vector_free(x);
  gsl_vector_free(d);
  gsl_matrix_free(lq);
  gsl_matrix_free(q);
  gsl_matrix_free(l);
  gsl_vector_free(rhs);

  return s;
}

int
test_LQ_LQsolve(void)
{
  int f;
  int s = 0;

  f = test_LQ_LQsolve_dim(hilb2, hilb2_solution, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_LQsolve hilbert(2)");
  s += f;

  f = test_LQ_LQsolve_dim(hilb3, hilb3_solution, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_LQsolve hilbert(3)");
  s += f;

  f = test_LQ_LQsolve_dim(hilb4, hilb4_solution, 4 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_LQsolve hilbert(4)");
  s += f;

  f = test_LQ_LQsolve_dim(hilb12, hilb12_solution, 0.5);
  gsl_test(f, "  LQ_LQsolve hilbert(12)");
  s += f;

  f = test_LQ_LQsolve_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_LQsolve vander(2)");
  s += f;

  f = test_LQ_LQsolve_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_LQsolve vander(3)");
  s += f;

  f = test_LQ_LQsolve_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_LQsolve vander(4)");
  s += f;

  f = test_LQ_LQsolve_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  LQ_LQsolve vander(12)");
  s += f;

  return s;
}

int
test_LQ_lssolve_T_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  unsigned long i, M = m->size1, N = m->size2;

  gsl_vector * rhs = gsl_vector_alloc(M);
  gsl_matrix * lq  = gsl_matrix_alloc(N,M);
  gsl_vector * d = gsl_vector_alloc(N);
  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * r = gsl_vector_alloc(M);
  gsl_vector * res = gsl_vector_alloc(M);

  gsl_matrix_transpose_memcpy(lq,m);
  for(i=0; i<M; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_linalg_LQ_decomp(lq, d);
  s += gsl_linalg_LQ_lssolve_T(lq, d, rhs, x, res);

  for(i=0; i<N; i++) {
    int foo = check(gsl_vector_get(x, i), actual[i], eps);
    if(foo) {
      printf("(%3lu,%3lu)[%lu]: %22.18g   %22.18g\n", M, N, i, gsl_vector_get(x, i), actual[i]);
    }
    s += foo;
  }


   /* compute residual r = b - m x */
  if (M == N) {
    gsl_vector_set_zero(r);
  } else {
    gsl_vector_memcpy(r, rhs);
    gsl_blas_dgemv(CblasNoTrans, -1.0, m, x, 1.0, r);
  };

  for(i=0; i<N; i++) {
    int foo = check(gsl_vector_get(res, i), gsl_vector_get(r,i), sqrt(eps));
    if(foo) {
      printf("(%3lu,%3lu)[%lu]: %22.18g   %22.18g\n", M, N, i, gsl_vector_get(res, i), gsl_vector_get(r,i));
    }
    s += foo;
  }

  gsl_vector_free(r);
  gsl_vector_free(res);
  gsl_vector_free(x);
  gsl_vector_free(d);
  gsl_matrix_free(lq);
  gsl_vector_free(rhs);

  return s;
}

int
test_LQ_lssolve_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  unsigned long i, M = m->size1, N = m->size2;

  gsl_vector * rhs = gsl_vector_alloc(M);
  gsl_matrix * lq  = gsl_matrix_alloc(M, N);
  gsl_vector * tau = gsl_vector_alloc(N);
  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * r = gsl_vector_alloc(M);
  gsl_vector * res = gsl_vector_alloc(M);

  for (i = 0; i < M; i++)
    gsl_vector_set(rhs, i, i + 1.0);

  gsl_matrix_memcpy(lq, m);
  gsl_linalg_LQ_decomp(lq, tau);
  gsl_linalg_LQ_lssolve(lq, tau, rhs, x, res);

  for (i = 0; i < N; i++)
    {
      int foo = check(gsl_vector_get(x, i), actual[i], eps);
      if(foo) {
        printf("(%3lu,%3lu)[%lu]: %22.18g   %22.18g\n", M, N, i, gsl_vector_get(x, i), actual[i]);
      }
      s += foo;
    }

  /* compute residual r = b - m x */
  if (M == N) {
    gsl_vector_set_zero(r);
  } else {
    gsl_vector_memcpy(r, rhs);
    gsl_blas_dgemv(CblasNoTrans, -1.0, m, x, 1.0, r);
  };

  for (i = 0; i < N; i++)
    {
      int foo = check(gsl_vector_get(res, i), gsl_vector_get(r,i), sqrt(eps));
      if(foo) {
        printf("(%3lu,%3lu)[%lu]: %22.18g   %22.18g\n", M, N, i, gsl_vector_get(res, i), gsl_vector_get(r,i));
      }
      s += foo;
    }

  gsl_vector_free(r);
  gsl_vector_free(res);
  gsl_vector_free(x);
  gsl_vector_free(tau);
  gsl_matrix_free(lq);
  gsl_vector_free(rhs);

  return s;
}

int
test_LQ_lssolve_T(void)
{
  int f;
  int s = 0;

  f = test_LQ_lssolve_T_dim(m53, m53_lssolution, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve_T m(5,3)");
  s += f;

  f = test_LQ_lssolve_T_dim(hilb2, hilb2_solution, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve_T hilbert(2)");
  s += f;

  f = test_LQ_lssolve_T_dim(hilb3, hilb3_solution, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve_T hilbert(3)");
  s += f;

  f = test_LQ_lssolve_T_dim(hilb4, hilb4_solution, 4 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve_T hilbert(4)");
  s += f;

  f = test_LQ_lssolve_T_dim(hilb12, hilb12_solution, 0.5);
  gsl_test(f, "  LQ_lssolve_T hilbert(12)");
  s += f;

  f = test_LQ_lssolve_T_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve_T vander(2)");
  s += f;

  f = test_LQ_lssolve_T_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve_T vander(3)");
  s += f;

  f = test_LQ_lssolve_T_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve_T vander(4)");
  s += f;

  f = test_LQ_lssolve_T_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  LQ_lssolve_T vander(12)");
  s += f;

  return s;
}

int
test_LQ_lssolve(void)
{
  int f;
  int s = 0;

  f = test_LQ_lssolve_dim(hilb2, hilb2_solution, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve hilbert(2)");
  s += f;

  f = test_LQ_lssolve_dim(hilb3, hilb3_solution, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve hilbert(3)");
  s += f;

  f = test_LQ_lssolve_dim(hilb4, hilb4_solution, 4 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve hilbert(4)");
  s += f;

  f = test_LQ_lssolve_dim(hilb12, hilb12_solution, 0.5);
  gsl_test(f, "  LQ_lssolve hilbert(12)");
  s += f;

  f = test_LQ_lssolve_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve vander(2)");
  s += f;

  f = test_LQ_lssolve_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve vander(3)");
  s += f;

  f = test_LQ_lssolve_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_lssolve vander(4)");
  s += f;

  f = test_LQ_lssolve_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  LQ_lssolve vander(12)");
  s += f;

  return s;
}

int
test_LQ_decomp_dim(const gsl_matrix * m, double eps)
{
  int s = 0;
  unsigned long i,j, M = m->size1, N = m->size2;

  gsl_matrix * lq = gsl_matrix_alloc(M,N);
  gsl_matrix * a  = gsl_matrix_alloc(M,N);
  gsl_matrix * q  = gsl_matrix_alloc(N,N);
  gsl_matrix * l  = gsl_matrix_alloc(M,N);
  gsl_vector * d = gsl_vector_alloc(GSL_MIN(M,N));

  gsl_matrix_memcpy(lq,m);

  s += gsl_linalg_LQ_decomp(lq, d);
  s += gsl_linalg_LQ_unpack(lq, d, q, l);
  
   /* compute a = q r */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, l, q, 0.0, a);

  for(i=0; i<M; i++) {
    for(j=0; j<N; j++) {
      double aij = gsl_matrix_get(a, i, j);
      double mij = gsl_matrix_get(m, i, j);
      int foo = check(aij, mij, eps);
      if(foo) {
        printf("(%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n", M, N, i,j, aij, mij);
      }
      s += foo;
    }
  }

  gsl_vector_free(d);
  gsl_matrix_free(lq);
  gsl_matrix_free(a);
  gsl_matrix_free(q);
  gsl_matrix_free(l);

  return s;
}

int
test_LQ_decomp(void)
{
  int f;
  int s = 0;

  f = test_LQ_decomp_dim(m35, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp m(3,5)");
  s += f;

  f = test_LQ_decomp_dim(m53, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp m(5,3)");
  s += f;

  f = test_LQ_decomp_dim(hilb2, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp hilbert(2)");
  s += f;

  f = test_LQ_decomp_dim(hilb3, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp hilbert(3)");
  s += f;

  f = test_LQ_decomp_dim(hilb4, 4 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp hilbert(4)");
  s += f;

  f = test_LQ_decomp_dim(hilb12, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp hilbert(12)");
  s += f;

  f = test_LQ_decomp_dim(vander2, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp vander(2)");
  s += f;

  f = test_LQ_decomp_dim(vander3, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp vander(3)");
  s += f;

  f = test_LQ_decomp_dim(vander4, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_decomp vander(4)");
  s += f;

  f = test_LQ_decomp_dim(vander12, 0.0005);  /* FIXME: bad accuracy */
  gsl_test(f, "  LQ_decomp vander(12)");
  s += f;

  return s;
}

int
test_PTLQ_solve_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  int signum;
  unsigned long i, dim = m->size1;

  gsl_permutation * perm = gsl_permutation_alloc(dim);
  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * lq  = gsl_matrix_alloc(dim,dim);
  gsl_vector * d = gsl_vector_alloc(dim);
  gsl_vector * x = gsl_vector_alloc(dim);
  gsl_vector * norm = gsl_vector_alloc(dim);

  gsl_matrix_transpose_memcpy(lq,m);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_linalg_PTLQ_decomp(lq, d, perm, &signum, norm);
  s += gsl_linalg_PTLQ_solve_T(lq, d, perm, rhs, x);
  for(i=0; i<dim; i++) {
    int foo = check(gsl_vector_get(x, i), actual[i], eps);
    if(foo) {
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(x, i), actual[i]);
    }
    s += foo;
  }

  gsl_vector_free(norm);
  gsl_vector_free(x);
  gsl_vector_free(d);
  gsl_matrix_free(lq);
  gsl_vector_free(rhs);
  gsl_permutation_free(perm);

  return s;
}

int test_PTLQ_solve(void)
{
  int f;
  int s = 0;

  f = test_PTLQ_solve_dim(hilb2, hilb2_solution, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_solve hilbert(2)");
  s += f;

  f = test_PTLQ_solve_dim(hilb3, hilb3_solution, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_solve hilbert(3)");
  s += f;

  f = test_PTLQ_solve_dim(hilb4, hilb4_solution, 2 * 2048.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_solve hilbert(4)");
  s += f;

  f = test_PTLQ_solve_dim(hilb12, hilb12_solution, 0.5);
  gsl_test(f, "  PTLQ_solve hilbert(12)");
  s += f;

  f = test_PTLQ_solve_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_solve vander(2)");
  s += f;

  f = test_PTLQ_solve_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_solve vander(3)");
  s += f;

  f = test_PTLQ_solve_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_solve vander(4)");
  s += f;

  f = test_PTLQ_solve_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  PTLQ_solve vander(12)");
  s += f;

  return s;
}

int
test_PTLQ_LQsolve_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  int signum;
  unsigned long i, dim = m->size1;

  gsl_permutation * perm = gsl_permutation_alloc(dim);
  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * lq  = gsl_matrix_alloc(dim,dim);
  gsl_matrix * q  = gsl_matrix_alloc(dim,dim);
  gsl_matrix * l  = gsl_matrix_alloc(dim,dim);
  gsl_vector * d = gsl_vector_alloc(dim);
  gsl_vector * x = gsl_vector_alloc(dim);
  gsl_vector * norm = gsl_vector_alloc(dim);

  gsl_matrix_transpose_memcpy(lq,m);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_linalg_PTLQ_decomp2(lq, q, l, d, perm, &signum, norm);
  s += gsl_linalg_PTLQ_LQsolve_T(q, l, perm, rhs, x);
  for(i=0; i<dim; i++) {
    int foo = check(gsl_vector_get(x, i), actual[i], eps);
    if(foo) {
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(x, i), actual[i]);
    }
    s += foo;
  }

  gsl_vector_free(norm);
  gsl_vector_free(x);
  gsl_vector_free(d);
  gsl_matrix_free(lq);
  gsl_vector_free(rhs);
  gsl_permutation_free(perm);

  return s;
}

int test_PTLQ_LQsolve(void)
{
  int f;
  int s = 0;

  f = test_PTLQ_LQsolve_dim(hilb2, hilb2_solution, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_LQsolve hilbert(2)");
  s += f;

  f = test_PTLQ_LQsolve_dim(hilb3, hilb3_solution, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_LQsolve hilbert(3)");
  s += f;

  f = test_PTLQ_LQsolve_dim(hilb4, hilb4_solution, 2 * 2048.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_LQsolve hilbert(4)");
  s += f;

  f = test_PTLQ_LQsolve_dim(hilb12, hilb12_solution, 0.5);
  gsl_test(f, "  PTLQ_LQsolve hilbert(12)");
  s += f;

  f = test_PTLQ_LQsolve_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_LQsolve vander(2)");
  s += f;

  f = test_PTLQ_LQsolve_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_LQsolve vander(3)");
  s += f;

  f = test_PTLQ_LQsolve_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_LQsolve vander(4)");
  s += f;

  f = test_PTLQ_LQsolve_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  PTLQ_LQsolve vander(12)");
  s += f;

  return s;
}

int
test_PTLQ_decomp_dim(const gsl_matrix * m, double eps)
{
  int s = 0, signum;
  unsigned long i,j, M = m->size1, N = m->size2;

  gsl_matrix * lq = gsl_matrix_alloc(N,M);
  gsl_matrix * a  = gsl_matrix_alloc(N,M);
  gsl_matrix * q  = gsl_matrix_alloc(M,M);
  gsl_matrix * l  = gsl_matrix_alloc(N,M);
  gsl_vector * d = gsl_vector_alloc(GSL_MIN(M,N));
  gsl_vector * norm = gsl_vector_alloc(N);

  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_matrix_transpose_memcpy(lq,m);

  s += gsl_linalg_PTLQ_decomp(lq, d, perm, &signum, norm);
  s += gsl_linalg_LQ_unpack(lq, d, q, l);

   /* compute a = l q */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, l, q, 0.0, a);


   /* Compute P LQ  by permuting the rows of LQ */

  for (i = 0; i < M; i++) {
    gsl_vector_view col = gsl_matrix_column (a, i);
    gsl_permute_vector_inverse (perm, &col.vector);
  }

  for(i=0; i<M; i++) {
    for(j=0; j<N; j++) {
      double aij = gsl_matrix_get(a, j, i);
      double mij = gsl_matrix_get(m, i, j);
      int foo = check(aij, mij, eps);
      if(foo) {
        printf("(%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n", M, N, i,j, aij, mij);
      }
      s += foo;
    }
  }

  gsl_permutation_free (perm);
  gsl_vector_free(norm);
  gsl_vector_free(d);
  gsl_matrix_free(lq);
  gsl_matrix_free(a);
  gsl_matrix_free(q);
  gsl_matrix_free(l);

  return s;
}

int test_PTLQ_decomp(void)
{
  int f;
  int s = 0;

  f = test_PTLQ_decomp_dim(m35, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp m(3,5)");
  s += f;

  f = test_PTLQ_decomp_dim(m53, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp m(5,3)");
  s += f;

  f = test_PTLQ_decomp_dim(s35, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp s(3,5)");
  s += f;

  f = test_PTLQ_decomp_dim(s53, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp s(5,3)");
  s += f;

  f = test_PTLQ_decomp_dim(hilb2, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp hilbert(2)");
  s += f;

  f = test_PTLQ_decomp_dim(hilb3, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp hilbert(3)");
  s += f;

  f = test_PTLQ_decomp_dim(hilb4, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp hilbert(4)");
  s += f;

  f = test_PTLQ_decomp_dim(hilb12, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp hilbert(12)");
  s += f;

  f = test_PTLQ_decomp_dim(vander2, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp vander(2)");
  s += f;

  f = test_PTLQ_decomp_dim(vander3, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp vander(3)");
  s += f;

  f = test_PTLQ_decomp_dim(vander4, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  PTLQ_decomp vander(4)");
  s += f;

  f = test_PTLQ_decomp_dim(vander12, 0.0005);  /* FIXME: bad accuracy */
  gsl_test(f, "  PTLQ_decomp vander(12)");
  s += f;

  return s;
}

int
test_LQ_update_dim(const gsl_matrix * m, double eps)
{
  int s = 0;
  unsigned long i,j, M = m->size1, N = m->size2;

  gsl_matrix * lq1  = gsl_matrix_alloc(N,M);
  gsl_matrix * lq2  = gsl_matrix_alloc(N,M);
  gsl_matrix * q1  = gsl_matrix_alloc(M,M);
  gsl_matrix * l1  = gsl_matrix_alloc(N,M);
  gsl_matrix * q2  = gsl_matrix_alloc(M,M);
  gsl_matrix * l2  = gsl_matrix_alloc(N,M);
  gsl_vector * d2 = gsl_vector_alloc(GSL_MIN(M,N));
  gsl_vector * u = gsl_vector_alloc(M);
  gsl_vector * v = gsl_vector_alloc(N);
  gsl_vector * w = gsl_vector_alloc(M);

  gsl_matrix_transpose_memcpy(lq1,m);
  gsl_matrix_transpose_memcpy(lq2,m);
  for(i=0; i<M; i++) gsl_vector_set(u, i, sin(i+1.0));
  for(i=0; i<N; i++) gsl_vector_set(v, i, cos(i+2.0) + sin(i*i+3.0));

  /* lq1 is updated */

  gsl_blas_dger(1.0, v, u, lq1);

  /* lq2 is first decomposed, updated later */

  s += gsl_linalg_LQ_decomp(lq2, d2);
  s += gsl_linalg_LQ_unpack(lq2, d2, q2, l2);

  /* compute w = Q^T u */

  gsl_blas_dgemv(CblasNoTrans, 1.0, q2, u, 0.0, w);

  /* now lq2 is updated */

  s += gsl_linalg_LQ_update(q2, l2, v, w);

  /* multiply q2*l2 */

  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,l2,q2,0.0,lq2);

  /*  check lq1==lq2 */

  for(i=0; i<N; i++) {
    for(j=0; j<M; j++) {
      double s1 = gsl_matrix_get(lq1, i, j);
      double s2 = gsl_matrix_get(lq2, i, j);
      
      int foo = check(s1, s2, eps);
#if 0
      if(foo) {
	  printf("LQ:(%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n", M, N, i,j, s1, s2);
      }
#endif
      s += foo;
    }
  }

  gsl_vector_free(d2);
  gsl_vector_free(u);
  gsl_vector_free(v);
  gsl_vector_free(w);
  gsl_matrix_free(lq1);
  gsl_matrix_free(lq2);
  gsl_matrix_free(q1);
  gsl_matrix_free(l1);
  gsl_matrix_free(q2);
  gsl_matrix_free(l2);

  return s;
}

int
test_LQ_update(void)
{
  int f;
  int s = 0;

  f = test_LQ_update_dim(m35, 2 * 512.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update m(3,5)");
  s += f;

  f = test_LQ_update_dim(m53, 2 * 2048.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update m(5,3)");
  s += f;

  f = test_LQ_update_dim(hilb2,  2 * 512.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update hilbert(2)");
  s += f;

  f = test_LQ_update_dim(hilb3,  2 * 512.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update hilbert(3)");
  s += f;

  f = test_LQ_update_dim(hilb4, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update hilbert(4)");
  s += f;

  f = test_LQ_update_dim(hilb12, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update hilbert(12)");
  s += f;

  f = test_LQ_update_dim(vander2, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update vander(2)");
  s += f;

  f = test_LQ_update_dim(vander3, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update vander(3)");
  s += f;

  f = test_LQ_update_dim(vander4, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  LQ_update vander(4)");
  s += f;

  f = test_LQ_update_dim(vander12, 0.0005);  /* FIXME: bad accuracy */
  gsl_test(f, "  LQ_update vander(12)");
  s += f;

  return s;
}
