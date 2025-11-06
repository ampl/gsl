/* bspline/test_ls.c
 *
 * Copyright (C) 2019, 2020 Patrick Alken
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

/*
test_LS()
  Test ordinary least squares solution on a set of data (x,y)

Inputs: order  - spline order
        nbreak - number of break points
        x      - x data
        y      - y data
        wts    - data weights
*/

static int
test_LS_eps(const size_t order, const size_t nbreak, const gsl_vector *x,
            const gsl_vector *y, const gsl_vector *wts, const double tol)
{
  int status = GSL_SUCCESS;
  gsl_bspline_workspace *w = gsl_bspline_alloc(order, nbreak);
  const size_t n = x->size;
  const size_t p = gsl_bspline_ncontrol(w);
  gsl_multifit_linear_workspace *multifit_p = gsl_multifit_linear_alloc(n, p);
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);         /* covariance matrices */
  gsl_matrix *cov_svd = gsl_matrix_alloc(p, p);
  gsl_vector *c = gsl_vector_alloc(p);              /* coefficient vectors */
  gsl_vector *c_svd = gsl_vector_alloc(p);
  gsl_vector *r = gsl_vector_alloc(n);              /* residual vectors */
  gsl_vector *r_svd = gsl_vector_alloc(n);
  double chisq, chisq_svd;
  double xmin, xmax;
  size_t i, j;

  /* uniform knots */
  gsl_vector_minmax(x, &xmin, &xmax);
  gsl_bspline_init_uniform(xmin, xmax, w);

  for (i = 0; i < n; ++i)
    {
      double xi = gsl_vector_get(x, i);
      gsl_vector_view v = gsl_matrix_row(X, i);
      gsl_bspline_eval_basis(xi, &v.vector, w);
    }

  /* solve least squares problem with banded Cholesky */
  gsl_bspline_wlssolve(x, y, wts, c, &chisq, w);

  /* solve least squares problem with dense svd */
  gsl_multifit_wlinear(X, wts, y, c_svd, cov_svd, &chisq_svd, multifit_p);

  /* test chi^2 values */
  gsl_test_rel(chisq, chisq_svd, tol,
               "least-squares: order=%zu, nbreak=%zu chisq", order, nbreak);

  /* test coefficient vectors */
  for (i = 0; i < p; ++i)
    {
      double ci = gsl_vector_get(c, i);
      double ci_expected = gsl_vector_get(c_svd, i);

      gsl_test_rel(ci, ci_expected, tol, "least-squares: order=%zu, nbreak=%zu coef[%zu]",
                   order, nbreak, i);
    }

  /* compute residual vectors and test */
  gsl_bspline_residuals(x, y, c, r, w);
  gsl_multifit_linear_residuals(X, y, c_svd, r_svd);

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);
      double ri_expected = gsl_vector_get(r_svd, i);

      gsl_test_rel(ri, ri_expected, tol, "least-squares: order=%zu, nbreak=%zu residual[%zu]",
                   order, nbreak, i);
    }

  /* compare covariance matrices */
  gsl_bspline_covariance(w->XTX, cov, w);

  for (i = 0; i < p; ++i)
    {
      for (j = 0; j <= i; ++j)
        {
          double covij = gsl_matrix_get(cov, i, j);
          double covij_expected = gsl_matrix_get(cov_svd, i, j);

          gsl_test_rel(covij, covij_expected, tol, "least-squares: order=%zu, nbreak=%zu covariance[%zu,%zu]",
                       order, nbreak, i, j);
        }
    }

  gsl_vector_free(c);
  gsl_vector_free(c_svd);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_matrix_free(cov_svd);
  gsl_vector_free(r);
  gsl_vector_free(r_svd);
  gsl_multifit_linear_free(multifit_p);
  gsl_bspline_free(w);

  return status;
}

static int
test_LS(gsl_rng * rng_p)
{
  int status = 0;
  const size_t n = 500;
  const double tol = 1.0e-10;
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *wts = gsl_vector_alloc(n);
  size_t i;

  for (i = 0; i < 10; ++i)
    {
      random_vector(-1.0, 1.0, x, rng_p);
      random_vector(-5.0, 5.0, y, rng_p);
      random_vector(0.0, 1.0, wts, rng_p);

      status += test_LS_eps(1, 7, x, y, wts, tol);
      status += test_LS_eps(2, 8, x, y, wts, tol);
      status += test_LS_eps(3, 9, x, y, wts, tol);
      status += test_LS_eps(4, 10, x, y, wts, tol);
      status += test_LS_eps(5, 11, x, y, wts, tol);
      status += test_LS_eps(6, 15, x, y, wts, tol);
    }

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(wts);

  return status;
}

/* construct LS design matrix for periodic spline fitting */
static int
test_PLS_design(const gsl_vector * x, gsl_matrix * A, gsl_bspline_workspace * w)
{
  const size_t N = x->size;
  const size_t ncontrol = gsl_bspline_ncontrol(w);
  const size_t spline_order = gsl_bspline_order(w);
  const size_t ncoef = ncontrol - spline_order + 1;
  size_t i;

  gsl_matrix_set_zero(A);

  for (i = 0; i < N; ++i)
    {
      double xi = gsl_vector_get(x, i);
      size_t jstart, j;
      size_t ncopy, nwrap;

      gsl_bspline_basis(xi, w->B, &jstart, w);

      /*
       * determine how many elements of B can be directly copied into row A(i,:) and
       * how many need to be added (wrapped) to the front of row A(i,:)
       */
      if (jstart >= ncoef)
        ncopy = 0;
      else if (jstart + spline_order < ncoef)
        ncopy = spline_order;
      else
        ncopy = ncoef - jstart;

      nwrap = spline_order - ncopy;

      for (j = 0; j < ncopy; ++j)
        {
          double Bj = gsl_vector_get(w->B, j);
          gsl_matrix_set(A, i, jstart + j, Bj);
        }

      for (j = 0; j < nwrap; ++j)
        {
          double Bj = gsl_vector_get(w->B, j + ncopy);
          double * Aij = gsl_matrix_ptr(A, i, j);
          *Aij += Bj;
        }
    }

  return GSL_SUCCESS;
}

/*
test_PLS_eps()
  Test ordinary least squares solution on a set of data (x,y)

Inputs: order  - spline order
        nbreak - number of break points
        x      - x data
        y      - y data
        wts    - data weights
        tol    - error tolerance
*/

static int
test_PLS_eps(const size_t order, const size_t nbreak, const gsl_vector *x,
             const gsl_vector *y, const gsl_vector *wts, const double tol)
{
  int status = GSL_SUCCESS;
  gsl_bspline_workspace *w = gsl_bspline_alloc(order, nbreak);
  const size_t n = x->size;
  const size_t ncontrol = gsl_bspline_ncontrol(w);
  const size_t p = ncontrol - order + 1;
  gsl_vector *x_sorted = gsl_vector_alloc(n);
  gsl_vector *y_sorted = gsl_vector_alloc(n);
  gsl_vector *wts_sorted = gsl_vector_alloc(n);
  gsl_permutation *perm = gsl_permutation_alloc(n);
  gsl_multifit_linear_workspace *multifit_p = gsl_multifit_linear_alloc(n, p);
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_matrix *cov_svd = gsl_matrix_alloc(p, p);
  gsl_vector *c = gsl_vector_alloc(ncontrol);         /* coefficient vectors */
  gsl_vector *c_svd = gsl_vector_alloc(p);
  double chisq, chisq_svd;
  double xmin, xmax;
  size_t i;

  /* sort the data vectors */
  gsl_sort_vector_index(perm, x);

  gsl_vector_memcpy(x_sorted, x);
  gsl_vector_memcpy(y_sorted, y);
  gsl_vector_memcpy(wts_sorted, wts);

  gsl_permute_vector(perm, x_sorted);
  gsl_permute_vector(perm, y_sorted);
  gsl_permute_vector(perm, wts_sorted);

  /* periodic knots */
  gsl_vector_minmax(x, &xmin, &xmax);
  gsl_bspline_init_periodic(xmin, xmax, w);

  /* solve least squares problem with periodic QR decomposition */
  gsl_bspline_pwlssolve(x_sorted, y_sorted, wts_sorted, c, &chisq, w);

  /* test matching endpoint derivatives */
  for (i = 0; i < order - 1; ++i)
    {
      double result0, result1;

      gsl_bspline_calc_deriv(xmin, c, i, &result0, w);
      gsl_bspline_calc_deriv(xmax, c, i, &result1, w);

      gsl_test_rel(result0, result1, tol,
                   "periodic-least-squares: order=%zu, nbreak=%zu deriv[%zu]",
                   order, nbreak, i);
    }

  /* solve least squares problem with dense svd */
  test_PLS_design(x, X, w);
  gsl_multifit_wlinear(X, wts, y, c_svd, cov_svd, &chisq_svd, multifit_p);

  /* test chi^2 values */
  gsl_test_rel(chisq, chisq_svd, tol,
               "periodic-least-squares: order=%zu, nbreak=%zu chisq", order, nbreak);

  /* test coefficient vectors */
  for (i = 0; i < p; ++i)
    {
      double ci = gsl_vector_get(c, i);
      double ci_expected = gsl_vector_get(c_svd, i);

      gsl_test_rel(ci, ci_expected, tol,
                   "periodic-least-squares: order=%zu, nbreak=%zu coef[%zu]",
                   order, nbreak, i);
    }

  gsl_permutation_free(perm);
  gsl_vector_free(x_sorted);
  gsl_vector_free(y_sorted);
  gsl_vector_free(wts_sorted);
  gsl_vector_free(c);
  gsl_vector_free(c_svd);
  gsl_matrix_free(X);
  gsl_matrix_free(cov_svd);
  gsl_multifit_linear_free(multifit_p);
  gsl_bspline_free(w);

  return status;
}

static int
test_PLS(gsl_rng * rng_p)
{
  int status = 0;
  const double tol = 1.0e-10;
  const size_t n = 500;
  size_t order, nbreak;
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *wts = gsl_vector_alloc(n);

  for (order = 1; order <= 10; ++order)
    {
      for (nbreak = GSL_MAX(order, 2); nbreak <= 20; ++nbreak)
        {
          random_vector(-1.0, 1.0, x, rng_p);
          random_vector(-5.0, 5.0, y, rng_p);
          random_vector(0.0, 1.0, wts, rng_p);

          status += test_PLS_eps(order, nbreak, x, y, wts, tol);
        }
    }

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(wts);

  return status;
}
