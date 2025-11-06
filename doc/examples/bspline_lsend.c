#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double
f(double x)
{
  return (1.0 / (1.0 + 25.0*x*x));
}

int
main (void)
{
  const size_t n = 500;                                     /* number of data points to fit */
  const size_t k = 10;                                      /* spline order */
  const double a = -1.0;                                    /* data interval [a,b] */
  const double b = 1.0;
  const double sigma = 0.03;                                /* noise */
  const double lambda_sq = 10.0;                            /* regularization parameter */
  const size_t nderiv = 1;                                  /* derivative order to regularize at endpoints */
  gsl_bspline_workspace *work = gsl_bspline_alloc(k, 20);   /* 20 breakpoints */
  const size_t p = gsl_bspline_ncontrol(work);              /* number of control points */
  gsl_vector *c = gsl_vector_alloc(p);                      /* control points for unregularized model */
  gsl_vector *creg = gsl_vector_alloc(p);                   /* control points for regularized model */
  gsl_vector *x = gsl_vector_alloc(n);                      /* x data */
  gsl_vector *y = gsl_vector_alloc(n);                      /* y data */
  gsl_vector *w = gsl_vector_alloc(n);                      /* data weights */
  gsl_matrix *XTX = gsl_matrix_alloc(p, k);                 /* normal equations matrix */
  gsl_vector *XTy = gsl_vector_alloc(p);                    /* normal equations rhs */
  gsl_matrix *cov = gsl_matrix_alloc(p, p);                 /* covariance matrix */
  gsl_matrix *cov_reg = gsl_matrix_alloc(p, p);             /* covariance matrix for regularized model */
  gsl_matrix *A1 = gsl_matrix_alloc(p, k);                  /* left regularization matrix */
  gsl_matrix *A2 = gsl_matrix_alloc(p, k);                  /* right regularization matrix */
  size_t i;
  gsl_rng *r;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);

  /* this is the data to be fitted */
  i = 0;
  while (i < n)
    {
      double ui = gsl_rng_uniform(r);
      double xi = (b - a) * ui + a;
      double yi, dyi;

      dyi = gsl_ran_gaussian(r, sigma);
      yi = f(xi) + dyi;

      gsl_vector_set(x, i, xi);
      gsl_vector_set(y, i, yi);
      gsl_vector_set(w, i, 1.0 / (sigma * sigma));

      printf("%f %f\n", xi, yi);

      ++i;
    }

  printf("\n\n");

  /* use uniform breakpoints on [a, b] */
  gsl_bspline_init_uniform(a, b, work);

  gsl_bspline_lsnormal(x, y, w, XTy, XTX, work); /* form normal equations matrix */
  gsl_linalg_cholesky_band_decomp(XTX);          /* banded Cholesky decomposition */
  gsl_linalg_cholesky_band_solve(XTX, XTy, c);   /* solve for unregularized solution */
  gsl_bspline_covariance(XTX, cov, work);        /* compute covariance matrix */

  /* compute regularization matrices */
  gsl_bspline_oprod(nderiv, a, A1, work);
  gsl_bspline_oprod(nderiv, b, A2, work);

  /* multiply by lambda^2 */
  gsl_matrix_scale(A1, lambda_sq);
  gsl_matrix_scale(A2, lambda_sq);

  gsl_bspline_lsnormal(x, y, w, XTy, XTX, work);  /* form normal equations matrix */
  gsl_matrix_add(XTX, A1);                        /* add regularization terms */
  gsl_matrix_add(XTX, A2);
  gsl_linalg_cholesky_band_decomp(XTX);           /* banded Cholesky decomposition */
  gsl_linalg_cholesky_band_solve(XTX, XTy, creg); /* solve for regularized solution */
  gsl_bspline_covariance(XTX, cov_reg, work);     /* compute covariance matrix */

  /* output the spline curves */
  {
    double xi;

    for (xi = a; xi <= b; xi += 0.001)
      {
        double result_unreg, result_reg;
        double err_unreg, err_reg;

        /* compute unregularized spline value and error at xi */
        gsl_bspline_calc(xi, c, &result_unreg, work);
        gsl_bspline_err(xi, 0, cov, &err_unreg, work);

        /* compute regularized spline value and error at xi */
        gsl_bspline_calc(xi, creg, &result_reg, work);
        gsl_bspline_err(xi, 0, cov_reg, &err_reg, work);

        printf("%f %e %e %e %e %e\n",
               xi,
               f(xi),
               result_unreg,
               err_unreg,
               result_reg,
               err_reg);
      }
  }

  {
    double result0, result1;

    gsl_bspline_calc_deriv(a, c, nderiv, &result0, work);
    gsl_bspline_calc_deriv(b, c, nderiv, &result1, work);
    fprintf(stderr, "unregularized endpoint deriv %zu: [%14.6e, %14.6e]\n", nderiv, result0, result1);

    gsl_bspline_calc_deriv(a, creg, nderiv, &result0, work);
    gsl_bspline_calc_deriv(b, creg, nderiv, &result1, work);
    fprintf(stderr, "  regularized endpoint deriv %zu: [%14.6e, %14.6e]\n", nderiv, result0, result1);
  }

  gsl_rng_free(r);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_vector_free(creg);
  gsl_bspline_free(work);
  gsl_matrix_free(XTX);
  gsl_vector_free(XTy);
  gsl_matrix_free(cov);
  gsl_matrix_free(cov_reg);
  gsl_matrix_free(A1);
  gsl_matrix_free(A2);

  return 0;
}
