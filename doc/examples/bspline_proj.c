#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

/* function to be projected */
double
f(double x, void * params)
{
  (void) params;
  return (x * (-7.0 + x*(-2.0 + 3.0*x)));
}

int
main (void)
{
  const size_t k = 4;                                       /* spline order */
  const double a = -2.0;                                    /* spline interval [a,b] */
  const double b = 2.0;
  gsl_bspline_workspace *work = gsl_bspline_alloc(k, 10);   /* 10 breakpoints */
  const size_t n = gsl_bspline_ncontrol(work);              /* number of control points */
  gsl_vector *c = gsl_vector_alloc(n);                      /* control points for spline */
  gsl_vector *y = gsl_vector_alloc(n);                      /* rhs vector */
  gsl_matrix *G = gsl_matrix_alloc(n, k);                   /* Gram matrix */
  gsl_function F;

  F.function = f;
  F.params = NULL;

  /* use uniform breakpoints on [a, b] */
  gsl_bspline_init_uniform(a, b, work);

  /* compute Gram matrix */
  gsl_bspline_gram(0, G, work);

  /* construct rhs vector */
  gsl_bspline_proj_rhs(&F, y, work);

  /* solve system */
  gsl_linalg_cholesky_band_decomp(G);
  gsl_linalg_cholesky_band_solve(G, y, c);

  /* output the result */
  {
    double x;

    for (x = a; x <= b; x += 0.01)
      {
        double result;
        gsl_bspline_calc(x, c, &result, work);
        printf("%f %f %f\n", x, GSL_FN_EVAL(&F, x), result);
      }
  }

  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(G);
  gsl_bspline_free(work);

  return 0;
}
