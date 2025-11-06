#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

void
print_basis(FILE *fp_knots, FILE *fp_spline, const size_t order, const size_t nbreak)
{
  gsl_bspline_workspace *w = gsl_bspline_alloc(order, nbreak);
  const size_t p = gsl_bspline_ncontrol(w);
  const size_t n = 300;
  const double a = 0.0;
  const double b = 1.0;
  const double dx = (b - a) / (n - 1.0);
  gsl_vector *B = gsl_vector_alloc(p);
  size_t i, j;

  /* use uniform breakpoints on [a, b] */
  gsl_bspline_init_uniform(a, b, w);

  gsl_vector_fprintf(fp_knots, w->knots, "%f");

  for (i = 0; i < n; ++i)
    {
      double xi = i * dx;

      gsl_bspline_eval_basis(xi, B, w);

      fprintf(fp_spline, "%f ", xi);

      for (j = 0; j < p; ++j)
        fprintf(fp_spline, "%f ", gsl_vector_get(B, j));

      fprintf(fp_spline, "\n");
    }

  fprintf(fp_knots, "\n\n");
  fprintf(fp_spline, "\n\n");

  gsl_vector_free(B);
  gsl_bspline_free(w);
}

int
main (void)
{
  const size_t nbreak = 11; /* number of breakpoints */
  FILE *fp_knots = fopen("bspline1_knots.txt", "w");
  FILE *fp_spline = fopen("bspline1.txt", "w");

  print_basis(fp_knots, fp_spline, 2, nbreak); /* linear splines */
  print_basis(fp_knots, fp_spline, 3, nbreak); /* quadratic splines */
  print_basis(fp_knots, fp_spline, 4, nbreak); /* cubic splines */
  print_basis(fp_knots, fp_spline, 5, nbreak); /* quartic splines */

  fclose(fp_knots);
  fclose(fp_spline);

  return 0;
}
