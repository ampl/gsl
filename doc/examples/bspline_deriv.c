#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_bspline.h>

int
main (void)
{
  const size_t nbreak = 6;
  const size_t spline_order = 4;
  gsl_bspline_workspace *w = gsl_bspline_alloc(spline_order, nbreak);
  const size_t p = gsl_bspline_ncontrol(w);
  const size_t n = 300;
  const double a = 0.0;
  const double b = 1.0;
  const double dx = (b - a) / (n - 1.0);
  gsl_matrix *dB = gsl_matrix_alloc(p, spline_order);
  size_t i, j, k;

  /* uniform breakpoints on [a, b] */
  gsl_bspline_init_uniform(a, b, w);

  /* output knot vector */
  gsl_vector_fprintf(stdout, w->knots, "%f");
  printf("\n\n");

  for (i = 0; i < spline_order; ++i)
    {
      for (j = 0; j < n; ++j)
        {
          double xj = j * dx;

          gsl_bspline_eval_deriv_basis(xj, i, dB, w);

          printf("%f ", xj);

          for (k = 0; k < p; ++k)
            printf("%f ", gsl_matrix_get(dB, k, i));

          printf("\n");
        }

      printf("\n\n");
    }

  gsl_matrix_free(dB);
  gsl_bspline_free(w);

  return 0;
}
