#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

int
main()
{
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  const size_t N = 100;             /* number of points to interpolate */
  const double xa[] = { 0.0, 1.0 }; /* define unit square */
  const double ya[] = { 0.0, 1.0 };
  const size_t nx = sizeof(xa) / sizeof(double); /* x grid points */
  const size_t ny = sizeof(ya) / sizeof(double); /* y grid points */
  double *za = malloc(nx * ny * sizeof(double));
  gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  size_t i, j;

  /* set z grid values */
  gsl_spline2d_set(spline, za, 0, 0, 0.0);
  gsl_spline2d_set(spline, za, 0, 1, 1.0);
  gsl_spline2d_set(spline, za, 1, 1, 0.5);
  gsl_spline2d_set(spline, za, 1, 0, 1.0);

  /* initialize interpolation */
  gsl_spline2d_init(spline, xa, ya, za, nx, ny);

  /* interpolate N values in x and y and print out grid for plotting */
  for (i = 0; i < N; ++i)
    {
      double xi = i / (N - 1.0);

      for (j = 0; j < N; ++j)
        {
          double yj = j / (N - 1.0);
          double zij = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);

          printf("%f %f %f\n", xi, yj, zij);
        }
      printf("\n");
    }

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  free(za);

  return 0;
}
