#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

int
main (void)
{
  const size_t n = 9;                  /* number of data to interpolate */
  const double x_data[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
  const double y_data[] = { 3.0, 2.9, 2.5, 1.0, 0.9, 0.8, 0.5, 0.2, 0.1 };
  gsl_vector_const_view xv = gsl_vector_const_view_array(x_data, n);
  gsl_vector_const_view yv = gsl_vector_const_view_array(y_data, n);
  const size_t k = 4;                                 /* spline order */
  gsl_vector *c = gsl_vector_alloc(n);                /* control points for spline */
  gsl_bspline_workspace *work = gsl_bspline_alloc_ncontrol(k, n);
  gsl_matrix * XB = gsl_matrix_alloc(n, 3*(k-1) + 1); /* banded collocation matrix */
  gsl_vector_uint * ipiv = gsl_vector_uint_alloc(n);
  double t;
  size_t i;

  /* initialize knots for interpolation */
  gsl_bspline_init_interp(&xv.vector, work);

  /* compute collocation matrix for interpolation */
  gsl_bspline_col_interp(&xv.vector, XB, work);

  /* solve linear system with banded LU */
  gsl_linalg_LU_band_decomp(n, k - 1, k - 1, XB, ipiv);
  gsl_linalg_LU_band_solve(k - 1, k - 1, XB, ipiv, &yv.vector, c);

  /* output the data */
  for (i = 0; i < n; ++i)
    {
      double xi = x_data[i];
      double yi = y_data[i];
      printf("%f %f\n", xi, yi);
    }

  printf("\n\n");

  /* output the spline */
  for (t = x_data[0]; t <= x_data[n-1]; t += 0.005)
    {
      double result;
      gsl_bspline_calc(t, c, &result, work);
      printf("%f %f\n", t, result);
    }

  gsl_vector_free(c);
  gsl_bspline_free(work);
  gsl_vector_uint_free(ipiv);
  gsl_matrix_free(XB);

  return 0;
}
