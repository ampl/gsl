/* bspline/test_greville.c
 *
 * Copyright (C) 2019, 2020, 2021 Patrick Alken
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

static int
test_greville()
{
  int status = GSL_SUCCESS;

  /* Check Greville abscissae functionality on a non-uniform k=1 */
  {
    size_t i; /* looping */

    /* Test parameters */
    const size_t k = 1;
    const double bpoint_data[]    = { 0.0, 0.2, 0.5, 0.75, 1.0 };
    const size_t nbreak           = sizeof(bpoint_data)/sizeof(bpoint_data[0]);

    /* Expected results */
    const double abscissae_data[] = { 0.1, 0.35, 0.625, 0.875 };
    const size_t nabscissae       = sizeof(abscissae_data)/sizeof(abscissae_data[0]);

    gsl_vector_const_view bpoints = gsl_vector_const_view_array(bpoint_data, nbreak);
    gsl_bspline_workspace *w = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_init_augment(&bpoints.vector, w);

    gsl_test_int(nabscissae, gsl_bspline_ncontrol(w),
                 "b-spline k=%d number of abscissae", k);

    for (i = 0; i < nabscissae; ++i)
      {
        gsl_test_abs(gsl_bspline_greville_abscissa(i, w),
                     abscissae_data[i], 2*k*GSL_DBL_EPSILON,
                     "b-spline k=%d Greville abscissa #%d at x = %f",
                     k, i, abscissae_data[i]);
      }

    gsl_bspline_free(w);
  }

  /* Check Greville abscissae functionality on a non-uniform k=2 */
  {
    size_t i; /* looping */

    /* Test parameters */
    const size_t k = 2;
    const double bpoint_data[]    = { 0.0, 0.2, 0.5, 0.75, 1.0 };
    const size_t nbreak           = sizeof(bpoint_data)/sizeof(bpoint_data[0]);

    /* Expected results */
    const double abscissae_data[] = { 0.0, 0.2, 0.5, 0.75, 1.0 };
    const size_t nabscissae       = sizeof(abscissae_data)/sizeof(abscissae_data[0]);

    gsl_vector_const_view bpoints = gsl_vector_const_view_array(bpoint_data, nbreak);
    gsl_bspline_workspace *w = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_init_augment(&bpoints.vector, w);

    gsl_test_int(nabscissae, gsl_bspline_ncontrol(w),
                 "b-spline k=%d number of abscissae", k);

    for (i = 0; i < nabscissae; ++i)
      {
        gsl_test_abs(gsl_bspline_greville_abscissa(i, w),
                     abscissae_data[i], 2*k*GSL_DBL_EPSILON,
                     "b-spline k=%d Greville abscissa #%d at x = %f",
                     k, i, abscissae_data[i]);
      }

    gsl_bspline_free(w);
  }

  /* Check Greville abscissae functionality on non-uniform k=3 */
  {
    size_t i; /* looping */

    /* Test parameters */
    const size_t k = 3;
    const double bpoint_data[]    = { 0.0, 0.2, 0.5, 0.75, 1.0 };
    const size_t nbreak           = sizeof(bpoint_data)/sizeof(bpoint_data[0]);

    /* Expected results */
    const double abscissae_data[] = {      0.0, 1.0/10.0, 7.0/20.0,
                                      5.0/ 8.0, 7.0/ 8.0,      1.0 };
    const size_t nabscissae       = sizeof(abscissae_data)/sizeof(abscissae_data[0]);

    gsl_vector_const_view bpoints = gsl_vector_const_view_array(bpoint_data, nbreak);
    gsl_bspline_workspace *w = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_init_augment(&bpoints.vector, w);

    gsl_test_int(nabscissae, gsl_bspline_ncontrol(w),
                 "b-spline k=%d number of abscissae", k);

    for (i = 0; i < nabscissae; ++i)
      {
        gsl_test_abs(gsl_bspline_greville_abscissa(i, w),
                     abscissae_data[i], 2*k*GSL_DBL_EPSILON,
                     "b-spline k=%d Greville abscissa #%d at x = %f",
                     k, i, abscissae_data[i]);
      }

    gsl_bspline_free(w);
  }

  /* Check Greville abscissae functionality on non-uniform k=4 */
  {
    size_t i; /* looping */

    /* Test parameters */
    const size_t k = 4;
    const double bpoint_data[]    = { 0.0, 0.2, 0.5, 0.75, 1.0 };
    const size_t nbreak           = sizeof(bpoint_data)/sizeof(bpoint_data[0]);

    /* Expected results */
    const double abscissae_data[] = { 0.0,  1.0/15.0,  7.0/30.0,  29.0/60.0,
                                            3.0/ 4.0, 11.0/12.0,        1.0 };
    const size_t nabscissae       = sizeof(abscissae_data)/sizeof(abscissae_data[0]);

    gsl_vector_const_view bpoints = gsl_vector_const_view_array(bpoint_data, nbreak);
    gsl_bspline_workspace *w = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_init_augment(&bpoints.vector, w);

    gsl_test_int(nabscissae, gsl_bspline_ncontrol(w),
                 "b-spline k=%d number of abscissae", k);

    for (i = 0; i < nabscissae; ++i)
      {
        gsl_test_abs(gsl_bspline_greville_abscissa(i, w),
                     abscissae_data[i], 2*k*GSL_DBL_EPSILON,
                     "b-spline k=%d Greville abscissa #%d at x = %f",
                     k, i, abscissae_data[i]);
      }

    gsl_bspline_free(w);
  }

  /* Knots computed from prescribed Greville abscissae for k = 4 */
  {
    size_t i; /* looping */

    /* Test parameters */
    const size_t k = 4;
    const double abscissae_data[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
    const size_t nabscissae       = sizeof(abscissae_data)/sizeof(abscissae_data[0]);

    /* Expected results */
    const double bpoint_data[]    = { 1.0, 4.0, 4.0, 4.0, 7.0 };
    const size_t nbreak           = sizeof(bpoint_data)/sizeof(bpoint_data[0]);

    /* Compute knots from Greville abscissae */
    double abserr;
    gsl_vector_const_view abscissae
        = gsl_vector_const_view_array(abscissae_data, nabscissae);
    gsl_bspline_workspace *w = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_init_greville(&abscissae.vector, w, &abserr);

    for (i = 0; i < nbreak; ++i)
      {
        gsl_test_abs(gsl_bspline_breakpoint(i,w), bpoint_data[i], GSL_DBL_EPSILON*50,
                     "b-spline k=%d init_greville breakpoint #%d", k, i);
      }

    gsl_test_abs(abserr, 0.0, GSL_DBL_EPSILON*15,
                 "b-spline k=%d nbreak=%d init_greville abserr", k, nbreak);

    gsl_bspline_free(w);
  }

  /* Knots computed from prescribed Greville abscissae for k = 8 */
  {
    size_t i; /* looping */

    /* Test parameters */
    const size_t k = 8;

    const double abscissae_data[] = { 1.0, 10.0/7, 13.0/7, 16.0/7, 22.0/7,
                                      4.0, 34.0/7, 40.0/7, 43.0/7, 46.0/7, 7.0 };
    const size_t nabscissae       = sizeof(abscissae_data)/sizeof(abscissae_data[0]);

    /* Expected results */
    const double bpoint_data[]    = { 1.0, 4.0, 4.0, 4.0, 7.0 };
    const size_t nbreak           = sizeof(bpoint_data)/sizeof(bpoint_data[0]);

    /* Compute knots from Greville abscissae */
    double abserr;
    gsl_vector_const_view abscissae
        = gsl_vector_const_view_array(abscissae_data, nabscissae);
    gsl_bspline_workspace *w = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_init_greville(&abscissae.vector, w, &abserr);

    for (i = 0; i < nbreak; ++i)
      {
        gsl_test_abs(gsl_bspline_breakpoint(i,w), bpoint_data[i], GSL_DBL_EPSILON*50,
                     "b-spline k=%d init_greville breakpoint #%d", k, i);
      }

    gsl_test_abs(abserr, 0.0, GSL_DBL_EPSILON*15,
                 "b-spline k=%d nbreak=%d init_greville abserr", k, nbreak);

    gsl_bspline_free(w);
  }

  /* Knots computed from prescribed Greville abscissae for k = 2 */
  /* Not an interesting calculation but checks the k = 2 edge case */
  {
    size_t i; /* looping */

    /* Test parameters */
    const size_t k = 2;
    const double abscissae_data[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
    const size_t nabscissae       = sizeof(abscissae_data)/sizeof(abscissae_data[0]);

    /* Expected results */
    const double bpoint_data[]    = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
    const size_t nbreak           = sizeof(bpoint_data)/sizeof(bpoint_data[0]);

    /* Compute knots from Greville abscissae */
    double abserr;
    gsl_vector_const_view abscissae
        = gsl_vector_const_view_array(abscissae_data, nabscissae);
    gsl_bspline_workspace *w = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_init_greville(&abscissae.vector, w, &abserr);

    for (i = 0; i < nbreak; ++i)
      {
        gsl_test_abs(gsl_bspline_breakpoint(i,w), bpoint_data[i], GSL_DBL_EPSILON,
                     "b-spline k=%d init_greville breakpoint #%d", k, i);
      }

    gsl_test_abs(abserr, 0.0, GSL_DBL_EPSILON,
                 "b-spline k=%d nbreak=%d init_greville abserr", k, nbreak);

    gsl_bspline_free(w);
  }

  /* Knots computed from prescribed abscissae for edge case when nbreak = 2 */
  {
    size_t i; /* looping */

    /* Test parameters */
    const size_t k = 4;
    double abscissae_data[] = { 1.0, 3.0, 5.0, 7.0 };
    const size_t nabscissae = sizeof(abscissae_data)/sizeof(abscissae_data[0]);

    /* Expected results */
    const double bpoint_data[] = { 1.0, 7.0 };
    const size_t nbreak        = sizeof(bpoint_data)/sizeof(bpoint_data[0]);

    /* Compute knots from Greville abscissae where abscissae are recoverable */
    double abserr;
    gsl_vector_view abscissae
        = gsl_vector_view_array(abscissae_data, nabscissae);
    gsl_bspline_workspace *w = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_init_greville(&abscissae.vector, w, &abserr);

    /* Check recovery of breakpoints and abscissae */
    for (i = 0; i < nbreak; ++i)
      {
        gsl_test_abs(gsl_bspline_breakpoint(i,w), bpoint_data[i], GSL_DBL_EPSILON,
                     "b-spline k=%d init_greville breakpoint #%d", k, i);
      }

    gsl_test_abs(abserr, 0.0, GSL_DBL_EPSILON,
                 "b-spline k=%d nbreak=%d init_greville abserr", k, nbreak);

    /* Modify interior abscissae so they cannot be recovered with nbreak = 2 */
    /* Then recompute breakpoints and check that abserr is as expected */
    abscissae_data[1] -= 1;
    abscissae_data[2] += 1;
    gsl_bspline_init_greville(&abscissae.vector, w, &abserr);
    for (i = 0; i < nbreak; ++i)
      {
        gsl_test_abs(gsl_bspline_breakpoint(i,w), bpoint_data[i], GSL_DBL_EPSILON,
                     "b-spline k=%d init_greville breakpoint #%d", k, i);
      }

    gsl_test_abs(abserr, /* deliberate error */ 2.0, GSL_DBL_EPSILON,
                 "b-spline k=%d nbreak=%d init_greville abserr large", k, nbreak);

    gsl_bspline_free(w);
  }

  return status;
}
