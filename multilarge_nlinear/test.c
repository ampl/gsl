/* multilarge_nlinear/test.c
 * 
 * Copyright (C) 2015, 2016 Patrick Alken
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

/* These tests are based on the NIST Statistical Reference Datasets
   See http://www.nist.gov/itl/div898/strd/index.html for more
   information. */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>

#include "test_fdf.c"

static void
test_proc(const gsl_multilarge_nlinear_trs *trs,
          const gsl_multilarge_nlinear_scale *scale,
          const int fdtype)
{
  gsl_multilarge_nlinear_parameters fdf_params =
    gsl_multilarge_nlinear_default_parameters();

  fdf_params.trs = trs;
  fdf_params.scale = scale;
  fdf_params.fdtype = fdtype;

  if (trs == gsl_multilarge_nlinear_trs_lm ||
      trs == gsl_multilarge_nlinear_trs_lmaccel)
    fdf_params.solver = gsl_multilarge_nlinear_solver_cholesky;
  else
    fdf_params.solver = gsl_multilarge_nlinear_solver_mcholesky;

  test_fdf_main(&fdf_params);
}

int
main (void)
{
  const gsl_multilarge_nlinear_trs **nlinear_trs[7];
  const gsl_multilarge_nlinear_scale **nlinear_scales[3];
  const gsl_multilarge_nlinear_trs **trs;
  const gsl_multilarge_nlinear_scale **scale;
  int fdtype;
  size_t i = 0;

  gsl_ieee_env_setup();

  /* initialize arrays */

  nlinear_trs[0] = &gsl_multilarge_nlinear_trs_lm;
  nlinear_trs[1] = &gsl_multilarge_nlinear_trs_lmaccel;
  nlinear_trs[2] = &gsl_multilarge_nlinear_trs_dogleg;
  nlinear_trs[3] = &gsl_multilarge_nlinear_trs_ddogleg;
  nlinear_trs[4] = &gsl_multilarge_nlinear_trs_subspace2D;
  nlinear_trs[5] = &gsl_multilarge_nlinear_trs_cgst;
  nlinear_trs[6] = NULL;

  nlinear_scales[0] = &gsl_multilarge_nlinear_scale_levenberg;
  nlinear_scales[1] = &gsl_multilarge_nlinear_scale_more;
  nlinear_scales[2] = NULL;

  /* run testsuite over all parameter combinations */

  for (trs = nlinear_trs[i]; trs != NULL; trs = nlinear_trs[++i])
    {
      size_t j = 0;

      fprintf(stderr, "trs = %s\n", (*trs)->name);

      for (scale = nlinear_scales[j]; scale != NULL; scale = nlinear_scales[++j])
        {
          for (fdtype = GSL_MULTILARGE_NLINEAR_FWDIFF;
               fdtype <= GSL_MULTILARGE_NLINEAR_CTRDIFF; ++fdtype)
            {
              test_proc(*trs, *scale, fdtype);
            }
        }
    }

  exit (gsl_test_summary ());
}
