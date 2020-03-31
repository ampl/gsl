/* multifit_nlinear/test.c
 * 
 * Copyright (C) 2007, 2013, 2015, 2016 Brian Gough, Patrick Alken
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
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>

#include "test_fdf.c"

static void
test_proc(const gsl_multifit_nlinear_trs *trs,
          const gsl_multifit_nlinear_scale *scale,
          const gsl_multifit_nlinear_solver *solver,
          const int fdtype)
{
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();

  fdf_params.trs = trs;
  fdf_params.scale = scale;
  fdf_params.solver = solver;
  fdf_params.fdtype = fdtype;

  test_fdf_main(&fdf_params);
}

int
main (void)
{
  const gsl_multifit_nlinear_trs **nlinear_trs[6];
  const gsl_multifit_nlinear_solver **nlinear_solvers[5];
  const gsl_multifit_nlinear_scale **nlinear_scales[3];
  const gsl_multifit_nlinear_trs **trs;
  const gsl_multifit_nlinear_solver **solver;
  const gsl_multifit_nlinear_scale **scale;
  size_t i = 0;

  gsl_ieee_env_setup();

  /* initialize arrays */

  nlinear_trs[0] = &gsl_multifit_nlinear_trs_lm;
  nlinear_trs[1] = &gsl_multifit_nlinear_trs_lmaccel;
  nlinear_trs[2] = &gsl_multifit_nlinear_trs_dogleg;
  nlinear_trs[3] = &gsl_multifit_nlinear_trs_ddogleg;
  nlinear_trs[4] = &gsl_multifit_nlinear_trs_subspace2D;
  nlinear_trs[5] = NULL;

  nlinear_solvers[0] = &gsl_multifit_nlinear_solver_cholesky;
  nlinear_solvers[1] = &gsl_multifit_nlinear_solver_mcholesky;
  nlinear_solvers[2] = &gsl_multifit_nlinear_solver_qr;
  nlinear_solvers[3] = &gsl_multifit_nlinear_solver_svd;
  nlinear_solvers[4] = NULL;

  /* skip Marquardt scaling since it won't pass */
  nlinear_scales[0] = &gsl_multifit_nlinear_scale_levenberg;
  nlinear_scales[1] = &gsl_multifit_nlinear_scale_more;
  nlinear_scales[2] = NULL;

  /* run testsuite over all parameter combinations */

  for (trs = nlinear_trs[i]; trs != NULL; trs = nlinear_trs[++i])
    {
      size_t j = 0;

      fprintf(stderr, "trs = %s\n", (*trs)->name);

      for (solver = nlinear_solvers[j]; solver != NULL; solver = nlinear_solvers[++j])
        {
          size_t k = 0;

          /* don't use Cholesky solver with dogleg methods */
          if (i > 1 && *solver == gsl_multifit_nlinear_solver_cholesky)
            continue;

          fprintf(stderr, "solver = %s\n", (*solver)->name);
          for (scale = nlinear_scales[k]; scale != NULL; scale = nlinear_scales[++k])
            {
              test_proc(*trs, *scale, *solver, GSL_MULTIFIT_NLINEAR_FWDIFF);
              test_proc(*trs, *scale, *solver, GSL_MULTIFIT_NLINEAR_CTRDIFF);
            }
        }
    }

  exit (gsl_test_summary ());
}
