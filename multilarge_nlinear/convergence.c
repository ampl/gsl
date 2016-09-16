/* multilargenlin/convergence.c
 * 
 * Copyright (C) 2015 Patrick Alken
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multilarge_nlinear.h>

static double scaled_infnorm(const gsl_vector *x, const gsl_vector *g);

/*
gsl_multilarge_nlinear_test()
  Convergence tests for nonlinear least squares minimization

(1) |dx_i| <= xtol * (1 + |x_i|) for all i
(2) || g .* x ||_inf <= gtol ||f||^2
(3) ||f(x+dx) - f(x)|| <= ftol * max(||f(x)||, 1)

Inputs: xtol - tolerance for step size
        gtol - tolerance for gradient vector
        ftol - tolerance for residual vector
        info - (output)
               1 - stopped by small x step
               2 - stopped by small gradient
               3 - stopped by small residual vector change
        w    - workspace
*/

int
gsl_multilarge_nlinear_test (const double xtol, const double gtol,
                             const double ftol, int *info,
                             const gsl_multilarge_nlinear_workspace * w)
{
  int status;
  double gnorm, fnorm, phi;

  *info = 0;

  status = gsl_multifit_test_delta(w->dx, w->x, xtol*xtol, xtol);
  if (status == GSL_SUCCESS)
    {
      *info = 1;
      return GSL_SUCCESS;
    }

  /* compute gnorm = max_i( g_i * max(x_i, 1) ) */
  gnorm = scaled_infnorm(w->x, w->g);

  /* compute fnorm = ||f|| */
  fnorm = gsl_blas_dnrm2(w->f);
  phi = 0.5 * fnorm * fnorm;

#if 0
  fprintf(stderr, "gnorm = %.12e fnorm = %.12e gnorm/phi = %.12e\n", gnorm, fnorm, gnorm / phi);
#endif

  if (gnorm <= gtol * GSL_MAX(phi, 1.0))
    {
      *info = 2;
      return GSL_SUCCESS;
    }

#if 0
  if (dfnorm <= ftol * GSL_MAX(fnorm, 1.0))
    {
      *info = 3;
      return GSL_SUCCESS;
    }
#endif

  return GSL_CONTINUE;
}

static double
scaled_infnorm(const gsl_vector *x, const gsl_vector *g)
{
  const size_t n = x->size;
  size_t i;
  double norm = 0.0;

  for (i = 0; i < n; ++i)
    {
      double xi = GSL_MAX(gsl_vector_get(x, i), 1.0);
      double gi = gsl_vector_get(g, i);
      double tmp = fabs(xi * gi);

      if (tmp > norm)
        norm = tmp;
    }

  return norm;
}
