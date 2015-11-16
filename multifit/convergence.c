/* multifit/convergence.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

static double scaled_infnorm(const gsl_vector *x, const gsl_vector *g);

/*
gsl_multifit_fdfsolver_test()
  Convergence tests for nonlinear minimization

(1) |dx_i| <= xtol * (1 + |x_i|) for all i
(2) || g .* x ||_inf <= gtol ||f||^2
(3) ||f(x+dx) - f(x)|| <= ftol * max(||f(x)||, 1)

Inputs: s - fdfsolver
        xtol - tolerance for step size
        gtol - tolerance for gradient vector
        ftol - tolerance for residual vector
        info - (output)
          1 - stopped by small x step
          2 - stopped by small gradient
          3 - stopped by small residual vector change
*/

int
gsl_multifit_fdfsolver_test (const gsl_multifit_fdfsolver * s,
                             const double xtol, const double gtol,
                             const double ftol, int *info)
{
  int status;
  double gnorm, fnorm, phi;

  *info = 0;

  status = gsl_multifit_test_delta(s->dx, s->x, xtol*xtol, xtol);
  if (status == GSL_SUCCESS)
    {
      *info = 1;
      return GSL_SUCCESS;
    }

  /* compute gradient g = J^T f */
  (s->type->gradient) (s->state, s->g);

  /* compute gnorm = max_i( g_i * max(x_i, 1) ) */
  gnorm = scaled_infnorm(s->x, s->g);

  /* compute fnorm = ||f|| */
  fnorm = gsl_blas_dnrm2(s->f);
  phi = 0.5 * fnorm * fnorm;

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
} /* gsl_multifit_fdfsolver_test() */

int
gsl_multifit_test_delta (const gsl_vector * dx, const gsl_vector * x, 
                         double epsabs, double epsrel)
{
  size_t i;
  int ok = 1;
  const size_t n = x->size ;

  if (epsrel < 0.0)
    {
      GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL);
    }

  for (i = 0 ; i < n ; i++)
    {
      double xi = gsl_vector_get(x,i);
      double dxi = gsl_vector_get(dx,i);
      double tolerance = epsabs + epsrel * fabs(xi)  ;

      if (fabs(dxi) < tolerance)
        {
          ok = 1;
        }
      else
        {
          ok = 0;
          break;
        }
    }

  if (ok)
    return GSL_SUCCESS ;

  return GSL_CONTINUE;
}

int
gsl_multifit_test_gradient (const gsl_vector * g, double epsabs)
{
  size_t i;

  double residual = 0;

  const size_t n = g->size;

  if (epsabs < 0.0)
    {
      GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);
    }
 
  for (i = 0 ; i < n ; i++)
    {
      double gi = gsl_vector_get(g, i);
      
      residual += fabs(gi);
    }


  if (residual < epsabs)
    {
      return GSL_SUCCESS;
    }
  
  return GSL_CONTINUE ;
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
