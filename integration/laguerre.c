/* integration/laguerre.c
 * 
 * Copyright (C) 2017 Patrick Alken
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

/*
 * The code in this module is based on IQPACK, specifically the LGPL
 * implementation found in HERMITE_RULE:
 * https://people.sc.fsu.edu/~jburkardt/c_src/hermite_rule/hermite_rule.html
 */

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

static int
laguerre_check(const size_t n, const gsl_integration_fixed_params * params)
{
  (void) n;

  if (params->b <= 0.0)
    {
      GSL_ERROR("b must be positive", GSL_EDOM);
    }
  else if (params->alpha <= -1.0)
    {
      GSL_ERROR("alpha must be > -1", GSL_EDOM);
    }
  else
    {
      return GSL_SUCCESS;
    }
}

static int
laguerre_init(const size_t n, double * diag, double * subdiag, gsl_integration_fixed_params * params)
{
  size_t i;

  /* construct the diagonal and subdiagonal elements of Jacobi matrix */
  for (i = 0; i < n; i++)
    {
      diag[i] = 2.0 * i + 1.0 + params->alpha;
      subdiag[i] = sqrt ((i + 1.0) * (params->alpha + i + 1.0));
    }

  params->zemu = gsl_sf_gamma(params->alpha + 1.0);
  params->shft = params->a;
  params->slp = 1.0 / params->b;
  params->al = params->alpha;
  params->be = 0.0;

  return GSL_SUCCESS;
}

static const gsl_integration_fixed_type laguerre_type =
{
  laguerre_check,
  laguerre_init
};

const gsl_integration_fixed_type *gsl_integration_fixed_laguerre = &laguerre_type;
