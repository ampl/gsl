/* integration/rational.c
 * 
 * Copyright (C) 2017 Konrad Griessinger, Patrick Alken
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
rational_check(const size_t n, const gsl_integration_fixed_params * params)
{
  if (fabs(params->b - params->a) <= GSL_DBL_EPSILON)
    {
      GSL_ERROR("|b - a| too small", GSL_EDOM);
    }
  else if (params->alpha <= -1.0)
    {
      GSL_ERROR("alpha must be > -1", GSL_EDOM);
    }
  else if (params->beta >= 0.0 || params->alpha+params->beta+2*n >= 0.0 || 0.0 >= params->alpha+2*n)
    {
      GSL_ERROR("beta < alpha + beta + 2n < 0 is required", GSL_EDOM);
    }
  else if (params->a + params->b <= 0.0)
    {
      GSL_ERROR("a + b <= 0 is not allowed", GSL_EDOM);
    }
  else
    {
      return GSL_SUCCESS;
    }
}

static int
rational_init(const size_t n, double * diag, double * subdiag, gsl_integration_fixed_params * params)
{
  const double absum = params->beta + params->alpha;
  const double a1 = params->alpha + 1.0;
  const double aba1 = absum*a1;
  double ab2i = absum + 2.0;
  size_t i;

  /* construct the diagonal and subdiagonal elements of Jacobi matrix */
  
  diag[0] = -a1/(absum + 2.0);
  subdiag[0] = sqrt( -diag[0] * ( params->beta + 1.0 ) / ( (absum + 2.0)*(absum + 3.0) ) );
  
  for (i = 1; i < n-1; i++)
    {
      ab2i += 2.0;
      diag[i] = ( -aba1 - 2.0 * i * ( absum + i + 1.0 ) ) / ( ab2i * ( ab2i - 2.0 ) );
      subdiag[i] = sqrt( (i+1.0) * ( params->alpha + i + 1.0 ) / ( ab2i - 1.0 ) * ( params->beta + i + 1.0 ) / ( ab2i * ab2i ) * ( absum + i + 1.0 ) / ( ab2i + 1.0 ) );
    }
  
  diag[n-1] = ( -aba1 - 2.0 * (n-1.0) * ( absum + n ) ) / ( (absum + 2.0*n) * ( absum + 2.0*n - 2.0 ) );
  subdiag[n-1] = 0.0;
    
  params->zemu = gsl_sf_gamma(params->alpha + 1.0) * gsl_sf_gamma(-absum - 1.0) / gsl_sf_gamma(-params->beta);
  params->shft = params->a;
  params->slp = params->b + params->a;
  params->al = params->alpha;
  params->be = params->beta;

  return GSL_SUCCESS;
}

static const gsl_integration_fixed_type rational_type =
{
  rational_check,
  rational_init
};

const gsl_integration_fixed_type *gsl_integration_fixed_rational = &rational_type;
