/* multifit_nlinear/nielsen.c
 * 
 * Copyright (C) 2016 Patrick Alken
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
 * This module contains routines for updating the Levenberg-Marquardt
 * damping parameter on each iteration using Nielsen's method:
 *
 * [1] H. B. Nielsen, K. Madsen, Introduction to Optimization and
 *     Data Fitting, Informatics and Mathematical Modeling,
 *     Technical University of Denmark (DTU), 2010.
 *
 * 5 routines are needed to implement the update procedure:
 *
 * 1. accept - update parameter after a step has been accepted
 * 2. reject - update parameter after a step has been rejected
 * 3. free   - free workspace state
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

#define LM_ONE_THIRD         (0.333333333333333)

static int nielsen_init(const gsl_matrix * J, const gsl_vector * diag,
                        double * mu, long * nu);
static int nielsen_accept(const double rho, double * mu, long * nu);
static int nielsen_reject(double * mu, long * nu);

static int
nielsen_init(const gsl_matrix * J, const gsl_vector * diag,
             double * mu, long * nu)
{
  const double mu0 = 1.0e-3;
  const size_t p = J->size2;
  size_t j;
  double max = -1.0;

  *nu = 2;

  /* set mu = mu0 * max(diag(J~^T J~)), with J~ = J D^{-1} */

  for (j = 0; j < p; ++j)
    {
      gsl_vector_const_view v = gsl_matrix_const_column(J, j);
      double dj = gsl_vector_get(diag, j);
      double norm = gsl_blas_dnrm2(&v.vector) / dj;
      max = GSL_MAX(max, norm);
    }

  *mu = mu0 * max * max;

  return GSL_SUCCESS;
}

static int
nielsen_accept(const double rho, double * mu, long * nu)
{
  double  b;
  
  /* reset nu */
  *nu = 2;

  b = 2.0 * rho - 1.0;
  b = 1.0 - b*b*b;
  *mu *= GSL_MAX(LM_ONE_THIRD, b);

  return GSL_SUCCESS;
}

static int
nielsen_reject(double * mu, long * nu)
{
  *mu *= (double) *nu;

  /* nu := 2*nu */
  *nu <<= 1;

  return GSL_SUCCESS;
}
