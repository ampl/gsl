/* multilarge_nlinear/dummy.c
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

/* dummy linear solver */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>

static void *dummy_alloc (const size_t n, const size_t p);
static int dummy_init(const void * vtrust_state, void * vstate);
static int dummy_presolve(const double mu, const void * vtrust_state, void * vstate);
static int dummy_solve(const gsl_vector * g, gsl_vector *x,
                       const  void * vtrust_state, void *vstate);
static int dummy_rcond(double * rcond, const gsl_matrix * JTJ, void * vstate);
static int dummy_covar(const gsl_matrix * JTJ, gsl_matrix * covar, void * vstate);

static void *
dummy_alloc (const size_t n, const size_t p)
{
  (void) n;
  (void) p;
  return NULL;
}

static void
dummy_free(void *vstate)
{
  (void) vstate;
}

static int
dummy_init(const void * vtrust_state, void * vstate)
{
  (void) vtrust_state;
  (void) vstate;
  return GSL_SUCCESS;
}

static int
dummy_presolve(const double mu, const void * vtrust_state, void * vstate)
{
  (void) mu;
  (void) vtrust_state;
  (void) vstate;
  return GSL_SUCCESS;
}

static int
dummy_solve(const gsl_vector * g, gsl_vector *x,
            const void * vtrust_state, void *vstate)
{
  (void) g;
  (void) x;
  (void) vtrust_state;
  (void) vstate;
  return GSL_SUCCESS;
}

static int
dummy_rcond(double * rcond, const gsl_matrix * JTJ, void * vstate)
{
  (void) vstate;
  (void) rcond;
  (void) JTJ;
  *rcond = 0.0;
  return GSL_SUCCESS;
}

static int
dummy_covar(const gsl_matrix * JTJ, gsl_matrix * covar, void * vstate)
{
  (void) vstate;
  (void) JTJ;
  gsl_matrix_set_zero(covar);
  return GSL_SUCCESS;
}

static const gsl_multilarge_nlinear_solver dummy_type =
{
  "dummy",
  dummy_alloc,
  dummy_init,
  dummy_presolve,
  dummy_solve,
  dummy_rcond,
  dummy_covar,
  dummy_free
};

const gsl_multilarge_nlinear_solver *gsl_multilarge_nlinear_solver_none = &dummy_type;
