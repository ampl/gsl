/* itersolve.c
 * 
 * Copyright (C) 2014 Patrick Alken
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_splinalg.h>

gsl_splinalg_itersolve *
gsl_splinalg_itersolve_alloc(const gsl_splinalg_itersolve_type *T,
                             const size_t n, const size_t m)
{
  gsl_splinalg_itersolve *w;

  w = calloc(1, sizeof(gsl_splinalg_itersolve));
  if (w == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for itersolve struct",
                     GSL_ENOMEM);
    }

  w->type = T;
  w->normr = 0.0;

  w->state = w->type->alloc(n, m);
  if (w->state == NULL)
    {
      gsl_splinalg_itersolve_free(w);
      GSL_ERROR_NULL("failed to allocate space for itersolve state",
                     GSL_ENOMEM);
    }

  return w;
} /* gsl_splinalg_itersolve_alloc() */

void
gsl_splinalg_itersolve_free(gsl_splinalg_itersolve *w)
{
  RETURN_IF_NULL(w);

  if (w->state)
    w->type->free(w->state);

  free(w);
}

const char *
gsl_splinalg_itersolve_name(const gsl_splinalg_itersolve *w)
{
  return w->type->name;
}

int
gsl_splinalg_itersolve_iterate(const gsl_spmatrix *A, const gsl_vector *b,
                               const double tol, gsl_vector *x,
                               gsl_splinalg_itersolve *w)
{
  int status = w->type->iterate(A, b, tol, x, w->state);

  /* store current residual */
  w->normr = w->type->normr(w->state);

  return status;
}

double
gsl_splinalg_itersolve_normr(const gsl_splinalg_itersolve *w)
{
  return w->normr;
}
