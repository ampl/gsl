/* multifit/work.c
 * 
 * Copyright (C) 2000, 2007, 2009 Brian Gough
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>

gsl_multifit_linear_workspace *
gsl_multifit_linear_alloc (const size_t nmax, const size_t pmax)
{
  gsl_multifit_linear_workspace *w;

  w = calloc (1, sizeof (gsl_multifit_linear_workspace));

  if (w == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for multifit_linear struct",
                     GSL_ENOMEM, 0);
    }

  w->nmax = nmax;                     /* max number of observations */
  w->pmax = pmax;                     /* max number of parameters */
  w->n = 0;
  w->p = 0;
  w->rcond = 0.0;

  w->A = gsl_matrix_alloc (nmax, pmax);

  if (w->A == 0)
    {
      gsl_multifit_linear_free(w);
      GSL_ERROR_VAL ("failed to allocate space for A", GSL_ENOMEM, 0);
    }

  w->Q = gsl_matrix_alloc (pmax, pmax);

  if (w->Q == 0)
    {
      gsl_multifit_linear_free(w);
      GSL_ERROR_VAL ("failed to allocate space for Q", GSL_ENOMEM, 0);
    }

  w->QSI = gsl_matrix_alloc (pmax, pmax);

  if (w->QSI == 0)
    {
      gsl_multifit_linear_free(w);
      GSL_ERROR_VAL ("failed to allocate space for QSI", GSL_ENOMEM, 0);
    }

  w->S = gsl_vector_alloc (pmax);

  if (w->S == 0)
    {
      gsl_multifit_linear_free(w);
      GSL_ERROR_VAL ("failed to allocate space for S", GSL_ENOMEM, 0);
    }

  w->t = gsl_vector_alloc (nmax);

  if (w->t == 0)
    {
      gsl_multifit_linear_free(w);
      GSL_ERROR_VAL ("failed to allocate space for t", GSL_ENOMEM, 0);
    }

  w->xt = gsl_vector_calloc (pmax);

  if (w->xt == 0)
    {
      gsl_multifit_linear_free(w);
      GSL_ERROR_VAL ("failed to allocate space for xt", GSL_ENOMEM, 0);
    }

  w->D = gsl_vector_calloc (pmax);

  if (w->D == 0)
    {
      gsl_multifit_linear_free(w);
      GSL_ERROR_VAL ("failed to allocate space for D", GSL_ENOMEM, 0);
    }

  return w;
}

void
gsl_multifit_linear_free (gsl_multifit_linear_workspace * w)
{
  RETURN_IF_NULL (w);

  if (w->A)
    gsl_matrix_free (w->A);

  if (w->Q)
    gsl_matrix_free (w->Q);

  if (w->QSI)
    gsl_matrix_free (w->QSI);

  if (w->S)
    gsl_vector_free (w->S);

  if (w->t)
    gsl_vector_free (w->t);

  if (w->xt)
    gsl_vector_free (w->xt);

  if (w->D)
    gsl_vector_free (w->D);

  free (w);
}

