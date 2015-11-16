/* multifit/lmmisc.c
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

/* compute step dx by solving (J^T J + mu*I) dx = -J^T f */
static int
lmniel_calc_dx(const double mu, const gsl_matrix *A, const gsl_vector *rhs,
               gsl_vector *dx, lmniel_state_t *state)
{
  int status;
  gsl_matrix *A_copy = state->A_copy;
  gsl_vector_view diag = gsl_matrix_diagonal(A_copy);

  /* make a copy of J^T J matrix */
  gsl_matrix_memcpy(A_copy, A);

  /* augment normal equations with LM parameter: A -> A + mu*I */
  gsl_vector_add_constant(&diag.vector, mu);

  status = gsl_linalg_QR_decomp(A_copy, state->work);
  if (status)
    return status;

  status = gsl_linalg_QR_solve(A_copy, state->work, rhs, dx);
  if (status)
    return status;

  return GSL_SUCCESS;
} /* lmniel_calc_dx() */

/* compute x_trial = x + dx */
static void
lmniel_trial_step(const gsl_vector * x, const gsl_vector * dx,
                  gsl_vector * x_trial)
{
  size_t i, N = x->size;

  for (i = 0; i < N; i++)
    {
      double dxi = gsl_vector_get (dx, i);
      double xi = gsl_vector_get (x, i);
      gsl_vector_set (x_trial, i, xi + dxi);
    }
} /* lmniel_trial_step() */

/*
lmniel_calc_dF()
  Compute dF = F(x) - F(x + dx) = 1/2 (f - f_new)^T (f + f_new)
*/
static double
lmniel_calc_dF(const gsl_vector *f, const gsl_vector *f_new)
{
  const size_t N = f->size;
  size_t i;
  double dF = 0.0;

  for (i = 0; i < N; ++i)
    {
      double fi = gsl_vector_get(f, i);
      double fnewi = gsl_vector_get(f_new, i);

      dF += (fi - fnewi) * (fi + fnewi);
    }

  dF *= 0.5;

  return dF;
} /* lmniel_calc_dF() */

/*
lmniel_calc_dL()
  Compute dL = L(0) - L(dx) = 1/2 dx^T (mu * D^T D dx - g)
Here, the mg input is -g
*/

static double
lmniel_calc_dL(const double mu, const gsl_vector *diag,
               const gsl_vector *dx, const gsl_vector *mg)
{
  const size_t p = dx->size;
  size_t i;
  double dL = 0.0;

  for (i = 0; i < p; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double di = gsl_vector_get(diag, i);
      double mgi = gsl_vector_get(mg, i); /* -g_i */

      dL += dxi * (mu * di * di * dxi + mgi);
    }

  dL *= 0.5;

  return dL;
} /* lmniel_calc_dL() */
