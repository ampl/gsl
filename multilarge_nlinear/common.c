/* multilarge_nlinear/common.c
 * 
 * Copyright (C) 2015, 2016 Patrick Alken
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

static double scaled_enorm (const gsl_vector * d, const gsl_vector * f);
static void scaled_addition (const double alpha, const gsl_vector * x,
                             const double beta, const gsl_vector * y,
                             gsl_vector * z);
static double quadratic_preduction(const gsl_multilarge_nlinear_trust_state * trust_state,
                                   const gsl_vector * dx, gsl_vector * work);

/* compute || diag(d) f || */
static double
scaled_enorm (const gsl_vector * d, const gsl_vector * f)
{
  double e2 = 0;
  size_t i, n = f->size;
  for (i = 0; i < n; i++)
    {
      double fi = gsl_vector_get (f, i);
      double di = gsl_vector_get (d, i);
      double u = di * fi;
      e2 += u * u;
    }
  return sqrt (e2);
}

/* compute z = alpha*x + beta*y */
static void
scaled_addition (const double alpha, const gsl_vector * x,
                 const double beta, const gsl_vector * y, gsl_vector * z)
{
  const size_t N = z->size;
  size_t i;

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get (x, i);
      double yi = gsl_vector_get (y, i);
      gsl_vector_set (z, i, alpha * xi + beta * yi);
    }
}

/*
quadratic_preduction()
  Calculate predicted reduction based on standard
quadratic model:

m_k(dx) = Phi(x_k) + dx' g + 1/2 dx' B_k dx

predicted_reduction = m_k(0) - m_k(dx)
                    = -2 g^T dx / ||f||^2 - ( ||J*dx|| / ||f|| )^2
                    = -2 fhat . beta - ||beta||^2

where: beta = J*dx / ||f||

Inputs: trust_state - trust state
        dx          - proposed step, size p
        work        - workspace, size n

Return: predicted reduction
*/

static double
quadratic_preduction(const gsl_multilarge_nlinear_trust_state * trust_state,
                     const gsl_vector * dx, gsl_vector * work)
{
  const gsl_vector * f = trust_state->f;
  const gsl_multilarge_nlinear_parameters * params = trust_state->params;
  const double normf = gsl_blas_dnrm2(f);
  double gTdx;      /* g^T dx */
  gsl_multilarge_nlinear_fdf * fdf = trust_state->fdf;
  double pred_reduction, u;

  /* compute g^T dx */
  gsl_blas_ddot(trust_state->g, dx, &gTdx);

  /* first term: -2 g^T dx / ||f||^2 */
  pred_reduction = -2.0 * gTdx / (normf * normf);

  if (params->solver == gsl_multilarge_nlinear_solver_cholesky ||
      params->solver == gsl_multilarge_nlinear_solver_mcholesky)
    {
      const size_t p = fdf->p;
      gsl_vector_view workp = gsl_vector_subvector(work, 0, p);

      /* compute workp = J^T J dx */
      gsl_blas_dsymv(CblasLower, 1.0, trust_state->JTJ, dx, 0.0, &workp.vector);

      /* compute u = dx^T J^T J dx = ||J dx||^2 */
      gsl_blas_ddot(&workp.vector, dx, &u);

      pred_reduction -= u / (normf * normf);
    }
  else
    {
      int status;
      const gsl_vector * x = trust_state->x;
      const gsl_vector * swts = trust_state->sqrt_wts;

      /* compute work = J*dx */
      status = gsl_multilarge_nlinear_eval_df(CblasNoTrans, x, f, dx,
                                              swts, params->h_df, params->fdtype,
                                              fdf, work, NULL, NULL);
      if (status)
        {
          GSL_ERROR_VAL("error computing preduction", status, 0.0);
        }

      /* compute u = ||J*dx|| / ||f|| */
      u = gsl_blas_dnrm2(work) / normf;

      pred_reduction -= u * u;
    }

  return pred_reduction;
}
