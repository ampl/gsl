/* multilarge_nlinear/scaling.c
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

/*
 * This module handles the updating of the scaling matrix D_k in the
 * trust region subproblem:
 *
 * min m_k (dx), || D_k dx || <= Delta_k
 *
 * where m_k(dx) is a model which approximates the cost function
 * F(x_k + dx) near the current iteration point x_k
 *
 * D_k can be updated according to several different strategies.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

static int init_diag_levenberg(const gsl_matrix * JTJ, gsl_vector * diag);
static int update_diag_levenberg(const gsl_matrix * JTJ, gsl_vector * diag);

static int init_diag_marquardt(const gsl_matrix * JTJ, gsl_vector * diag);
static int update_diag_marquardt (const gsl_matrix * JTJ, gsl_vector * diag);

static int init_diag_more(const gsl_matrix * JTJ, gsl_vector * diag);
static int update_diag_more(const gsl_matrix * JTJ, gsl_vector * diag);

/* Levenberg scaling, D = I */
static int
init_diag_levenberg(const gsl_matrix * JTJ, gsl_vector * diag)
{
  (void)JTJ; /* avoid unused parameter warning */
  gsl_vector_set_all(diag, 1.0);
  return GSL_SUCCESS;
}

static int
update_diag_levenberg(const gsl_matrix * JTJ, gsl_vector * diag)
{
  (void)JTJ;  /* avoid unused parameter warning */
  (void)diag; /* avoid unused parameter warning */

  /* nothing to do */
  return GSL_SUCCESS;
}

/* initialize diagonal scaling matrix D according to Marquardt method */
static int
init_diag_marquardt(const gsl_matrix * JTJ, gsl_vector * diag)
{
  return update_diag_marquardt(JTJ, diag);
}

/* update diagonal scaling matrix D according to Marquardt method */
static int
update_diag_marquardt (const gsl_matrix * JTJ, gsl_vector * diag)
{
  const size_t p = JTJ->size2;
  size_t j;

  for (j = 0; j < p; j++)
    {
      double Jjj = gsl_matrix_get(JTJ, j, j);
      double norm;

      if (Jjj <= 0.0)
        norm = 1.0;
      else
        norm = sqrt(Jjj);

      gsl_vector_set(diag, j, norm);
    }

  return GSL_SUCCESS;
}

/* initialize diagonal scaling matrix D according to Eq 6.3 of
 * More, 1978 */
static int
init_diag_more(const gsl_matrix * JTJ, gsl_vector * diag)
{
  int status;

  gsl_vector_set_zero(diag);
  status = update_diag_more(JTJ, diag);

  return status;
}

/* update diagonal scaling matrix D according to Eq. 6.3 of
 * More, 1978 */
static int
update_diag_more (const gsl_matrix * JTJ, gsl_vector * diag)
{
  const size_t p = JTJ->size2;
  size_t j;

  for (j = 0; j < p; j++)
    {
      double Jjj = gsl_matrix_get(JTJ, j, j);
      double *diagj = gsl_vector_ptr(diag, j);
      double norm;

      if (Jjj <= 0.0)
        norm = 1.0;
      else
        norm = sqrt(Jjj);

      *diagj = GSL_MAX(*diagj, norm);
    }

  return GSL_SUCCESS;
}

static const gsl_multilarge_nlinear_scale levenberg_type =
{
  "levenberg",
  init_diag_levenberg,
  update_diag_levenberg
};

static const gsl_multilarge_nlinear_scale marquardt_type =
{
  "marquardt",
  init_diag_marquardt,
  update_diag_marquardt
};

static const gsl_multilarge_nlinear_scale more_type =
{
  "more",
  init_diag_more,
  update_diag_more
};

const gsl_multilarge_nlinear_scale *gsl_multilarge_nlinear_scale_levenberg = &levenberg_type;
const gsl_multilarge_nlinear_scale *gsl_multilarge_nlinear_scale_marquardt = &marquardt_type;
const gsl_multilarge_nlinear_scale *gsl_multilarge_nlinear_scale_more = &more_type;
