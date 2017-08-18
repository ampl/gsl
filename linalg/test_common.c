/* linalg/test_common.c
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

static int create_random_vector(gsl_vector * v, gsl_rng * r);
static int create_posdef_matrix(gsl_matrix * m, gsl_rng * r);
static int create_hilbert_matrix2(gsl_matrix * m);

static int
create_random_vector(gsl_vector * v, gsl_rng * r)
{
  const size_t N = v->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double vi = gsl_rng_uniform(r);
      gsl_vector_set(v, i, vi);
    }

  return GSL_SUCCESS;
}

static int
create_symm_matrix(gsl_matrix * m, gsl_rng * r)
{
  const size_t N = m->size1;
  size_t i, j;

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j <= i; ++j)
        {
          double mij = gsl_rng_uniform(r);
          gsl_matrix_set(m, i, j, mij);
        }
    }

  /* copy lower triangle to upper */
  gsl_matrix_transpose_tricpy('L', 0, m, m);

  return GSL_SUCCESS;
}

static int
create_posdef_matrix(gsl_matrix * m, gsl_rng * r)
{
  const size_t N = m->size1;
  const double alpha = 10.0 * N;
  size_t i;

  /* The idea is to make a symmetric diagonally dominant
   * matrix. Make a symmetric matrix and add alpha*I to
   * its diagonal */

  create_symm_matrix(m, r);

  for (i = 0; i < N; ++i)
    {
      double mii = gsl_matrix_get(m, i, i);
      gsl_matrix_set(m, i, i, mii + alpha);
    }

  return GSL_SUCCESS;
}

static int
create_hilbert_matrix2(gsl_matrix * m)
{
  const size_t N = m->size1;
  size_t i, j;

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          gsl_matrix_set(m, i, j, 1.0/(i+j+1.0));
        }
    }

  return GSL_SUCCESS;
}
