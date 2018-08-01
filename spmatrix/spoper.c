/* spoper.c
 * 
 * Copyright (C) 2012 Patrick Alken
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
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

int
gsl_spmatrix_scale(gsl_spmatrix *m, const double x)
{
  size_t i;

  for (i = 0; i < m->nz; ++i)
    m->data[i] *= x;

  return GSL_SUCCESS;
} /* gsl_spmatrix_scale() */

int
gsl_spmatrix_minmax(const gsl_spmatrix *m, double *min_out, double *max_out)
{
  double min, max;
  size_t n;

  if (m->nz == 0)
    {
      GSL_ERROR("matrix is empty", GSL_EINVAL);
    }

  min = m->data[0];
  max = m->data[0];

  for (n = 1; n < m->nz; ++n)
    {
      double x = m->data[n];

      if (x < min)
        min = x;

      if (x > max)
        max = x;
    }

  *min_out = min;
  *max_out = max;

  return GSL_SUCCESS;
} /* gsl_spmatrix_minmax() */

/*
gsl_spmatrix_add()
  Add two sparse matrices

Inputs: c - (output) a + b
        a - (input) sparse matrix
        b - (input) sparse matrix

Return: success or error
*/

int
gsl_spmatrix_add(gsl_spmatrix *c, const gsl_spmatrix *a,
                 const gsl_spmatrix *b)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (b->size1 != M || b->size2 != N || c->size1 != M || c->size2 != N)
    {
      GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
    }
  else if (a->sptype != b->sptype || a->sptype != c->sptype)
    {
      GSL_ERROR("matrices must have same sparse storage format",
                GSL_EINVAL);
    }
  else if (GSL_SPMATRIX_ISTRIPLET(a))
    {
      GSL_ERROR("triplet format not yet supported", GSL_EINVAL);
    }
  else
    {
      int status = GSL_SUCCESS;
      size_t *w = a->work_sze;
      double *x = c->work_dbl;
      size_t *Cp, *Ci;
      double *Cd;
      size_t j, p;
      size_t nz = 0; /* number of non-zeros in c */
      size_t inner_size, outer_size;

      if (GSL_SPMATRIX_ISCCS(a))
        {
          inner_size = M;
          outer_size = N;
        }
      else if (GSL_SPMATRIX_ISCRS(a))
        {
          inner_size = N;
          outer_size = M;
        }
      else
        {
          GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
        }

      if (c->nzmax < a->nz + b->nz)
        {
          status = gsl_spmatrix_realloc(a->nz + b->nz, c);
          if (status)
            return status;
        }

      /* initialize w = 0 */
      for (j = 0; j < inner_size; ++j)
        w[j] = 0;

      Ci = c->i;
      Cp = c->p;
      Cd = c->data;

      for (j = 0; j < outer_size; ++j)
        {
          Cp[j] = nz;

          /* CCS: x += A(:,j); CRS: x += A(j,:) */
          nz = gsl_spblas_scatter(a, j, 1.0, w, x, j + 1, c, nz);

          /* CCS: x += B(:,j); CRS: x += B(j,:) */
          nz = gsl_spblas_scatter(b, j, 1.0, w, x, j + 1, c, nz);

          for (p = Cp[j]; p < nz; ++p)
            Cd[p] = x[Ci[p]];
        }

      /* finalize last column of c */
      Cp[j] = nz;
      c->nz = nz;

      return status;
    }
} /* gsl_spmatrix_add() */

/*
gsl_spmatrix_d2sp()
  Convert a dense gsl_matrix to sparse (triplet) format

Inputs: S - (output) sparse matrix in triplet format
        A - (input) dense matrix to convert
*/

int
gsl_spmatrix_d2sp(gsl_spmatrix *S, const gsl_matrix *A)
{
  int s = GSL_SUCCESS;
  size_t i, j;

  gsl_spmatrix_set_zero(S);
  S->size1 = A->size1;
  S->size2 = A->size2;

  for (i = 0; i < A->size1; ++i)
    {
      for (j = 0; j < A->size2; ++j)
        {
          double x = gsl_matrix_get(A, i, j);

          if (x != 0.0)
            gsl_spmatrix_set(S, i, j, x);
        }
    }

  return s;
} /* gsl_spmatrix_d2sp() */

/*
gsl_spmatrix_sp2d()
  Convert a sparse matrix to dense format
*/

int
gsl_spmatrix_sp2d(gsl_matrix *A, const gsl_spmatrix *S)
{
  if (A->size1 != S->size1 || A->size2 != S->size2)
    {
      GSL_ERROR("matrix sizes do not match", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_set_zero(A);

      if (GSL_SPMATRIX_ISTRIPLET(S))
        {
          size_t n;

          for (n = 0; n < S->nz; ++n)
            {
              size_t i = S->i[n];
              size_t j = S->p[n];
              double x = S->data[n];

              gsl_matrix_set(A, i, j, x);
            }
        }
      else
        {
          GSL_ERROR("non-triplet formats not yet supported", GSL_EINVAL);
        }

      return GSL_SUCCESS;
    }
} /* gsl_spmatrix_sp2d() */
