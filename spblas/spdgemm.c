/* spdgemm.c
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
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_errno.h>

/*
gsl_spblas_dgemm()
  Multiply two sparse matrices

Inputs: alpha - scalar factor
        A     - sparse matrix
        B     - sparse matrix
        C     - (output) C = alpha * A * B

Return: success or error

Notes:
1) based on CSparse routine cs_multiply
*/

int
gsl_spblas_dgemm(const double alpha, const gsl_spmatrix *A,
                 const gsl_spmatrix *B, gsl_spmatrix *C)
{
  if (A->size2 != B->size1 || A->size1 != C->size1 || B->size2 != C->size2)
    {
      GSL_ERROR("matrix dimensions do not match", GSL_EBADLEN);
    }
  else if (A->sptype != B->sptype || A->sptype != C->sptype)
    {
      GSL_ERROR("matrix storage formats do not match", GSL_EINVAL);
    }
  else if (!GSL_SPMATRIX_ISCCS(A))
    {
      GSL_ERROR("compressed column format required", GSL_EINVAL);
    }
  else
    {
      int status = GSL_SUCCESS;
      const size_t M = A->size1;
      const size_t N = B->size2;
      int *Bi = B->i;
      int *Bp = B->p;
      double *Bd = B->data;
      int *w = A->work.work_int;      /* workspace of length M */
      double *x = C->work.work_atomic; /* workspace of length M */
      int *Cp, *Ci;
      double *Cd;
      size_t j;
      int p;
      size_t nz = 0;

      if (C->nzmax < A->nz + B->nz)
        {
          status = gsl_spmatrix_realloc(A->nz + B->nz, C);
          if (status)
            {
              GSL_ERROR("unable to realloc matrix C", status);
            }
        }

      /* initialize workspace to 0 */
      for (j = 0; j < M; ++j)
        w[j] = 0;

      Cp = C->p;
      Ci = C->i;
      Cd = C->data;

      for (j = 0; j < N; ++j)
        {
          if (nz + M > C->nzmax)
            {
              status = gsl_spmatrix_realloc(2 * C->nzmax + M, C);
              if (status)
                {
                  GSL_ERROR("unable to realloc matrix C", status);
                }

              /* these pointers could have changed due to reallocation */
              Ci = C->i;
              Cd = C->data;
            }

          Cp[j] = nz; /* column j of C starts here */

          for (p = Bp[j]; p < Bp[j + 1]; ++p)
            {
              nz = gsl_spblas_scatter(A, Bi[p], Bd[p], w, x, (int) (j + 1), C, nz);
            }

          for (p = Cp[j]; p < (int) nz; ++p)
            Cd[p] = x[Ci[p]];
        }

      Cp[N] = nz;
      C->nz = nz;

      /* scale by alpha */
      gsl_spmatrix_scale(C, alpha);

      return status;
    }
} /* gsl_spblas_dgemm() */

/*
gsl_spblas_scatter()

  Keep a running total x -> x + alpha*A(:,j) for adding matrices together in CCS,
which will eventually be stored in C(:,j)

  When a new non-zero element with row index i is found, update C->i with
the row index. C->data is updated only by the calling function after all
matrices have been added via this function.

Inputs: A     - sparse matrix m-by-n
        j     - column index
        alpha - scalar factor
        w     - keeps track which rows of column j have been added to C;
                initialize to 0 prior to first call
        x     - column vector of length m
        mark  -
        C     - output matrix whose jth column will be added to A(:,j)
        nz    - (input/output) number of non-zeros in matrix C

Notes:
1) This function is designed to be called successively when adding multiple
matrices together. Column j of C is stored contiguously as per CCS but not
necessarily in order - ie: the row indices C->i may not be in ascending order.

2) based on CSparse routine cs_scatter
*/

size_t
gsl_spblas_scatter(const gsl_spmatrix *A, const size_t j, const double alpha,
                   int *w, double *x, const int mark, gsl_spmatrix *C,
                   size_t nz)
{
  int p;
  int *Ai = A->i;
  int *Ap = A->p;
  double *Ad = A->data;
  int *Ci = C->i;

  for (p = Ap[j]; p < Ap[j + 1]; ++p)
    {
      size_t i = Ai[p];          /* A(i,j) is nonzero */

      if (w[i] < mark)           /* check if row i has been stored in column j yet */
        {
          w[i] = mark;           /* i is new entry in column j */
          Ci[nz++] = i;          /* add i to pattern of C(:,j) */
          x[i] = alpha * Ad[p];  /* x(i) = alpha * A(i,j) */
        }
      else                       /* this (i,j) exists in C from a previous call */
        {
          x[i] += alpha * Ad[p]; /* add alpha*A(i,j) to C(i,j) */
        }
    }

  return (nz) ;
} /* gsl_spblas_scatter() */
