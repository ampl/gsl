/* spdgemv.c
 * 
 * Copyright (C) 2012-2014 Patrick Alken
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_blas.h>

/*
gsl_spblas_dgemv()
  Multiply a sparse matrix and a vector

Inputs: alpha - scalar factor
        A     - sparse matrix
        x     - dense vector
        beta  - scalar factor
        y     - (input/output) dense vector

Return: y = alpha*op(A)*x + beta*y
*/

int
gsl_spblas_dgemv(const CBLAS_TRANSPOSE_t TransA, const double alpha,
                 const gsl_spmatrix *A, const gsl_vector *x,
                 const double beta, gsl_vector *y)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if ((TransA == CblasNoTrans && N != x->size) ||
      (TransA == CblasTrans && M != x->size))
    {
      GSL_ERROR("invalid length of x vector", GSL_EBADLEN);
    }
  else if ((TransA == CblasNoTrans && M != y->size) ||
           (TransA == CblasTrans && N != y->size))
    {
      GSL_ERROR("invalid length of y vector", GSL_EBADLEN);
    }
  else
    {
      size_t j;
      size_t incX, incY;
      size_t lenX, lenY;
      double *X, *Y;
      double *Ad;
      int *Ap, *Ai, *Aj;
      int p;

      if (TransA == CblasNoTrans)
        {
          lenX = N;
          lenY = M;
        }
      else
        {
          lenX = M;
          lenY = N;
        }

      /* form y := beta*y */

      Y = y->data;
      incY = y->stride;

      if (beta == 0.0)
        {
          size_t jy = 0;
          for (j = 0; j < lenY; ++j)
            {
              Y[jy] = 0.0;
              jy += incY;
            }
        }
      else if (beta != 1.0)
        {
          size_t jy = 0;
          for (j = 0; j < lenY; ++j)
            {
              Y[jy] *= beta;
              jy += incY;
            }
        }

      if (alpha == 0.0)
        return GSL_SUCCESS;

      /* form y := alpha*op(A)*x + y */
      Ap = A->p;
      Ad = A->data;
      X = x->data;
      incX = x->stride;

      if ((GSL_SPMATRIX_ISCCS(A) && (TransA == CblasNoTrans)) ||
          (GSL_SPMATRIX_ISCRS(A) && (TransA == CblasTrans)))
        {
          Ai = A->i;

          for (j = 0; j < lenX; ++j)
            {
              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  Y[Ai[p] * incY] += alpha * Ad[p] * X[j * incX];
                }
            }
        }
      else if ((GSL_SPMATRIX_ISCCS(A) && (TransA == CblasTrans)) ||
               (GSL_SPMATRIX_ISCRS(A) && (TransA == CblasNoTrans)))
        {
          Ai = A->i;

          for (j = 0; j < lenY; ++j)
            {
              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  Y[j * incY] += alpha * Ad[p] * X[Ai[p] * incX];
                }
            }
        }
      else if (GSL_SPMATRIX_ISTRIPLET(A))
        {
          if (TransA == CblasNoTrans)
            {
              Ai = A->i;
              Aj = A->p;
            }
          else
            {
              Ai = A->p;
              Aj = A->i;
            }

          for (p = 0; p < (int) A->nz; ++p)
            {
              Y[Ai[p] * incY] += alpha * Ad[p] * X[Aj[p] * incX];
            }
        }
      else
        {
          GSL_ERROR("unsupported matrix type", GSL_EINVAL);
        }

      return GSL_SUCCESS;
    }
} /* gsl_spblas_dgemv() */
