/* multifit/multiwlinear.c
 * 
 * Copyright (C) 2015 Patrick Alken
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "linear_common.c"

int
gsl_multifit_wlinear (const gsl_matrix * X,
                      const gsl_vector * w,
                      const gsl_vector * y,
                      gsl_vector * c,
                      gsl_matrix * cov,
                      double *chisq, gsl_multifit_linear_workspace * work)
{
  size_t rank;
  int status = gsl_multifit_wlinear_tsvd(X, w, y, GSL_DBL_EPSILON, c, cov, chisq, &rank, work);

  return status;
}

int
gsl_multifit_wlinear_tsvd (const gsl_matrix * X,
                           const gsl_vector * w,
                           const gsl_vector * y,
                           const double tol,
                           gsl_vector * c,
                           gsl_matrix * cov,
                           double * chisq,
                           size_t * rank,
                           gsl_multifit_linear_workspace * work)
{
  const size_t n = X->size1;
  const size_t p = X->size2;

  if (y->size != n)
    {
      GSL_ERROR("number of observations in y does not match matrix",
                GSL_EBADLEN);
    }
  else if (w->size != n)
    {
      GSL_ERROR("number of weights in w does not match matrix",
                GSL_EBADLEN);
    }
  else if (p != c->size)
    {
      GSL_ERROR ("number of parameters c does not match matrix",
                 GSL_EBADLEN);
    }
  else if (tol <= 0)
    {
      GSL_ERROR ("tolerance must be positive", GSL_EINVAL);
    }
  else
    {
      int status;
      double rnorm, snorm;
      gsl_matrix_view A = gsl_matrix_submatrix(work->A, 0, 0, n, p);
      gsl_vector_view b = gsl_vector_subvector(work->t, 0, n);

      /* compute A = sqrt(W) X, b = sqrt(W) y */
      status = gsl_multifit_linear_applyW(X, w, y, &A.matrix, &b.vector);
      if (status)
        return status;

      /* compute SVD of A */
      status = gsl_multifit_linear_bsvd(&A.matrix, work);
      if (status)
        return status;

      status = multifit_linear_solve(X, &b.vector, tol, 0.0, rank,
                                     c, &rnorm, &snorm, work);
      if (status)
        return status;

      *chisq = rnorm * rnorm;

      /* variance-covariance matrix cov = s2 * (Q S^-1) (Q S^-1)^T */
      {
        const size_t p = X->size2;
        size_t i, j;
        gsl_matrix_view QSI = gsl_matrix_submatrix(work->QSI, 0, 0, p, p);
        gsl_vector_view D = gsl_vector_subvector(work->D, 0, p);

        for (i = 0; i < p; i++)
          {
            gsl_vector_view row_i = gsl_matrix_row (&QSI.matrix, i);
            double d_i = gsl_vector_get (&D.vector, i);

            for (j = i; j < p; j++)
              {
                gsl_vector_view row_j = gsl_matrix_row (&QSI.matrix, j);
                double d_j = gsl_vector_get (&D.vector, j);
                double s;

                gsl_blas_ddot (&row_i.vector, &row_j.vector, &s);

                gsl_matrix_set (cov, i, j, s / (d_i * d_j));
                gsl_matrix_set (cov, j, i, s / (d_i * d_j));
              }
          }
      }
    }

  return GSL_SUCCESS;
}

int
gsl_multifit_wlinear_svd (const gsl_matrix * X,
                          const gsl_vector * w,
                          const gsl_vector * y,
                          double tol,
                          size_t * rank,
                          gsl_vector * c,
                          gsl_matrix * cov,
                          double *chisq, gsl_multifit_linear_workspace * work)
{
  int status = gsl_multifit_wlinear_tsvd(X, w, y, tol, c, cov, chisq, rank, work);
  return status;
}

int
gsl_multifit_wlinear_usvd (const gsl_matrix * X,
                           const gsl_vector * w,
                           const gsl_vector * y,
                           double tol,
                           size_t * rank,
                           gsl_vector * c,
                           gsl_matrix * cov,
                           double *chisq, gsl_multifit_linear_workspace * work)
{
  /* FIXME: this call does actually perform balancing, but this function is deprecated anyway */
  int status = gsl_multifit_wlinear_tsvd(X, w, y, tol, c, cov, chisq, rank, work);
  return status;
}
