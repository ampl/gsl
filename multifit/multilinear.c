/* multifit/multilinear.c
 * 
 * Copyright (C) 2000, 2007, 2010 Brian Gough
 * Copyright (C) 2013, 2015 Patrick Alken
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

static int multifit_linear_svd (const gsl_matrix * X,
                                const int balance,
                                gsl_multifit_linear_workspace * work);

int
gsl_multifit_linear (const gsl_matrix * X,
                     const gsl_vector * y,
                     gsl_vector * c,
                     gsl_matrix * cov,
                     double *chisq, gsl_multifit_linear_workspace * work)
{
  size_t rank;
  int status = gsl_multifit_linear_tsvd(X, y, GSL_DBL_EPSILON, c, cov, chisq, &rank, work);

  return status;
}

/*
gsl_multifit_linear_tsvd()
  Solve linear least squares system with truncated SVD

Inputs: X     - least squares matrix, n-by-p
        y     - right hand side vector, n-by-1
        tol   - tolerance for singular value truncation; if
                s_j <= tol * s_0
                then it is discarded from series expansion
        c     - (output) solution vector, p-by-1
        cov   - (output) covariance matrix, p-by-p
        chisq - (output) cost function chi^2
        rank  - (output) effective rank (number of singular values
                used in solution)
        work  - workspace
*/

int
gsl_multifit_linear_tsvd (const gsl_matrix * X,
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
      double rnorm = 0.0, snorm;

      /* compute balanced SVD */
      status = gsl_multifit_linear_bsvd (X, work);
      if (status)
        return status;

      status = multifit_linear_solve (X, y, tol, -1.0, rank,
                                      c, &rnorm, &snorm, work);

      *chisq = rnorm * rnorm;

      /* variance-covariance matrix cov = s2 * (Q S^-1) (Q S^-1)^T */
      {
        double r2 = rnorm * rnorm;
        double s2 = r2 / (double)(n - *rank);
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

                gsl_matrix_set (cov, i, j, s * s2 / (d_i * d_j));
                gsl_matrix_set (cov, j, i, s * s2 / (d_i * d_j));
              }
          }
      }

      return status;
    }
}

/*
gsl_multifit_linear_svd()
  Perform SVD decomposition of the matrix X and store in work without
balancing
*/

int
gsl_multifit_linear_svd (const gsl_matrix * X,
                         gsl_multifit_linear_workspace * work)
{
  /* do not balance by default */
  int status = multifit_linear_svd(X, 0, work);

  return status;
}

/*
gsl_multifit_linear_bsvd()
  Perform SVD decomposition of the matrix X and store in work with
balancing
*/

int
gsl_multifit_linear_bsvd (const gsl_matrix * X,
                          gsl_multifit_linear_workspace * work)
{
  int status = multifit_linear_svd(X, 1, work);

  return status;
}

size_t
gsl_multifit_linear_rank(const double tol, const gsl_multifit_linear_workspace * work)
{
  double s0 = gsl_vector_get (work->S, 0);
  size_t rank = 0;
  size_t j;

  for (j = 0; j < work->p; j++)
    {
      double sj = gsl_vector_get (work->S, j);

      if (sj > tol * s0)
        ++rank;
    }

  return rank;
}

/* Estimation of values for given x */

int
gsl_multifit_linear_est (const gsl_vector * x,
                         const gsl_vector * c,
                         const gsl_matrix * cov, double *y, double *y_err)
{

  if (x->size != c->size)
    {
      GSL_ERROR ("number of parameters c does not match number of observations x",
         GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {
      GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    }
  else if (c->size != cov->size1)
    {
      GSL_ERROR ("number of parameters c does not match size of covariance matrix cov",
         GSL_EBADLEN);
    }
  else
    {
      size_t i, j;
      double var = 0;
      
      gsl_blas_ddot(x, c, y);       /* y = x.c */

      /* var = x' cov x */

      for (i = 0; i < x->size; i++)
        {
          const double xi = gsl_vector_get (x, i);
          var += xi * xi * gsl_matrix_get (cov, i, i);

          for (j = 0; j < i; j++)
            {
              const double xj = gsl_vector_get (x, j);
              var += 2 * xi * xj * gsl_matrix_get (cov, i, j);
            }
        }

      *y_err = sqrt (var);

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_rcond()
  Return reciprocal condition number of LS matrix;
gsl_multifit_linear_svd() must first be called to
compute the SVD of X and its reciprocal condition number
*/

double
gsl_multifit_linear_rcond (const gsl_multifit_linear_workspace * w)
{
  return w->rcond;
}

/*
gsl_multifit_linear_residuals()
  Compute vector of residuals from fit

Inputs: X - design matrix
        y - rhs vector
        c - fit coefficients
        r - (output) where to store residuals
*/

int
gsl_multifit_linear_residuals (const gsl_matrix *X, const gsl_vector *y,
                               const gsl_vector *c, gsl_vector *r)
{
  if (X->size1 != y->size)
    {
      GSL_ERROR
        ("number of observations in y does not match rows of matrix X",
         GSL_EBADLEN);
    }
  else if (X->size2 != c->size)
    {
      GSL_ERROR ("number of parameters c does not match columns of matrix X",
                 GSL_EBADLEN);
    }
  else if (y->size != r->size)
    {
      GSL_ERROR ("number of observations in y does not match number of residuals",
                 GSL_EBADLEN);
    }
  else
    {
      /* r = y - X c */
      gsl_vector_memcpy(r, y);
      gsl_blas_dgemv(CblasNoTrans, -1.0, X, c, 1.0, r);

      return GSL_SUCCESS;
    }
} /* gsl_multifit_linear_residuals() */

/* Perform a SVD decomposition on the least squares matrix X = U S Q^T
 *
 * Inputs: X       - least squares matrix
 *         balance - 1 to perform column balancing
 *         work    - workspace
 *
 * Notes:
 * 1) On output,
 *    work->A contains the matrix U
 *    work->Q contains the matrix Q
 *    work->S contains the vector of singular values
 * 2) The matrix X may have smaller dimensions than the workspace
 *    in the case of stdform2() - but the dimensions cannot be larger
 * 3) On output, work->n and work->p are set to the dimensions of X
 * 4) On output, work->rcond is set to the reciprocal condition number of X
 */

static int
multifit_linear_svd (const gsl_matrix * X,
                     const int balance,
                     gsl_multifit_linear_workspace * work)
{
  const size_t n = X->size1;
  const size_t p = X->size2;

  if (n > work->nmax || p > work->pmax)
    {
      GSL_ERROR("observation matrix larger than workspace", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_view A = gsl_matrix_submatrix(work->A, 0, 0, n, p);
      gsl_matrix_view Q = gsl_matrix_submatrix(work->Q, 0, 0, p, p);
      gsl_matrix_view QSI = gsl_matrix_submatrix(work->QSI, 0, 0, p, p);
      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);
      gsl_vector_view xt = gsl_vector_subvector(work->xt, 0, p);
      gsl_vector_view D = gsl_vector_subvector(work->D, 0, p);

      /* Copy X to workspace,  A <= X */

      gsl_matrix_memcpy (&A.matrix, X);

      /* Balance the columns of the matrix A if requested */

      if (balance) 
        {
          gsl_linalg_balance_columns (&A.matrix, &D.vector);
        }
      else
        {
          gsl_vector_set_all (&D.vector, 1.0);
        }

      /* decompose A into U S Q^T */
      gsl_linalg_SV_decomp_mod (&A.matrix, &QSI.matrix, &Q.matrix,
                                &S.vector, &xt.vector);

      /* compute reciprocal condition number rcond = smin / smax */
      {
        double smin, smax;
        gsl_vector_minmax(&S.vector, &smin, &smax);
        work->rcond = smin / smax;
      }

      work->n = n;
      work->p = p;

      return GSL_SUCCESS;
    }
}
