/* randist/wishart.c
 *
 * Copyright (C) 2017 Timothée Flutre
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>


/* Generate a random matrix from a Wishart distribution using the Bartlett
 * decomposition, following Smith and Hocking, Journal of the Royal Statistical
 * Society. Series C (Applied Statistics), Vol. 21, No. 3 (1972), pp. 341-345.
 *
 * df      degrees of freedom
 * L       matrix resulting from the Cholesky decomposition of
 *         the scale matrix V = L L^T (dimension d x d)
 * result  output matrix (dimension d x d)
 * work    matrix used for intermediate computations (dimension d x d)
 */
int
gsl_ran_wishart (const gsl_rng * r,
                 const double df,
                 const gsl_matrix * L,
                 gsl_matrix * result,
                 gsl_matrix * work)
{
  if (L->size1 != L->size2)
    {
      GSL_ERROR("L should be a square matrix", GSL_ENOTSQR);
    }
  else if (result->size1 != result->size2)
    {
      GSL_ERROR("result should be a square matrix", GSL_ENOTSQR);
    }
  else if (work->size1 != work->size2)
    {
      GSL_ERROR("work should be a square matrix", GSL_ENOTSQR);
    }
  else if (result->size1 != L->size1)
    {
      GSL_ERROR("incompatible dimensions of result matrix", GSL_EBADLEN);
    }
  else if (work->size1 != L->size1)
    {
      GSL_ERROR("incompatible dimensions of work matrix", GSL_EBADLEN);
    }
  else if (df <= L->size1 - 1)
    {
      GSL_ERROR("incompatible degrees of freedom", GSL_EDOM);
    }
  else
    {
      /* result: X = L A A^T L^T */

      size_t d = L->size1, i, j;

      /* insure the upper part of A is zero before filling its lower part */
      gsl_matrix_set_zero(work);
      for (i = 0; i < d; ++i)
        {
          gsl_matrix_set(work, i, i, sqrt(gsl_ran_chisq(r, df - i)));

          for (j = 0; j < i; ++j)
            {
              gsl_matrix_set(work, i, j, gsl_ran_ugaussian(r));
            }
        }

      /* compute L * A */
      gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0,
                     L, work);

      /* compute (L * A) * (L * A)^T */
      gsl_blas_dsyrk(CblasUpper, CblasNoTrans, 1.0, work, 0.0, result);
      for (i = 0; i < d; ++i)
        {
          for (j = 0; j < i; ++j)
            {
              gsl_matrix_set(result, i, j, gsl_matrix_get(result, j, i));
            }
        }

      return GSL_SUCCESS;
    }
}

/* Compute the log of the probability density function at a given quantile
 * matrix for a Wishart distribution using the Cholesky decomposition of the
 * scale matrix.
 *
 * X       quantile matrix at which to evaluate the PDF (dimension d x d)
 * L_X     matrix resulting from the Cholesky decomposition of
 *         of the quantile matrix at which to evaluate the PDF
 *         X = L_X L_X^T (dimension d x d)
 * df      degrees of freedom
 * L       matrix resulting from the Cholesky decomposition of
 *         the scale matrix V = L L^T (dimension d x d)
 * result  output of the density (dimension 1)
 * work    matrix used for intermediate computations (dimension d x d)
 */
int
gsl_ran_wishart_log_pdf (const gsl_matrix * X,
                         const gsl_matrix * L_X,
                         const double df,
                         const gsl_matrix * L,
                         double * result,
                         gsl_matrix * work)
{
  if (L->size1 != L->size2)
    {
      GSL_ERROR("L should be a square matrix", GSL_ENOTSQR);
    }
  else if (X->size1 != X->size2)
    {
      GSL_ERROR("X should be a square matrix", GSL_ENOTSQR);
    }
  else if (L_X->size1 != L_X->size2)
    {
      GSL_ERROR("L_X should be a square matrix", GSL_ENOTSQR);
    }
  else if (X->size1 != L->size1)
    {
      GSL_ERROR("incompatible dimensions of X matrix", GSL_EBADLEN);
    }
  else if (L_X->size1 != L->size1)
    {
      GSL_ERROR("incompatible dimensions of L_X matrix", GSL_EBADLEN);
    }
  else if (df <= L->size1 - 1)
    {
      GSL_ERROR("incompatible degrees of freedom", GSL_EDOM);
    }
  else
    {
      size_t d = L->size1, i;
      int status;
      double log_mv_Ga, log_det_V, log_det_X, tr_Vinv_X;

      /* compute the log of the multivariate Gamma */
      log_mv_Ga = d * (d-1) * 0.25 * log(M_PI);
      for (i = 0; i < d; ++i)
        {
          log_mv_Ga += gsl_sf_lngamma((df - i + 1) * 0.5);
        }

      /* compute the log of the determinant of the scale matrix */
      log_det_V = log(gsl_matrix_get(L, 0, 0));
      for (i = 1; i < d; ++i)
        {
          log_det_V += log(gsl_matrix_get(L, i, i));
        }
      log_det_V = 2 * log_det_V;

      /* compute the log of the determinant of the quantile matrix */
      log_det_X = log(gsl_matrix_get(L_X, 0, 0));
      for (i = 1; i < d; ++i)
        {
          log_det_X += log(gsl_matrix_get(L_X, i, i));
        }
      log_det_X = 2 * log_det_X;

      /* compute the trace of V^(-1) X */
      status = gsl_linalg_cholesky_solve_mat(L, X, work);
      if (status)
        return status;
      tr_Vinv_X = gsl_matrix_get(work, 0, 0);
      for (i = 1; i < d; ++i)
        {
          tr_Vinv_X += gsl_matrix_get(work, i, i);
        }

      /* add all to get the log of the PDF */
      *result = - (0.5 * df * d) * log(2.0)
                - (0.5 * df) * log_det_V
                - log_mv_Ga
                + 0.5 * (df - d - 1) * log_det_X
                - 0.5 * tr_Vinv_X;

      return GSL_SUCCESS;
    }
}

int
gsl_ran_wishart_pdf (const gsl_matrix * X,
                     const gsl_matrix * L_X,
                     const double df,
                     const gsl_matrix * L,
                     double * result,
                     gsl_matrix * work)
{
  double logpdf;
  int status = gsl_ran_wishart_log_pdf(X, L_X, df, L, &logpdf, work);

  if (status == GSL_SUCCESS)
    *result = exp(logpdf);

  return status;
}
