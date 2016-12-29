#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

/* Fit
 *
 * y = X c
 *
 * where X is an n x p matrix of n observations for p variables.
 *
 * The solution includes a possible standard form Tikhonov regularization:
 *
 * c = (X^T X + lambda^2 I)^{-1} X^T y
 *
 * where lambda^2 is the Tikhonov regularization parameter.
 *
 * The function multifit_linear_svd() must first be called to
 * compute the SVD decomposition of X
 *
 * Inputs: X        - least squares matrix
 *         y        - right hand side vector
 *         tol      - singular value tolerance
 *         lambda   - Tikhonov regularization parameter lambda;
 *                    ignored if <= 0
 *         rank     - (output) effective rank
 *         c        - (output) model coefficient vector
 *         rnorm    - (output) residual norm ||y - X c||
 *         snorm    - (output) solution norm ||c||
 *         work     - workspace
 *
 * Notes:
 * 1) The dimensions of X must match work->n and work->p which are set
 *    by multifit_linear_svd()
 * 2) On input:
 *    work->A contains U
 *    work->Q contains Q
 *    work->S contains singular values
 * 3) If this function is called from gsl_multifit_wlinear(), then
 *    the input y points to work->t, which contains sqrt(W) y. Since
 *    work->t is also used as scratch workspace by this function, we
 *    do the necessary computations with y first to avoid problems.
 * 4) When lambda <= 0, singular values are truncated when:
 *    s_j <= tol * s_0
 */

static int
multifit_linear_solve (const gsl_matrix * X,
                       const gsl_vector * y,
                       const double tol,
                       const double lambda,
                       size_t * rank,
                       gsl_vector * c,
                       double *rnorm,
                       double *snorm,
                       gsl_multifit_linear_workspace * work)
{
  const size_t n = X->size1;
  const size_t p = X->size2;

  if (n != work->n || p != work->p)
    {
      GSL_ERROR("observation matrix does not match workspace", GSL_EBADLEN);
    }
  else if (n != y->size)
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
      const double lambda_sq = lambda * lambda;

      double rho_ls = 0.0;     /* contribution to rnorm from OLS */

      size_t j, p_eff;

      /* these inputs are previously computed by multifit_linear_svd() */
      gsl_matrix_view A = gsl_matrix_submatrix(work->A, 0, 0, n, p);
      gsl_matrix_view Q = gsl_matrix_submatrix(work->Q, 0, 0, p, p);
      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);

      /* workspace */
      gsl_matrix_view QSI = gsl_matrix_submatrix(work->QSI, 0, 0, p, p);
      gsl_vector_view xt = gsl_vector_subvector(work->xt, 0, p);
      gsl_vector_view D = gsl_vector_subvector(work->D, 0, p);
      gsl_vector_view t = gsl_vector_subvector(work->t, 0, n);

      /*
       * Solve y = A c for c
       * c = Q diag(s_i / (s_i^2 + lambda_i^2)) U^T y
       */

      /* compute xt = U^T y */
      gsl_blas_dgemv (CblasTrans, 1.0, &A.matrix, y, 0.0, &xt.vector);

      if (n > p)
        {
          /*
           * compute OLS residual norm = || y - U U^T y ||;
           * for n = p, U U^T = I, so no need to calculate norm
           */
          gsl_vector_memcpy(&t.vector, y);
          gsl_blas_dgemv(CblasNoTrans, -1.0, &A.matrix, &xt.vector, 1.0, &t.vector);
          rho_ls = gsl_blas_dnrm2(&t.vector);
        }

      if (lambda > 0.0)
        {
          /* xt <-- [ s(i) / (s(i)^2 + lambda^2) ] .* U^T y */
          for (j = 0; j < p; ++j)
            {
              double sj = gsl_vector_get(&S.vector, j);
              double f = (sj * sj) / (sj * sj + lambda_sq);
              double *ptr = gsl_vector_ptr(&xt.vector, j);

              /* use D as workspace for residual norm */
              gsl_vector_set(&D.vector, j, (1.0 - f) * (*ptr));

              *ptr *= sj / (sj*sj + lambda_sq);
            }

          /* compute regularized solution vector */
          gsl_blas_dgemv (CblasNoTrans, 1.0, &Q.matrix, &xt.vector, 0.0, c);

          /* compute solution norm */
          *snorm = gsl_blas_dnrm2(c);

          /* compute residual norm */
          *rnorm = gsl_blas_dnrm2(&D.vector);

          if (n > p)
            {
              /* add correction to residual norm (see eqs 6-7 of [1]) */
              *rnorm = sqrt((*rnorm) * (*rnorm) + rho_ls * rho_ls);
            }

          /* reset D vector */
          gsl_vector_set_all(&D.vector, 1.0);
        }
      else
        {
          /* Scale the matrix Q, QSI = Q S^{-1} */

          gsl_matrix_memcpy (&QSI.matrix, &Q.matrix);

          {
            double s0 = gsl_vector_get (&S.vector, 0);
            p_eff = 0;

            for (j = 0; j < p; j++)
              {
                gsl_vector_view column = gsl_matrix_column (&QSI.matrix, j);
                double sj = gsl_vector_get (&S.vector, j);
                double alpha;

                if (sj <= tol * s0)
                  {
                    alpha = 0.0;
                  }
                else
                  {
                    alpha = 1.0 / sj;
                    p_eff++;
                  }

                gsl_vector_scale (&column.vector, alpha);
              }

            *rank = p_eff;
          }

          gsl_blas_dgemv (CblasNoTrans, 1.0, &QSI.matrix, &xt.vector, 0.0, c);

          /* Unscale the balancing factors */
          gsl_vector_div (c, &D.vector);

          *snorm = gsl_blas_dnrm2(c);
          *rnorm = rho_ls;
        }

      return GSL_SUCCESS;
    }
}
