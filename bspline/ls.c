/* bspline/ls.c
 *
 * Copyright (C) 2018, 2019, 2020 Patrick Alken
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

/* routines related to least squares approximations and B-splines */

/*
gsl_bspline_lssolve()
  Fit B-spline model to (x,y) data in a least squares sense,

min_c || y - X c ||^2

The normal equations matrix X^T X has band structure, with only
the diagonal, k-1 sub-diagonals, and k-1 super-diagonals nonzero,
where k is the spline order. This function builds the normal equations
matrix, computing only the nonzero elements and storing the diagonals
in the w->XTX matrix. The format of w->XTX is:

XTX = [ d11   d21   d31   ...   dk1  ]
      [ d22   d32   d42   ... dk-1,2 ]
      [ d33   d43   d53   ... dk-2,3 ]
      [ ...   ...   ...   ...   ...  ]
      [     dn,n-1   *     *     *   ]
      [ dnn    *     *     *     *   ]

So the first column is the diagonal of X^T X. The second column is the first
sub-diagonal, and so on. The * entries are not referenced. This is the transpose
of the format used by LAPACK for band matrices.

Inputs: x     - points at which to evaluate B-splines, length N
        y     - right hand side vector, length N
        c     - (output) solution of least squares system, length w->ncontrol
        chisq - (output) cost function chi^2
        w     - workspace

Return: success/error
*/

int
gsl_bspline_lssolve(const gsl_vector * x, const gsl_vector * y,
                    gsl_vector * c, double * chisq, gsl_bspline_workspace * w)
{
  int status = gsl_bspline_wlssolve(x, y, NULL, c, chisq, w);
  return status;
}

/*
gsl_bspline_wlssolve()
  Fit B-spline model to (x,y) data in a weighted least squares sense,

min_c || y - X c ||_W^2

The normal equations matrix X^T W X has band structure, with only
the diagonal, k-1 sub-diagonals, and k-1 super-diagonals nonzero,
where k is the spline order. This function builds the normal equations
matrix, computing only the nonzero elements and storing the diagonals
in the w->XTX matrix. The format of w->XTX is:

XTX = [ d11   d21   d31   ...   dk1  ]
      [ d22   d32   d42   ... dk-1,2 ]
      [ d33   d43   d53   ... dk-2,3 ]
      [ ...   ...   ...   ...   ...  ]
      [     dn,n-1   *     *     *   ]
      [ dnn    *     *     *     *   ]

So the first column is the diagonal of X^T W X. The second column is the first
sub-diagonal, and so on. The * entries are not referenced. This is the transpose
of the format used by LAPACK for band matrices.

Inputs: x     - points at which to evaluate B-splines, length N
        y     - right hand side vector, length N
        wts   - weight vector, length N (can be NULL)
        c     - (output) solution of least squares system, length w->ncontrol
        chisq - (output) cost function chi^2
        w     - workspace

Return: success/error

Notes:
1) On successful output, w->XTX contains the Cholesky factor for X^T W X
*/

int
gsl_bspline_wlssolve(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts,
                     gsl_vector * c, double * chisq, gsl_bspline_workspace * w)
{
  const size_t N = x->size;

  if (y->size != N)
    {
      GSL_ERROR("x and y vectors have different lengths", GSL_EBADLEN);
    }
  else if (wts != NULL && wts->size != N)
    {
      GSL_ERROR("x and weight vectors have different lengths", GSL_EBADLEN);
    }
  else if (c->size != w->ncontrol)
    {
      GSL_ERROR("coefficient vector does not match workspace", GSL_EBADLEN);
    }
  else if (N < w->ncontrol)
    {
      GSL_ERROR("data vector has too few elements", GSL_EBADLEN);
    }
  else
    {
      int status;
      gsl_vector_view XTy = gsl_vector_subvector(w->work, 0, w->ncontrol);
      double ynorm2 = 0.0; /* y^T W y */
      double tmp;
      size_t i;

      /* compute X^T W X and X^T W y */
      status = gsl_bspline_lsnormal(x, y, wts, &XTy.vector, w->XTX, w);
      if (status)
        return status;

      /* compute ||y||^2 = y^T W y */
      for (i = 0; i < N; ++i)
        {
          double yi = gsl_vector_get(y, i);
          double wi = (wts != NULL) ? gsl_vector_get(wts, i) : 1.0;
          ynorm2 += (wi * yi) * yi;
        }

      /* perform banded Cholesky decomposition of X^T X matrix */
      status = gsl_linalg_cholesky_band_decomp(w->XTX);
      if (status)
        return status;

      /* solve: X^T W X c = X^T W y */
      status = gsl_linalg_cholesky_band_solve(w->XTX, &XTy.vector, c);
      if (status)
        return status;

      /* compute: -2 c^T X^T W y */
      gsl_blas_ddot(c, &XTy.vector, &tmp);
      *chisq = -2.0 * tmp;

      /* compute L^T c */
      gsl_vector_memcpy(&XTy.vector, c);
      cblas_dtbmv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit,
                  (int) w->XTX->size1, (int) (w->XTX->size2 - 1),
                  w->XTX->data, (int) w->XTX->tda, XTy.vector.data, XTy.vector.stride);

      tmp = gsl_blas_dnrm2(&XTy.vector);

      /* chi^2 = ||y||_W^2 - 2 c^T X^T W y + || L^T c ||^2 */
      *chisq += tmp * tmp + ynorm2;

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_lsnormal()
  Compute normal equations matrix X^T W X and right hand side vector X^T W y
The normal equations matrix is symmetric and banded, with lower bandwidth k - 1.
It is stored in symmetric banded format on output.

Inputs: x   - points at which to evaluate B-splines, length N
        y   - right hand side vector, length N
        wts - weight vector, length N (can be NULL)
        XTy - (output) X^T W y vector, length ncontrol
        XTX - (output) X^T W X matrix, banded symmetric format,
              ncontrol-by-spline_order
        w   - workspace
*/

int
gsl_bspline_lsnormal(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts,
                     gsl_vector * XTy, gsl_matrix * XTX, gsl_bspline_workspace * w)
{
  const size_t N = x->size;

  if (y->size != N)
    {
      GSL_ERROR("x and y vectors have different lengths", GSL_EBADLEN);
    }
  else if (wts != NULL && wts->size != N)
    {
      GSL_ERROR("x and weight vectors have different lengths", GSL_EBADLEN);
    }
  else if (N < w->ncontrol)
    {
      GSL_ERROR("data vector has too few elements", GSL_EBADLEN);
    }
  else if (XTX->size1 != w->ncontrol || XTX->size2 != w->spline_order)
    {
      GSL_ERROR("XTX matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (XTy->size != w->ncontrol)
    {
      GSL_ERROR("XTy vector has wrong length", GSL_EBADLEN);
    }
  else
    {
      int status;
      gsl_matrix_const_view Yv = gsl_matrix_const_view_vector(y, N, 1);
      gsl_matrix_view XTY = gsl_matrix_view_vector(XTy, w->ncontrol, 1);

      status = gsl_bspline_lsnormalm(x, &Yv.matrix, wts, &XTY.matrix, XTX, w);

      return status;
    }
}

/*
gsl_bspline_lsnormalm()
  Compute normal equations matrix X^T W X and right hand side matrix X^T W Y
The normal equations matrix is symmetric and banded, with lower bandwidth k - 1.
It is stored in symmetric banded format on output.

Inputs: x   - points at which to evaluate B-splines, length N
        Y   - matrix of right hand side vectors, N-by-nrhs
        wts - weight vector, length N (can be NULL)
        XTY - (output) X^T W Y vector, ncontrol-by-nrhs
        XTX - (output) X^T W X matrix, banded symmetric format,
              ncontrol-by-spline_order
        w   - workspace
*/

int
gsl_bspline_lsnormalm(const gsl_vector * x, const gsl_matrix * Y, const gsl_vector * wts,
                      gsl_matrix * XTY, gsl_matrix * XTX, gsl_bspline_workspace * w)
{
  const size_t N = x->size;
  const size_t nrhs = Y->size2;

  if (Y->size1 != N)
    {
      GSL_ERROR("x must match Y size1", GSL_EBADLEN);
    }
  else if (wts != NULL && wts->size != N)
    {
      GSL_ERROR("x and weight vectors have different lengths", GSL_EBADLEN);
    }
  else if (N < w->ncontrol)
    {
      GSL_ERROR("data vector has too few elements", GSL_EBADLEN);
    }
  else if (XTX->size1 != w->ncontrol || XTX->size2 != w->spline_order)
    {
      GSL_ERROR("XTX matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (XTY->size1 != w->ncontrol)
    {
      GSL_ERROR("XTY matrix has wrong size1", GSL_EBADLEN);
    }
  else if (XTY->size2 != nrhs)
    {
      GSL_ERROR("XTY matrix has wrong size2", GSL_EBADLEN);
    }
  else
    {
      int status;
      size_t i;

      gsl_matrix_set_zero(XTX);
      gsl_matrix_set_zero(XTY);

      for (i = 0; i < N; ++i)
        {
          double xi = gsl_vector_get(x, i);
          gsl_vector_const_view yi = gsl_matrix_const_row(Y, i);
          double wi = (wts != NULL) ? gsl_vector_get(wts, i) : 1.0;
          size_t istart, j, kk;

          /* compute non-zero B-spline functions for this x_i */
          status = gsl_bspline_basis(xi, w->B, &istart, w);
          if (status)
            return status;

          /* add w->B into the normal equations system */
          for (j = 0; j < w->spline_order; ++j)
            {
              double Bj = gsl_vector_get(w->B, j) * wi;
              gsl_vector_view v = gsl_matrix_row(XTY, istart + j);

              /* XTY(:,istart+j) += Bj * Y(:,i) */
              gsl_blas_daxpy(Bj, &yi.vector, &v.vector);

              for (kk = 0; kk <= j; ++kk)
                {
                  double Bkk = gsl_vector_get(w->B, kk);
                  double *dptr = gsl_matrix_ptr(XTX, istart + kk, j - kk);

                  *dptr += Bj * Bkk;
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_plssolve()
  Fit a periodic B-spline to (x,y) data in a least squares sense,

min_c || y - X c ||^2

Notes:
1) The periodic condition on the coefficients is:

c_j = c_{s + 1 + j}, j = 1, ..., k - 1

where s = ncontrol - k is the number of interior knots. This
means there are

nindep = s + 1 = ncontrol - k + 1

independent coefficients which must be determined.
*/

int
gsl_bspline_plssolve(const gsl_vector * x, const gsl_vector * y,
                     gsl_vector * c, double * chisq, gsl_bspline_workspace * w)
{
  int status = gsl_bspline_pwlssolve(x, y, NULL, c, chisq, w);
  return status;
}

/*
gsl_bspline_pwlssolve()
  Fit a periodic B-spline to weighted (x,y) data in a least squares sense,

min_c || y - X c ||_W^2

Inputs: x     - vector of input x values, length N
        y     - vector of right hand side y values, length N
        wts   - weight vector, length N
        c     - (output) spline coefficients (control points), length ncontrol
        chisq - || y - X c ||_W^2
        w     - workspace

Notes:
1) The periodic condition on the coefficients is:

c_j = c_{s + 1 + j}, j = 0, ..., k - 2

where s = ncontrol - k is the number of interior knots. This
means there are

ncoef = s + 1 = ncontrol - k + 1

independent coefficients which must be determined.

2) The elements of the N-by-ncoef X are:

X_{ij} = M_j(x_i)

where:

M_j(x) = { B_{i,k}(x) + B_{i+s+1,k}(x),  i = 0, ..., k - 2
         { B_{i,k}(x),                   i = k - 1, ..., s
*/

int
gsl_bspline_pwlssolve(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts,
                      gsl_vector * c, double * chisq, gsl_bspline_workspace * w)
{
  const size_t N = x->size;

  if (y->size != N)
    {
      GSL_ERROR("x and y vectors have different lengths", GSL_EBADLEN);
    }
  else if (wts != NULL && wts->size != N)
    {
      GSL_ERROR("x and weight vectors have different lengths", GSL_EBADLEN);
    }
  else if (c->size != w->ncontrol)
    {
      GSL_ERROR("coefficient vector does not match workspace", GSL_EBADLEN);
    }
  else if (N < w->ncontrol)
    {
      GSL_ERROR("data vector has too few elements", GSL_EBADLEN);
    }
  else
    {
      const size_t ncontrol = w->ncontrol;
      const size_t spline_order = w->spline_order;
      const size_t ncoef = ncontrol - spline_order + 1;        /* number of independent coefficients */
      gsl_matrix_view R = gsl_matrix_submatrix(w->R, 0, 0, ncoef, ncoef);
      gsl_vector_view QTy = gsl_vector_subvector(c, 0, ncoef);
      double rnorm;
      size_t i;

      /* compute R and Q^T y */
      gsl_bspline_plsqr(x, y, wts, &R.matrix, &QTy.vector, &rnorm, w);

      /* solve R c = QTy */
      gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, &R.matrix, &QTy.vector);

      /* compute || y - X c ||_W^2 */
      *chisq = rnorm * rnorm;

      /* set dependent coefficients according to: c_i = c_{i + s + 1} */
      for (i = ncoef; i < ncontrol; ++i)
        {
          double ci = gsl_vector_get(c, i - ncoef);
          gsl_vector_set(c, i, ci);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_plsqr()
  Construct R and Q^T y factors for periodic B-spline least squares system

Inputs: x     - vector of input x values, length N
        y     - vector of right hand side y values, length N
        wts   - weight vector, length N
        R a   - (output) upper triangular R factor, ncoef-by-ncoef, where
                ncoef = ncontrol - spline_order + 1
        QTy   - (output) Q^T y vector, length ncoef
        rnorm - (output) residual norm || y - X c ||
        w     - workspace
*/

int
gsl_bspline_plsqr(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts,
                  gsl_matrix * R, gsl_vector * QTy, double * rnorm,
                  gsl_bspline_workspace * w)
{
  const size_t N = x->size;
  const size_t ncontrol = w->ncontrol;
  const size_t spline_order = w->spline_order;
  const size_t ncoef = ncontrol - spline_order + 1; /* number of independent coefficients */

  if (y->size != N)
    {
      GSL_ERROR("x and y vectors have different lengths", GSL_EBADLEN);
    }
  else if (wts != NULL && wts->size != N)
    {
      GSL_ERROR("x and weight vectors have different lengths", GSL_EBADLEN);
    }
  else if (R->size1 != R->size2)
    {
      GSL_ERROR("R matrix must be square", GSL_ENOTSQR);
    }
  else if (R->size1 != ncoef)
    {
      GSL_ERROR("R matrix does not match workspace", GSL_EBADLEN);
    }
  else if (QTy->size != ncoef)
    {
      GSL_ERROR("QTy vector has wrong length", GSL_EBADLEN);
    }
  else if (N < w->ncontrol)
    {
      GSL_ERROR("data vector has too few elements", GSL_EBADLEN);
    }
  else if (w->nbreak < w->spline_order)
    {
      GSL_ERROR("number of breakpoints must be >= spline_order", GSL_EDOM);
    }
  else
    {
      gsl_vector_view work = gsl_vector_subvector(w->work, 0, ncoef);
      size_t i;

      gsl_matrix_set_zero(R);
      gsl_vector_set_zero(QTy);
      *rnorm = 0.0;

      for (i = 0; i < N; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double yi = gsl_vector_get(y, i);
          size_t jstart, j, j0, j1;
          size_t ncopy, nwrap;
          gsl_vector_view v1, v2;

          gsl_bspline_basis(xi, w->B, &jstart, w);

          if (wts != NULL)
           {
             double wi = gsl_vector_get(wts, i);
             double sqrt_wi = sqrt(wi);

             gsl_vector_scale(w->B, sqrt_wi);
             yi *= sqrt_wi;
           }

          /*
           * determine how many elements of B can be directly copied into row A(i,:) and
           * how many need to be added (wrapped) to the front of row A(i,:)
           */
          if (jstart >= ncoef)
            ncopy = 0;
          else if (jstart + spline_order < ncoef)
            ncopy = spline_order;
          else
            ncopy = ncoef - jstart;

          nwrap = spline_order - ncopy;

          if (nwrap == 0)
            {
              /*
               * No elements need to be wrapped to the front of the new row,
               * we can apply a reduced Givens rotation just on the B-spline elements
               * computed
               */
              v1 = gsl_vector_subvector(&work.vector, jstart, spline_order);
              gsl_vector_memcpy(&v1.vector, w->B);
              j0 = jstart;
              j1 = jstart + spline_order - 1;
            }
          else
            {
              /*
               * We have elements wrapped to the front of the row; we must perform
               * a full Givens rotation on the row of size ncoef
               */
              gsl_vector_set_zero(&work.vector);

              if (ncopy > 0)
                {
                  v1 = gsl_vector_subvector(w->B, 0, ncopy);
                  v2 = gsl_vector_subvector(&work.vector, jstart, ncopy);
                  gsl_vector_memcpy(&v2.vector, &v1.vector);
                }

              /* add the wrapped elements */
              v1 = gsl_vector_subvector(w->B, ncopy, nwrap);
              v2 = gsl_vector_subvector(&work.vector, 0, nwrap);
              gsl_vector_add(&v2.vector, &v1.vector);

              j0 = 0;
              j1 = ncoef - 1;
            }

          /*
           * now apply Givens transformations to the new row to update the R factor
           * and QTy vector
           */
          for (j = j0; j <= j1; ++j)
            {
              double *Rjj = gsl_matrix_ptr(R, j, j);
              double *Aj = gsl_vector_ptr(&work.vector, j);
              double *QTyj_ptr = gsl_vector_ptr(QTy, j);
              double QTyj = *QTyj_ptr;
              double c, s;

              /* compute Givens transformation */
              gsl_blas_drotg(Rjj, Aj, &c, &s);

              /* apply Givens transformation to R */
              if (j < j1)
                {
                  v1 = gsl_matrix_subrow(R, j, j + 1, j1 - j);
                  v2 = gsl_vector_subvector(&work.vector, j + 1, j1 - j);
                  gsl_blas_drot(&v1.vector, &v2.vector, c, s);
                }

              /* apply Givens transformation to rhs vector */
              *QTyj_ptr = c * QTyj + s * yi;
              yi = -s * QTyj + c * yi;
            }

          *rnorm = gsl_hypot(*rnorm, yi);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_residuals()
  Calculate residual vector,

r = y - X c

Inputs: x - points at which to evaluate B-splines, length N
        y - right hand side vector, length N
        c - solution of least squares system, length w->ncontrol
        r - (output) residual vector, length N
        w - workspace
*/

int
gsl_bspline_residuals(const gsl_vector * x, const gsl_vector * y, const gsl_vector * c, gsl_vector * r,
                      gsl_bspline_workspace * w)
{
  const size_t N = x->size;

  if (N != y->size)
    {
      GSL_ERROR("x and y vectors must be same size", GSL_EBADLEN);
    }
  else if (c->size != w->ncontrol)
    {
      GSL_ERROR("coefficient vector does not match workspace", GSL_EBADLEN);
    }
  else if (N != r->size)
    {
      GSL_ERROR("right hand side vector does not match residual vector", GSL_EBADLEN);
    }
  else
    {
      int status;
      size_t i;

      for (i = 0; i < N; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double yi = gsl_vector_get(y, i);
          double result;

          status = gsl_bspline_calc(xi, c, &result, w);
          if (status)
            return status;

          gsl_vector_set(r, i, yi - result);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_covariance()
  Calculate covariance matrix of parameters of least squares fit

Inputs: XTX - banded symmetric matrix with Cholesky factor stored
        cov - (output) covariance matrix, ncontrol-by-ncontrol
        w   - workspace

Notes:
1) gsl_bspline_lssolve() or gsl_bspline_wlssolve() must be called
first to compute normal equations matrix and Cholesky decomposition
*/

int
gsl_bspline_covariance(const gsl_matrix * XTX, gsl_matrix * cov, gsl_bspline_workspace * w)
{
  if (XTX->size1 != w->ncontrol || XTX->size2 != w->spline_order)
    {
      GSL_ERROR("XTX matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {
      GSL_ERROR("covariance matrix must be square", GSL_ENOTSQR);
    }
  else if (cov->size1 != w->ncontrol)
    {
      GSL_ERROR("covariance matrix does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status = gsl_linalg_cholesky_band_invert(XTX, cov);
      return status;
    }
}

/*
gsl_bspline_err()
  Compute 1-sigma spline error due to variances on model coefficients

f(x) = \sum_i c_i B_i(x)

From error propagation:

\delta f(x) = \sqrt B^T Cov_c B

and

\delta d^n/dx^n f(x) = \sqrt [d^n/dx^n B]^T Cov_c [d^n/dx^n B]

Inputs: x      - point at which to evaluate spline
        nderiv - derivative order
        cov    - parameter covariance matrix, ncontrol-by-ncontrol
        err    - (output) standard deviation of spline at x due to control point variances
        w      - workspace
*/

int
gsl_bspline_err(const double x, const size_t nderiv,
                const gsl_matrix * cov, double * err,
                gsl_bspline_workspace * w)
{
  if (cov->size1 != cov->size2)
    {
      GSL_ERROR("covariance matrix must be square", GSL_ENOTSQR);
    }
  else if (cov->size1 != w->ncontrol)
    {
      GSL_ERROR("covariance matrix does not match workspace", GSL_EBADLEN);
    }
  else if (nderiv >= w->spline_order)
    {
      /* quick return */
      *err = 0.0;
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t k = w->spline_order;
      gsl_vector_const_view B = gsl_matrix_const_column(w->dB, nderiv);
      size_t istart;

      status = gsl_bspline_basis_deriv(x, nderiv, w->dB, &istart, w);
      if (status)
        return status;

      {
        gsl_matrix_const_view C = gsl_matrix_const_submatrix(cov, istart, istart, k, k);
        gsl_vector * tmp = w->B;

        /* compute B^T C B, using w->B as temporary workspace */
        gsl_blas_dsymv(CblasLower, 1.0, &C.matrix, &B.vector, 0.0, tmp);
        gsl_blas_ddot(&B.vector, tmp, err);

        *err = sqrt(*err);
      }

      return GSL_SUCCESS;
    }
}

int
gsl_bspline_rcond(const gsl_matrix * XTX, double * rcond, gsl_bspline_workspace * w)
{
  if (XTX->size1 != w->ncontrol || XTX->size2 != w->spline_order)
    {
      GSL_ERROR("XTX matrix has wrong dimensions", GSL_EBADLEN);
    }
  else
    {
      int status = gsl_linalg_cholesky_band_rcond(XTX, rcond, w->work);

      if (status)
         return status;

      return GSL_SUCCESS;
    }
}
