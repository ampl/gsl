/* multifit/multireg.c
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

/*
 * References:
 *
 * [1] P. C. Hansen & D. P. O'Leary, "The use of the L-curve in
 * the regularization of discrete ill-posed problems",  SIAM J. Sci.
 * Comput. 14 (1993), pp. 1487-1503.
 *
 * [2] P. C. Hansen, "Discrete Inverse Problems: Insight and Algorithms,"
 * SIAM Press, 2010.
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "linear_common.c"

static int menger(const double a, const double b, const double c, double * kappa);

int
gsl_multifit_linear_solve (const double lambda,
                           const gsl_matrix * X,
                           const gsl_vector * y,
                           gsl_vector * c,
                           double *rnorm,
                           double *snorm,
                           gsl_multifit_linear_workspace * work)
{
  size_t rank;
  int status;

  status = multifit_linear_solve(X, y, GSL_DBL_EPSILON, lambda, &rank, c,
                                 rnorm, snorm, work);

  return status;
} /* gsl_multifit_linear_solve() */

/*
gsl_multifit_linear_applyW()
  Apply weight matrix to (X,y) LS system

Inputs: X    - least squares matrix n-by-p
        w    - weight vector n-by-1 or NULL for W = I
        y    - right hand side n-by-1
        WX   - (output) sqrt(W) X, n-by-p
        Wy   - (output) sqrt(W) y, n-by-1

Notes:
1) If w = NULL, on output WX = X and Wy = y
2) It is allowed for WX = X and Wy = y for in-place transform
*/

int
gsl_multifit_linear_applyW(const gsl_matrix * X,
                           const gsl_vector * w,
                           const gsl_vector * y,
                           gsl_matrix * WX,
                           gsl_vector * Wy)
{
  const size_t n = X->size1;
  const size_t p = X->size2;

  if (n != y->size)
    {
      GSL_ERROR("y vector does not match X", GSL_EBADLEN);
    }
  else if (w != NULL && n != w->size)
    {
      GSL_ERROR("weight vector does not match X", GSL_EBADLEN);
    }
  else if (n != WX->size1 || p != WX->size2)
    {
      GSL_ERROR("WX matrix dimensions do not match X", GSL_EBADLEN);
    }
  else if (n != Wy->size)
    {
      GSL_ERROR("Wy vector must be length n", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* copy WX = X; Wy = y if distinct pointers */
      if (WX != X)
        gsl_matrix_memcpy(WX, X);
      if (Wy != y)
        gsl_vector_memcpy(Wy, y);

      if (w != NULL)
        {
          /* construct WX = sqrt(W) X and Wy = sqrt(W) y */
          for (i = 0; i < n; ++i)
            {
              double wi = gsl_vector_get(w, i);
              double swi;
              gsl_vector_view row = gsl_matrix_row(WX, i);
              double *yi = gsl_vector_ptr(Wy, i);

              if (wi < 0.0)
                wi = 0.0;

              swi = sqrt(wi);
              gsl_vector_scale(&row.vector, swi);
              *yi *= swi;
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_wstdform1()
  Using regularization matrix
L = diag(l_1,l_2,...,l_p), transform to Tikhonov standard form:

X~ = sqrt(W) X L^{-1}
y~ = sqrt(W) y
c~ = L c

Inputs: L    - Tikhonov matrix as a vector of diagonal elements p-by-1;
               or NULL for L = I
        X    - least squares matrix n-by-p
        y    - right hand side vector n-by-1
        w    - weight vector n-by-1; or NULL for W = I
        Xs   - least squares matrix in standard form X~ n-by-p
        ys   - right hand side vector in standard form y~ n-by-1
        work - workspace

Return: success/error

Notes:
1) It is allowed for X = Xs and y = ys
*/

int
gsl_multifit_linear_wstdform1 (const gsl_vector * L,
                               const gsl_matrix * X,
                               const gsl_vector * w,
                               const gsl_vector * y,
                               gsl_matrix * Xs,
                               gsl_vector * ys,
                               gsl_multifit_linear_workspace * work)
{
  const size_t n = X->size1;
  const size_t p = X->size2;

  if (n > work->nmax || p > work->pmax)
    {
      GSL_ERROR("observation matrix larger than workspace", GSL_EBADLEN);
    }
  else if (L != NULL && p != L->size)
    {
      GSL_ERROR("L vector does not match X", GSL_EBADLEN);
    }
  else if (n != y->size)
    {
      GSL_ERROR("y vector does not match X", GSL_EBADLEN);
    }
  else if (w != NULL && n != w->size)
    {
      GSL_ERROR("weight vector does not match X", GSL_EBADLEN);
    }
  else if (n != Xs->size1 || p != Xs->size2)
    {
      GSL_ERROR("Xs matrix dimensions do not match X", GSL_EBADLEN);
    }
  else if (n != ys->size)
    {
      GSL_ERROR("ys vector must be length n", GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;

      /* compute Xs = sqrt(W) X and ys = sqrt(W) y */
      status = gsl_multifit_linear_applyW(X, w, y, Xs, ys);
      if (status)
        return status;

      if (L != NULL)
        {
          size_t j;

          /* construct X~ = sqrt(W) X * L^{-1} matrix */
          for (j = 0; j < p; ++j)
            {
              gsl_vector_view Xj = gsl_matrix_column(Xs, j);
              double lj = gsl_vector_get(L, j);

              if (lj == 0.0)
                {
                  GSL_ERROR("L matrix is singular", GSL_EDOM);
                }

              gsl_vector_scale(&Xj.vector, 1.0 / lj);
            }
        }

      return status;
    }
}

/*
gsl_multifit_linear_stdform1()
  Using regularization matrix L = diag(l_1,l_2,...,l_p),
and W = I, transform to Tikhonov standard form:

X~ = X L^{-1}
y~ = y
c~ = L c

Inputs: L    - Tikhonov matrix as a vector of diagonal elements p-by-1
        X    - least squares matrix n-by-p
        y    - right hand side vector n-by-1
        Xs   - least squares matrix in standard form X~ n-by-p
        ys   - right hand side vector in standard form y~ n-by-1
        work - workspace

Return: success/error

Notes:
1) It is allowed for X = Xs
*/

int
gsl_multifit_linear_stdform1 (const gsl_vector * L,
                              const gsl_matrix * X,
                              const gsl_vector * y,
                              gsl_matrix * Xs,
                              gsl_vector * ys,
                              gsl_multifit_linear_workspace * work)
{
  int status;

  status = gsl_multifit_linear_wstdform1(L, X, NULL, y, Xs, ys, work);

  return status;
}

int
gsl_multifit_linear_L_decomp (gsl_matrix * L, gsl_vector * tau)
{
  const size_t m = L->size1;
  const size_t p = L->size2;
  int status;

  if (tau->size != GSL_MIN(m, p))
    {
      GSL_ERROR("tau vector must be min(m,p)", GSL_EBADLEN);
    }
  else if (m >= p)
    {
      /* square or tall L matrix */
      status = gsl_linalg_QR_decomp(L, tau);
      return status;
    }
  else
    {
      /* more columns than rows, compute qr(L^T) */
      gsl_matrix_view LTQR = gsl_matrix_view_array(L->data, p, m);
      gsl_matrix *LT = gsl_matrix_alloc(p, m);

      /* XXX: use temporary storage due to difficulties in transforming
       * a rectangular matrix in-place */
      gsl_matrix_transpose_memcpy(LT, L);
      gsl_matrix_memcpy(&LTQR.matrix, LT);
      gsl_matrix_free(LT);

      status = gsl_linalg_QR_decomp(&LTQR.matrix, tau);

      return status;
    }
}

/*
gsl_multifit_linear_wstdform2()
  Using regularization matrix L which is m-by-p, transform to Tikhonov
standard form. This routine is separated into two cases:

Case 1: m >= p, here we can use the QR decomposition of L = QR, and note
that ||L c|| = ||R c|| where R is p-by-p. Therefore,

X~ = X R^{-1} is n-by-p
y~ = y is n-by-1
c~ is p-by-1
M is not used

Case 2: m < p

X~ is (n - p + m)-by-m
y~ is (n - p + m)-by-1
c~ is m-by-1
M is n-by-p (workspace)

Inputs: LQR  - output from gsl_multifit_linear_L_decomp()
        Ltau - output from gsl_multifit_linear_L_decomp()
        X    - least squares matrix n-by-p
        w    - weight vector n-by-1; or NULL for W = I
        y    - right hand side vector n-by-1
        Xs   - (output) least squares matrix in standard form
               case 1: n-by-p
               case 2: (n - p + m)-by-m
        ys   - (output) right hand side vector in standard form
               case 1: n-by-1
               case 2: (n - p + m)-by-1
        M    - (output) workspace matrix needed to reconstruct solution vector
               case 1: not used
               case 2: n-by-p
        work - workspace

Return: success/error

Notes:
1) If m >= p, on output:
   Xs = X R^{-1}
   ys = y

2) If m < p, on output:
   M(:,1:pm) contains QR decomposition of A * K_o, needed to reconstruct
   solution vector, where pm = p - m; M(:,p) contains Householder scalars
*/

int
gsl_multifit_linear_wstdform2 (const gsl_matrix * LQR,
                               const gsl_vector * Ltau,
                               const gsl_matrix * X,
                               const gsl_vector * w,
                               const gsl_vector * y,
                               gsl_matrix * Xs,
                               gsl_vector * ys,
                               gsl_matrix * M,
                               gsl_multifit_linear_workspace * work)
{
  const size_t m = LQR->size1;
  const size_t n = X->size1;
  const size_t p = X->size2;

  if (n > work->nmax || p > work->pmax)
    {
      GSL_ERROR("observation matrix larger than workspace", GSL_EBADLEN);
    }
  else if (p != LQR->size2)
    {
      GSL_ERROR("LQR and X matrices have different numbers of columns", GSL_EBADLEN);
    }
  else if (n != y->size)
    {
      GSL_ERROR("y vector does not match X", GSL_EBADLEN);
    }
  else if (w != NULL && n != w->size)
    {
      GSL_ERROR("weights vector must be length n", GSL_EBADLEN);
    }
  else if (m >= p) /* square or tall L matrix */
    {
      /* the sizes of Xs and ys depend on whether m >= p or m < p */
      if (n != Xs->size1 || p != Xs->size2)
        {
          GSL_ERROR("Xs matrix must be n-by-p", GSL_EBADLEN);
        }
      else if (n != ys->size)
        {
          GSL_ERROR("ys vector must have length n", GSL_EBADLEN);
        }
      else
        {
          int status;
          size_t i;
          gsl_matrix_const_view R = gsl_matrix_const_submatrix(LQR, 0, 0, p, p);

          /* compute Xs = sqrt(W) X and ys = sqrt(W) y */
          status = gsl_multifit_linear_applyW(X, w, y, Xs, ys);
          if (status)
            return status;

          /* compute X~ = X R^{-1} using QR decomposition of L */
          for (i = 0; i < n; ++i)
            {
              gsl_vector_view v = gsl_matrix_row(Xs, i);

              /* solve: R^T y = X_i */
              gsl_blas_dtrsv(CblasUpper, CblasTrans, CblasNonUnit, &R.matrix, &v.vector);
            }

          return GSL_SUCCESS;
        }
    }
  else /* L matrix with m < p */
    {
      const size_t pm = p - m;
      const size_t npm = n - pm;

      /*
       * This code closely follows section 2.6.1 of Hansen's
       * "Regularization Tools" manual
       */

      if (npm != Xs->size1 || m != Xs->size2)
        {
          GSL_ERROR("Xs matrix must be (n-p+m)-by-m", GSL_EBADLEN);
        }
      else if (npm != ys->size)
        {
          GSL_ERROR("ys vector must be of length (n-p+m)", GSL_EBADLEN);
        }
      else if (n != M->size1 || p != M->size2)
        {
          GSL_ERROR("M matrix must be n-by-p", GSL_EBADLEN);
        }
      else
        {
          int status;
          gsl_matrix_view A = gsl_matrix_submatrix(work->A, 0, 0, n, p);
          gsl_vector_view b = gsl_vector_subvector(work->t, 0, n);

          gsl_matrix_view LTQR = gsl_matrix_view_array(LQR->data, p, m);           /* qr(L^T) */
          gsl_matrix_view Rp = gsl_matrix_view_array(LQR->data, m, m);             /* R factor of L^T */
          gsl_vector_const_view LTtau = gsl_vector_const_subvector(Ltau, 0, m);

          /*
           * M(:,1:p-m) will hold QR decomposition of A K_o; M(:,p) will hold
           * Householder scalars
           */
          gsl_matrix_view MQR = gsl_matrix_submatrix(M, 0, 0, n, pm);
          gsl_vector_view Mtau = gsl_matrix_subcolumn(M, p - 1, 0, GSL_MIN(n, pm));

          gsl_matrix_view AKo, AKp, HqTAKp;
          gsl_vector_view v;
          size_t i;

          /* compute A = sqrt(W) X and b = sqrt(W) y */
          status = gsl_multifit_linear_applyW(X, w, y, &A.matrix, &b.vector);
          if (status)
            return status;

          /* compute: A <- A K = [ A K_p ; A K_o ] */
          gsl_linalg_QR_matQ(&LTQR.matrix, &LTtau.vector, &A.matrix);
          AKp = gsl_matrix_submatrix(&A.matrix, 0, 0, n, m); 
          AKo = gsl_matrix_submatrix(&A.matrix, 0, m, n, pm); 

          /* compute QR decomposition [H,T] = qr(A * K_o) and store in M */
          gsl_matrix_memcpy(&MQR.matrix, &AKo.matrix);
          gsl_linalg_QR_decomp(&MQR.matrix, &Mtau.vector);

          /* AKp currently contains A K_p; apply H^T from the left to get H^T A K_p */
          gsl_linalg_QR_QTmat(&MQR.matrix, &Mtau.vector, &AKp.matrix);

          /* the last npm rows correspond to H_q^T A K_p */
          HqTAKp = gsl_matrix_submatrix(&AKp.matrix, pm, 0, npm, m);

          /* solve: Xs R_p^T = H_q^T A K_p for Xs */
          gsl_matrix_memcpy(Xs, &HqTAKp.matrix);
          for (i = 0; i < npm; ++i)
            {
              gsl_vector_view x = gsl_matrix_row(Xs, i);
              gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, &Rp.matrix, &x.vector);
            }

          /*
           * compute: ys = H_q^T b; this is equivalent to computing
           * the last q elements of H^T b (q = npm)
           */
          v = gsl_vector_subvector(&b.vector, pm, npm);
          gsl_linalg_QR_QTvec(&MQR.matrix, &Mtau.vector, &b.vector);
          gsl_vector_memcpy(ys, &v.vector);

          return GSL_SUCCESS;
        }
    }
}

int
gsl_multifit_linear_stdform2 (const gsl_matrix * LQR,
                              const gsl_vector * Ltau,
                              const gsl_matrix * X,
                              const gsl_vector * y,
                              gsl_matrix * Xs,
                              gsl_vector * ys,
                              gsl_matrix * M,
                              gsl_multifit_linear_workspace * work)
{
  int status;

  status = gsl_multifit_linear_wstdform2(LQR, Ltau, X, NULL, y, Xs, ys, M, work);

  return status;
}

/*
gsl_multifit_linear_genform1()
  Backtransform regularized solution vector using matrix
L = diag(L)
*/

int
gsl_multifit_linear_genform1 (const gsl_vector * L,
                              const gsl_vector * cs,
                              gsl_vector * c,
                              gsl_multifit_linear_workspace * work)
{
  if (L->size > work->pmax)
    {
      GSL_ERROR("L vector does not match workspace", GSL_EBADLEN);
    }
  else if (L->size != cs->size)
    {
      GSL_ERROR("cs vector does not match L", GSL_EBADLEN);
    }
  else if (L->size != c->size)
    {
      GSL_ERROR("c vector does not match L", GSL_EBADLEN);
    }
  else
    {
      /* compute true solution vector c = L^{-1} c~ */
      gsl_vector_memcpy(c, cs);
      gsl_vector_div(c, L);

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_wgenform2()
  Backtransform regularized solution vector in standard form to recover
original vector

Inputs: LQR  - output from gsl_multifit_linear_L_decomp()
        Ltau - output from gsl_multifit_linear_L_decomp()
        X    - original least squares matrix n-by-p
        w    - original weight vector n-by-1 or NULL for W = I
        y    - original rhs vector n-by-1
        cs   - standard form solution vector
        c    - (output) original solution vector p-by-1
        M    - matrix computed by gsl_multifit_linear_wstdform2()
        work - workspace
*/

int
gsl_multifit_linear_wgenform2 (const gsl_matrix * LQR,
                               const gsl_vector * Ltau,
                               const gsl_matrix * X,
                               const gsl_vector * w,
                               const gsl_vector * y,
                               const gsl_vector * cs,
                               const gsl_matrix * M,
                               gsl_vector * c,
                               gsl_multifit_linear_workspace * work)
{
  const size_t m = LQR->size1;
  const size_t n = X->size1;
  const size_t p = X->size2;

  if (n > work->nmax || p > work->pmax)
    {
      GSL_ERROR("X matrix does not match workspace", GSL_EBADLEN);
    }
  else if (p != LQR->size2)
    {
      GSL_ERROR("LQR matrix does not match X", GSL_EBADLEN);
    }
  else if (p != c->size)
    {
      GSL_ERROR("c vector does not match X", GSL_EBADLEN);
    }
  else if (w != NULL && n != w->size)
    {
      GSL_ERROR("w vector does not match X", GSL_EBADLEN);
    }
  else if (n != y->size)
    {
      GSL_ERROR("y vector does not match X", GSL_EBADLEN);
    }
  else if (m >= p)                    /* square or tall L matrix */
    {
      if (p != cs->size)
        {
          GSL_ERROR("cs vector must be length p", GSL_EBADLEN);
        }
      else
        {
          int s;
          gsl_matrix_const_view R = gsl_matrix_const_submatrix(LQR, 0, 0, p, p); /* R factor of L */

          /* solve R c = cs for true solution c, using QR decomposition of L */
          gsl_vector_memcpy(c, cs);
          s = gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, &R.matrix, c);

          return s;
        }
    }
  else                                /* rectangular L matrix with m < p */
    {
      if (m != cs->size)
        {
          GSL_ERROR("cs vector must be length m", GSL_EBADLEN);
        }
      else if (n != M->size1 || p != M->size2)
        {
          GSL_ERROR("M matrix must be size n-by-p", GSL_EBADLEN);
        }
      else
        {
          int status;
          const size_t pm = p - m;
          gsl_matrix_view A = gsl_matrix_submatrix(work->A, 0, 0, n, p);
          gsl_vector_view b = gsl_vector_subvector(work->t, 0, n);
          gsl_matrix_view Rp = gsl_matrix_view_array(LQR->data, m, m); /* R_p */
          gsl_matrix_view LTQR = gsl_matrix_view_array(LQR->data, p, m);
          gsl_vector_const_view LTtau = gsl_vector_const_subvector(Ltau, 0, m);
          gsl_matrix_const_view MQR = gsl_matrix_const_submatrix(M, 0, 0, n, pm);
          gsl_vector_const_view Mtau = gsl_matrix_const_subcolumn(M, p - 1, 0, GSL_MIN(n, pm));
          gsl_matrix_const_view To = gsl_matrix_const_submatrix(&MQR.matrix, 0, 0, pm, pm);
          gsl_vector_view workp = gsl_vector_subvector(work->xt, 0, p);
          gsl_vector_view v1, v2;

          /* compute A = sqrt(W) X and b = sqrt(W) y */
          status = gsl_multifit_linear_applyW(X, w, y, &A.matrix, &b.vector);
          if (status)
            return status;

          /* initialize c to zero */
          gsl_vector_set_zero(c);

          /* compute c = L_inv cs = K_p R_p^{-T} cs */

          /* set c(1:m) = R_p^{-T} cs */
          v1 = gsl_vector_subvector(c, 0, m);
          gsl_vector_memcpy(&v1.vector, cs);
          gsl_blas_dtrsv(CblasUpper, CblasTrans, CblasNonUnit, &Rp.matrix, &v1.vector);

          /* c <- K R_p^{-T} cs = [ K_p R_p^{_T} cs ; 0 ] */
          gsl_linalg_QR_Qvec(&LTQR.matrix, &LTtau.vector, c);

          /* compute: b1 = b - A L_inv cs */
          gsl_blas_dgemv(CblasNoTrans, -1.0, &A.matrix, c, 1.0, &b.vector);

          /* compute: b2 = H^T b1 */
          gsl_linalg_QR_QTvec(&MQR.matrix, &Mtau.vector, &b.vector);

          /* compute: b3 = T_o^{-1} b2 */
          v1 = gsl_vector_subvector(&b.vector, 0, pm);
          gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, &To.matrix, &v1.vector);

          /* compute: b4 = K_o b3 */
          gsl_vector_set_zero(&workp.vector);
          v2 = gsl_vector_subvector(&workp.vector, m, pm);
          gsl_vector_memcpy(&v2.vector, &v1.vector);
          gsl_linalg_QR_Qvec(&LTQR.matrix, &LTtau.vector, &workp.vector);

          /* final solution vector */
          gsl_vector_add(c, &workp.vector);

          return GSL_SUCCESS;
        }
    }
}

int
gsl_multifit_linear_genform2 (const gsl_matrix * LQR,
                              const gsl_vector * Ltau,
                              const gsl_matrix * X,
                              const gsl_vector * y,
                              const gsl_vector * cs,
                              const gsl_matrix * M,
                              gsl_vector * c,
                              gsl_multifit_linear_workspace * work)
{
  int status;

  status = gsl_multifit_linear_wgenform2(LQR, Ltau, X, NULL, y, cs, M, c, work);

  return status;
}

/*
gsl_multifit_linear_lreg()
  Calculate regularization parameters to use in L-curve
analysis

Inputs: smin      - smallest singular value of LS system
        smax      - largest singular value of LS system > 0
        reg_param - (output) vector of regularization parameters
                    derived from singular values

Return: success/error

Notes:
1) reg_param is ordered in ascending order with:
reg_param(1) = smin
reg_param(N) = smax
*/

int
gsl_multifit_linear_lreg (const double smin, const double smax,
                          gsl_vector * reg_param)
{
  if (smax <= 0.0)
    {
      GSL_ERROR("smax must be positive", GSL_EINVAL);
    }
  else
    {
      const size_t N = reg_param->size;

      /* smallest regularization parameter */
      const double smin_ratio = 16.0 * GSL_DBL_EPSILON;
      const double new_smin = GSL_MAX(smin, smax*smin_ratio);
      double ratio;
      size_t i;

      /* reg_param(1) = smin */
      gsl_vector_set(reg_param, 0, new_smin);

      /* ratio so that reg_param(N) = smax */
      ratio = pow(smax / new_smin, 1.0 / ((double)N - 1.0));

      for (i = 1; i < N; ++i)
        {
          double rm1 = gsl_vector_get(reg_param, i - 1);
          gsl_vector_set(reg_param, i, ratio * rm1);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_lcurve()
  Calculate L-curve using regularization parameters estimated
from singular values of least squares matrix

Inputs: y         - right hand side vector
        reg_param - (output) vector of regularization parameters
                    derived from singular values
        rho       - (output) vector of residual norms ||y - X c||
        eta       - (output) vector of solution norms ||lambda c||
        work      - workspace

Return: success/error

Notes:
1) SVD of X must be computed first by calling multifit_linear_svd();
   work->n and work->p are initialized by this function
*/

int
gsl_multifit_linear_lcurve (const gsl_vector * y,
                            gsl_vector * reg_param,
                            gsl_vector * rho, gsl_vector * eta,
                            gsl_multifit_linear_workspace * work)
{
  const size_t n = y->size;
  const size_t N = rho->size; /* number of points on L-curve */

  if (n != work->n)
    {
      GSL_ERROR("y vector does not match workspace", GSL_EBADLEN);
    }
  else if (N < 3)
    {
      GSL_ERROR ("at least 3 points are needed for L-curve analysis",
                 GSL_EBADLEN);
    }
  else if (N != eta->size)
    {
      GSL_ERROR ("size of rho and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else if (reg_param->size != eta->size)
    {
      GSL_ERROR ("size of reg_param and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;
      const size_t p = work->p;

      size_t i, j;

      gsl_matrix_view A = gsl_matrix_submatrix(work->A, 0, 0, n, p);
      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);
      gsl_vector_view xt = gsl_vector_subvector(work->xt, 0, p);
      gsl_vector_view workp = gsl_matrix_subcolumn(work->QSI, 0, 0, p);
      gsl_vector_view workp2 = gsl_vector_subvector(work->D, 0, p); /* D isn't used for regularized problems */

      const double smax = gsl_vector_get(&S.vector, 0);
      const double smin = gsl_vector_get(&S.vector, p - 1);

      double dr; /* residual error from projection */
      double normy = gsl_blas_dnrm2(y);
      double normUTy;

      /* compute projection xt = U^T y */
      gsl_blas_dgemv (CblasTrans, 1.0, &A.matrix, y, 0.0, &xt.vector);

      normUTy = gsl_blas_dnrm2(&xt.vector);
      dr = normy*normy - normUTy*normUTy;

      /* calculate regularization parameters */
      gsl_multifit_linear_lreg(smin, smax, reg_param);

      for (i = 0; i < N; ++i)
        {
          double lambda = gsl_vector_get(reg_param, i);
          double lambda_sq = lambda * lambda;

          for (j = 0; j < p; ++j)
            {
              double sj = gsl_vector_get(&S.vector, j);
              double xtj = gsl_vector_get(&xt.vector, j);
              double f = sj / (sj*sj + lambda_sq);

              gsl_vector_set(&workp.vector, j, f * xtj);
              gsl_vector_set(&workp2.vector, j, (1.0 - sj*f) * xtj);
            }

          gsl_vector_set(eta, i, gsl_blas_dnrm2(&workp.vector));
          gsl_vector_set(rho, i, gsl_blas_dnrm2(&workp2.vector));
        }

      if (n > p && dr > 0.0)
        {
          /* add correction to residual norm (see eqs 6-7 of [1]) */
          for (i = 0; i < N; ++i)
            {
              double rhoi = gsl_vector_get(rho, i);
              double *ptr = gsl_vector_ptr(rho, i);

              *ptr = sqrt(rhoi*rhoi + dr);
            }
        }

      /* restore D to identity matrix */
      gsl_vector_set_all(work->D, 1.0);

      return status;
    }
} /* gsl_multifit_linear_lcurve() */

/*
gsl_multifit_linear_lcurvature()
  Calculate curvature of previously computed L-curve as a function of
regularization parameter

Inputs: y         - right hand side vector
        reg_param - vector of regularization parameters, length N
        rho       - vector of residual norms ||y - X c||, length N
        eta       - vector of solution norms ||c||, length N
        kappa     - (output) vector of curvature values, length N
        work      - workspace

Return: success/error

Notes:
1) SVD of X must be computed first by calling multifit_linear_svd();
   work->n and work->p are initialized by this function

2) Vectors reg_param, rho, eta must be computed by calling gsl_multifit_linear_lcurve()
*/

int
gsl_multifit_linear_lcurvature (const gsl_vector * y,
                                const gsl_vector * reg_param,
                                const gsl_vector * rho,
                                const gsl_vector * eta,
                                gsl_vector * kappa,
                                gsl_multifit_linear_workspace * work)
{
  const size_t n = y->size;
  const size_t N = rho->size; /* number of points on L-curve */

  if (n != work->n)
    {
      GSL_ERROR("y vector does not match workspace", GSL_EBADLEN);
    }
  else if (N != eta->size)
    {
      GSL_ERROR ("size of rho and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else if (reg_param->size != N)
    {
      GSL_ERROR ("size of reg_param and rho vectors do not match",
                 GSL_EBADLEN);
    }
  else if (kappa->size != N)
    {
      GSL_ERROR ("size of reg_param and rho vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;
      const size_t p = work->p;
      gsl_matrix_view U = gsl_matrix_submatrix(work->A, 0, 0, n, p);
      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);
      gsl_vector_view beta = gsl_vector_subvector(work->xt, 0, p);
      size_t i;

      /* compute projection beta = U^T y */
      gsl_blas_dgemv (CblasTrans, 1.0, &U.matrix, y, 0.0, &beta.vector);

      for (i = 0; i < N; ++i)
        {
          double lambda = gsl_vector_get(reg_param, i);
          double lambda_sq = lambda * lambda;
          double eta_i = gsl_vector_get(eta, i);
          double rho_i = gsl_vector_get(rho, i);
          double phi_i = 0.0, dphi_i = 0.0;
          double psi_i = 0.0, dpsi_i = 0.0;
          double deta_i, ddeta_i, drho_i, ddrho_i;
          double dlogeta_i, ddlogeta_i, dlogrho_i, ddlogrho_i;
          double kappa_i;
          size_t j;

          for (j = 0; j < p; ++j)
            {
              double beta_j = gsl_vector_get(&beta.vector, j);
              double s_j = gsl_vector_get(&S.vector, j);
              double sj_sq = s_j * s_j;
              double f_j = sj_sq / (sj_sq + lambda_sq); /* filter factor */
              double onemf_j = 1.0 - f_j;
              double f1_j = -2.0 * f_j * onemf_j / lambda;
              double f2_j = -f1_j * (3.0 - 4.0 * f_j) / lambda;
              double xi_j = beta_j / s_j;

              phi_i += f_j * f1_j * xi_j * xi_j;
              psi_i += onemf_j * f1_j * beta_j * beta_j;

              dphi_i += (f1_j * f1_j + f_j * f2_j) * xi_j * xi_j;
              dpsi_i += (-f1_j * f1_j + onemf_j * f2_j) * beta_j * beta_j;
            }

          deta_i = phi_i / eta_i;
          drho_i = -psi_i / rho_i;
          ddeta_i = dphi_i / eta_i - deta_i * (deta_i / eta_i);
          ddrho_i = -dpsi_i / rho_i - drho_i * (drho_i / rho_i);

          dlogeta_i = deta_i / eta_i;
          dlogrho_i = drho_i / rho_i;
          ddlogeta_i = ddeta_i / eta_i - dlogeta_i * dlogeta_i;
          ddlogrho_i = ddrho_i / rho_i - dlogrho_i * dlogrho_i;

          kappa_i = (dlogrho_i * ddlogeta_i - ddlogrho_i * dlogeta_i) /
                    pow(dlogrho_i * dlogrho_i + dlogeta_i * dlogeta_i, 1.5);
          gsl_vector_set(kappa, i, kappa_i);
        }

      return status;
    }
}

/*
gsl_multifit_linear_lcurvature_menger()
  Determine Menger curvature of each point on L-curve. For each
set of 3 points on the L-curve, the circle which passes through
the 3 points is computed. The radius of the circle is then used
as an estimate of the curvature at the middle point (Menger curvature).

Inputs: rho   - vector of residual norms ||A x - b||
        eta   - vector of solution norms ||L x||
        kappa - (output) vector of Menger curvature values

Return: success/error
*/

int
gsl_multifit_linear_lcurvature_menger(const gsl_vector *rho,
                                      const gsl_vector *eta,
                                      gsl_vector *kappa)
{
  const size_t n = rho->size;

  if (n < 3)
    {
      GSL_ERROR ("at least 3 points are needed for Menger curvature",
                 GSL_EBADLEN);
    }
  else if (n != eta->size)
    {
      GSL_ERROR ("size of rho and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else if (n != kappa->size)
    {
      GSL_ERROR ("size of rho and kappa vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      double x1, y1;      /* first point of triangle on L-curve */
      double x2, y2;      /* second point of triangle on L-curve */
      size_t i;

      /* initial values */
      x1 = log(gsl_vector_get(rho, 0));
      y1 = log(gsl_vector_get(eta, 0));

      x2 = log(gsl_vector_get(rho, 1));
      y2 = log(gsl_vector_get(eta, 1));

      for (i = 1; i < n - 1; ++i)
        {
          /*
           * The points (x1,y1), (x2,y2), (x3,y3) are the previous,
           * current, and next point on the L-curve. We will find
           * the circle which fits these 3 points and take its radius
           * as an estimate of the curvature at this point.
           */
          double x3 = log(gsl_vector_get(rho, i + 1));
          double y3 = log(gsl_vector_get(eta, i + 1));

          double a = gsl_hypot(x1 - x2, y1 - y2); /* distances between each pair of points */
          double b = gsl_hypot(x2 - x3, y2 - y3);
          double c = gsl_hypot(x1 - x3, y1 - y3);

          double kappa_i;

          /*
           * sometimes the triangle sides a,b,c do not
           * form a valid triangle, but in that case
           * menger() will set kappa to 0; ignore the
           * error value in these cases
           */
          menger(a, b, c, &kappa_i);

          gsl_vector_set(kappa, i, kappa_i);

          /* update previous/current L-curve values */
          x1 = x2;
          y1 = y2;
          x2 = x3;
          y2 = y3;
        }

      /* duplicate the first and last points */
      gsl_vector_set(kappa, 0, gsl_vector_get(kappa, 1));
      gsl_vector_set(kappa, n - 1, gsl_vector_get(kappa, n - 2));

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_lcorner()
  Determine point on L-curve of maximum curvature. For each
set of 3 points on the L-curve, the circle which passes through
the 3 points is computed. The radius of the circle is then used
as an estimate of the curvature at the middle point. The point
with maximum curvature is then selected.

Inputs: rho - vector of residual norms ||A x - b||
        eta - vector of solution norms ||L x||
        idx - (output) index i such that
              (log(rho(i)),log(eta(i))) is the point of
              maximum curvature

Return: success/error
*/

int
gsl_multifit_linear_lcorner(const gsl_vector *rho,
                            const gsl_vector *eta,
                            size_t *idx)
{
  const size_t n = rho->size;

  if (n < 3)
    {
      GSL_ERROR ("at least 3 points are needed for L-curve analysis",
                 GSL_EBADLEN);
    }
  else if (n != eta->size)
    {
      GSL_ERROR ("size of rho and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int s = GSL_SUCCESS;
      double x1, y1;      /* first point of triangle on L-curve */
      double x2, y2;      /* second point of triangle on L-curve */
      double kappa_max = -1.0;
      size_t i;

      /* initial values */
      x1 = log(gsl_vector_get(rho, 0));
      y1 = log(gsl_vector_get(eta, 0));

      x2 = log(gsl_vector_get(rho, 1));
      y2 = log(gsl_vector_get(eta, 1));

      for (i = 1; i < n - 1; ++i)
        {
          /*
           * The points (x1,y1), (x2,y2), (x3,y3) are the previous,
           * current, and next point on the L-curve. We will find
           * the circle which fits these 3 points and take its radius
           * as an estimate of the curvature at this point.
           */
          double x3 = log(gsl_vector_get(rho, i + 1));
          double y3 = log(gsl_vector_get(eta, i + 1));

          double a = gsl_hypot(x1 - x2, y1 - y2); /* distances between each pair of points */
          double b = gsl_hypot(x2 - x3, y2 - y3);
          double c = gsl_hypot(x1 - x3, y1 - y3);

          double kappa_i;

          s = menger(a, b, c, &kappa_i);
          if (s)
            continue;

          if (kappa_max < kappa_i)
            {
              kappa_max = kappa_i;
              *idx = i;
            }

          /* update previous/current L-curve values */
          x1 = x2;
          y1 = y2;
          x2 = x3;
          y2 = y3;
        }

      /* check if a maximum curvature was found */
      if (kappa_max <= 0.0)
        {
          /* possibly co-linear points */
          GSL_ERROR("failed to find maximum curvature", GSL_EINVAL);
        }

      return s;
    }
}

/*
gsl_multifit_linear_lcorner2()
  Determine point on L-curve (lambda^2, ||c||^2) of maximum curvature.
For each set of 3 points on the L-curve, the circle which passes through
the 3 points is computed. The radius of the circle is then used
as an estimate of the curvature at the middle point. The point
with maximum curvature is then selected.

This routine is based on the paper

M. Rezghi and S. M. Hosseini, "A new variant of L-curve for Tikhonov
regularization", J. Comp. App. Math., 231 (2009).

Inputs: reg_param - vector of regularization parameters
        eta       - vector of solution norms ||L x||
        idx       - (output) index i such that
                    (lambda(i)^2,eta(i)^2) is the point of
                    maximum curvature

Return: success/error
*/

int
gsl_multifit_linear_lcorner2(const gsl_vector *reg_param,
                             const gsl_vector *eta,
                             size_t *idx)
{
  const size_t n = reg_param->size;

  if (n < 3)
    {
      GSL_ERROR ("at least 3 points are needed for L-curve analysis",
                 GSL_EBADLEN);
    }
  else if (n != eta->size)
    {
      GSL_ERROR ("size of reg_param and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int s = GSL_SUCCESS;
      size_t i;
      double x1, y1;      /* first point of triangle on L-curve */
      double x2, y2;      /* second point of triangle on L-curve */
      double rmin = -1.0; /* minimum radius of curvature */

      /* initial values */
      x1 = gsl_vector_get(reg_param, 0);
      x1 *= x1;
      y1 = gsl_vector_get(eta, 0);
      y1 *= y1;

      x2 = gsl_vector_get(reg_param, 1);
      x2 *= x2;
      y2 = gsl_vector_get(eta, 1);
      y2 *= y2;

      for (i = 1; i < n - 1; ++i)
        {
          /*
           * The points (x1,y1), (x2,y2), (x3,y3) are the previous,
           * current, and next point on the L-curve. We will find
           * the circle which fits these 3 points and take its radius
           * as an estimate of the curvature at this point.
           */
          double lamip1 = gsl_vector_get(reg_param, i + 1);
          double etaip1 = gsl_vector_get(eta, i + 1);
          double x3 = lamip1 * lamip1;
          double y3 = etaip1 * etaip1;

          double x21 = x2 - x1;
          double y21 = y2 - y1;
          double x31 = x3 - x1;
          double y31 = y3 - y1;
          double h21 = x21*x21 + y21*y21;
          double h31 = x31*x31 + y31*y31;
          double d = fabs(2.0 * (x21*y31 - x31*y21));
          double r = sqrt(h21*h31*((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))) / d;

          /* if d =~ 0 then there are nearly colinear points */
          if (gsl_finite(r))
            {
              /* check for smallest radius of curvature */
              if (r < rmin || rmin < 0.0)
                {
                  rmin = r;
                  *idx = i;
                }
            }

          /* update previous/current L-curve values */
          x1 = x2;
          y1 = y2;
          x2 = x3;
          y2 = y3;
        }

      /* check if a minimum radius was found */
      if (rmin < 0.0)
        {
          /* possibly co-linear points */
          GSL_ERROR("failed to find minimum radius", GSL_EINVAL);
        }

      return s;
    }
} /* gsl_multifit_linear_lcorner2() */

#define GSL_MULTIFIT_MAXK      100

/*
gsl_multifit_linear_L()
  Compute discrete approximation to derivative operator of order
k on a regular grid of p points, ie: L is (p-k)-by-p
*/

int
gsl_multifit_linear_Lk(const size_t p, const size_t k, gsl_matrix *L)
{
  if (p <= k)
    {
      GSL_ERROR("p must be larger than derivative order", GSL_EBADLEN);
    }
  else if (k >= GSL_MULTIFIT_MAXK - 1)
    {
      GSL_ERROR("derivative order k too large", GSL_EBADLEN);
    }
  else if (p - k != L->size1 || p != L->size2)
    {
      GSL_ERROR("L matrix must be (p-k)-by-p", GSL_EBADLEN);
    }
  else
    {
      double c_data[GSL_MULTIFIT_MAXK];
      gsl_vector_view cv = gsl_vector_view_array(c_data, k + 1);
      size_t i, j;

      /* zeroth derivative */
      if (k == 0)
        {
          gsl_matrix_set_identity(L);
          return GSL_SUCCESS;
        }

      gsl_matrix_set_zero(L);
  
      gsl_vector_set_zero(&cv.vector);
      gsl_vector_set(&cv.vector, 0, -1.0);
      gsl_vector_set(&cv.vector, 1, 1.0);

      for (i = 1; i < k; ++i)
        {
          double cjm1 = 0.0;

          for (j = 0; j < k + 1; ++j)
            {
              double cj = gsl_vector_get(&cv.vector, j);

              gsl_vector_set(&cv.vector, j, cjm1 - cj);
              cjm1 = cj;
            }
        }

      /* build L, the c_i are the entries on the diagonals */
      for (i = 0; i < k + 1; ++i)
        {
          gsl_vector_view v = gsl_matrix_superdiagonal(L, i);
          double ci = gsl_vector_get(&cv.vector, i);

          gsl_vector_set_all(&v.vector, ci);
        }

      return GSL_SUCCESS;
    }
} /* gsl_multifit_linear_Lk() */

/*
gsl_multifit_linear_Lsobolev()
  Construct Sobolev smoothing norm operator

L = [ a_0 I; a_1 L_1; a_2 L_2; ...; a_k L_k ]

by computing the Cholesky factor of L^T L

Inputs: p     - number of columns of L
        kmax  - maximum derivative order (< p)
        alpha - vector of weights; alpha_k multiplies L_k, size kmax + 1
        L     - (output) upper triangular Sobolev matrix p-by-p,
                stored in upper triangle
        work  - workspace

Notes:
1) work->Q is used to store intermediate L_k matrices
*/

int
gsl_multifit_linear_Lsobolev(const size_t p, const size_t kmax,
                             const gsl_vector *alpha, gsl_matrix *L,
                             gsl_multifit_linear_workspace *work)
{
  if (p > work->pmax)
    {
      GSL_ERROR("p is larger than workspace", GSL_EBADLEN);
    }
  else if (p <= kmax)
    {
      GSL_ERROR("p must be larger than derivative order", GSL_EBADLEN);
    }
  else if (kmax + 1 != alpha->size)
    {
      GSL_ERROR("alpha must be size kmax + 1", GSL_EBADLEN);
    }
  else if (p != L->size1)
    {
      GSL_ERROR("L matrix is wrong size", GSL_EBADLEN);
    }
  else if (L->size1 != L->size2)
    {
      GSL_ERROR("L matrix is not square", GSL_ENOTSQR);
    }
  else
    {
      int s;
      size_t j, k;
      gsl_vector_view d = gsl_matrix_diagonal(L);
      const double alpha0 = gsl_vector_get(alpha, 0);

      /* initialize L to alpha0^2 I */
      gsl_matrix_set_zero(L);
      gsl_vector_add_constant(&d.vector, alpha0 * alpha0);

      for (k = 1; k <= kmax; ++k)
        {
          gsl_matrix_view Lk = gsl_matrix_submatrix(work->Q, 0, 0, p - k, p);
          double ak = gsl_vector_get(alpha, k);

          /* compute a_k L_k */
          s = gsl_multifit_linear_Lk(p, k, &Lk.matrix);
          if (s)
            return s;

          gsl_matrix_scale(&Lk.matrix, ak);

          /* LTL += L_k^T L_k */
          gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Lk.matrix, 1.0, L);
        }

      s = gsl_linalg_cholesky_decomp(L);
      if (s)
        return s;

      /* copy Cholesky factor to upper triangle and zero out bottom */
      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, L, L);

      for (j = 0; j < p; ++j)
        {
          for (k = 0; k < j; ++k)
            gsl_matrix_set(L, j, k, 0.0);
        }

      return GSL_SUCCESS;
    }
}

/*
menger()
  Calculate Menger curvature, which is the reciprocal radius
of a circle fitted to three points in the plane (x,y,z)

Inputs: a,b,c - Euclidean distance between the 3 points in the plane,
                a = | x - y |
                b = | x - z |
                c = | y - z |
        kappa - (output) Menger curvature

Return: if (a,b,c) are negative or represent colinear points, GSL_EDOM is returned;
otherwise GSL_SUCCESS

Notes:
1) See https://www.wikipedia.org/wiki/Menger_curvature for formula
*/

static int
menger(const double a, const double b, const double c, double * kappa)
{
  if (a <= 0.0 || b <= 0.0 || c <= 0.0)
    {
      /* triangle sides must be positive;
       * just set kappa to 0 and return */
      *kappa = 0.0;
      return GSL_EDOM;
    }
  else if (b + c < a || c + a < b || a + b < c)
    {
      /* the given sides do not form a valid triangle;
       * just set kappa to 0 and return */
      *kappa = 0.0;
      return GSL_EDOM;
    }
  else
    {
      double min, middle, max;
      double A;

      /* sort a,b,c */
      if (a < b)
        {
          if (a < c)
            {
              min = a;
              if (b < c)
                {
                  middle = b;
                  max = c;
                }
              else
                {
                  middle = c;
                  max = b;
                }
            }
          else
            {
              min = c;
              middle = a;
              max = b;
            }
        }
      else
        {
          if (b < c)
            {
              min = b;
              if (a < c)
                {
                  middle = a;
                  max = c;
                }
              else
                {
                  middle = c;
                  max = a;
                }
            }
          else
            {
              min = c;
              middle = b;
              max = a;
            }
        }

      /* check for quick return */
      if (min + middle == max)
        {
          /* colinear points */
          *kappa = 0.0;
          return GSL_SUCCESS;
        }

      /*
       * Compute area of triangle with sides a,b,c.
       * see: http://www.cs.berkeley.edu/~wkahan/Triangle.pdf
       * for a numerically stable formula
       */
      A = sqrt((max + (middle + min)) *
               (min - (max - middle)) *
               (min + (max - middle)) *
               (max + (middle - min)));

      *kappa = A / (a * b * c);

      return GSL_SUCCESS;
    }
}
