/* bspline/eval.c
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

/* routines related to evaluation of B-splines and their derivatives */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

static int bspline_process_interval_for_eval (const double x, size_t * i, const int flag,
                                              gsl_bspline_workspace * w);

static void bspline_pppack_bsplvb (const gsl_vector * t,
		                               const size_t jhigh,
		                               const size_t index,
		                               const double x,
		                               const size_t left,
		                               size_t * j,
		                               gsl_vector * deltal,
		                               gsl_vector * deltar, gsl_vector * biatx);

static void bspline_pppack_bsplvd (const gsl_vector * t,
		                               const size_t k,
		                               const double x,
		                               const size_t left,
		                               gsl_vector * deltal,
		                               gsl_vector * deltar,
		                               gsl_matrix * a,
		                               gsl_matrix * dbiatx, const size_t nderiv);

static int bspline_pppack_bvalue(const double x,
                                 const gsl_vector * c,
                                 const size_t nderiv,
                                 double * result,
                                 gsl_bspline_workspace * w);

/*
gsl_bspline_calc()
  Evalute B-spline at a given point

Inputs: x      - point for evaluation
        c      - coefficient vector, size ncontrol
        result - (output) B-spline value at x
        w      - workspace

Return: success/error
*/

int
gsl_bspline_calc(const double x, const gsl_vector * c,
                 double * result, gsl_bspline_workspace * w)
{
  int status = gsl_bspline_calc_deriv(x, c, 0, result, w);
  return status;
}

/*
gsl_bspline_calc_deriv()
  Evalute derivative of B-spline at a given point

Inputs: x      - point for evaluation
        c      - coefficient vector, size ncontrol
        nderiv - derivative order
        result - (output) B-spline value at x
        w      - workspace

Return: success/error

Notes:
1) Adapted from PPPACK routine bvalue
*/

int
gsl_bspline_calc_deriv(const double x, const gsl_vector * c,
                       const size_t nderiv, double * result,
                       gsl_bspline_workspace * w)
{
  if (c->size != w->ncontrol)
    {
      GSL_ERROR("coefficient vector does not match workspace", GSL_EBADLEN);
    }
  else if (nderiv >= w->spline_order)
    {
      /* quick return */
      *result = 0.0;
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const double t0 = gsl_vector_get(w->knots, 0);
      const double t1 = gsl_vector_get(w->knots, w->knots->size - 1);
      double *coef = w->B->data;          /* coefficients of Taylor expansion around t */
      double *deriv_result = w->dB->data; /* array of derivatives of Taylor expansion */

      if (x < t0 || x > t1)
        {
          /*
           * extrapolate for x < t_0 or x > t_1 by computing coefficients of a
           * Taylor series about the endpoint (t_0 or t_1) and using gsl_poly
           * to evaluate series
           */

          double fac = 1.0;
          double t = (x < t0) ? t0 : t1;
          double h = x - t;
          size_t i;

          /* FIXME: this is not efficient */
          for (i = 0; i < w->spline_order; ++i)
            {
              status = bspline_pppack_bvalue(t, c, i, &coef[i], w);
              if (status)
                return status;

              /* divide coef[i] by i! */
              coef[i] *= fac;
              fac /= (i + 1.0);
            }

          gsl_poly_eval_derivs(coef, w->spline_order, h, deriv_result, nderiv + 1);
          *result = deriv_result[nderiv];
        }
      else
        {
          /* t_0 <= x <= t_{n+k} */

          status = bspline_pppack_bvalue(x, c, nderiv, result, w);
          if (status)
            return status;
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_vector_calc()
  Evalute vector-valued B-spline at a given point

result = sum_{i=1}^{ncontrol} c_i B_{i,k}(x)

where c_i are vectors of length nvec

Inputs: x      - point for evaluation
        c      - spline coefficients, size nvec-by-ncontrol
        result - (output) B-spline value at x, size nvec
        w      - workspace

Return: success/error
*/

int
gsl_bspline_vector_calc(const double x, const gsl_matrix * c,
                        gsl_vector * result, gsl_bspline_workspace * w)
{
  int status = gsl_bspline_vector_calc_deriv(x, c, 0, result, w);
  return status;
}

/*
gsl_bspline_vector_calc_deriv()
  Evalute derivative of vector-valued B-spline at a given point

result = sum_{i=1}^{ncontrol} c_i d^(nderiv)/dx^(nderiv) B_{i,k}(x)

where c_i are vectors of length nvec

Inputs: x      - point for evaluation
        c      - spline coefficients, size nvec-by-ncontrol
        nderiv - derivative order
        result - (output) B-spline value at x, size nvec
        w      - workspace

Return: success/error

Notes:
1) Adapted from PPPACK routine bvalue
*/

int
gsl_bspline_vector_calc_deriv(const double x,
                              const gsl_matrix * c,
                              const size_t nderiv,
                              gsl_vector * result,
                              gsl_bspline_workspace * w)
{
  if (c->size1 != result->size)
    {
      GSL_ERROR("coefficient matrix does not match result vector", GSL_EBADLEN);
    }
  else if (c->size2 != w->ncontrol)
    {
      GSL_ERROR("coefficient matrix does not match workspace", GSL_EBADLEN);
    }
  else if (nderiv >= w->spline_order)
    {
      /* quick return */
      gsl_vector_set_zero(result);
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      size_t i;

      /* XXX: this is not efficient */
      for (i = 0; i < c->size1; ++i)
        {
          gsl_vector_const_view coef = gsl_matrix_const_row(c, i);
          double *resulti = gsl_vector_ptr(result, i);

          status = gsl_bspline_calc_deriv(x, &coef.vector, nderiv, resulti, w);
          if (status)
            return status;
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_eval_basis()
  Evaluate the basis functions B_i(x) for all i. This is
a wrapper function for gsl_bspline_basis() which
formats the output in a nice way.

Inputs: x - point for evaluation
        B - (output) where to store B_i(x) values
            the length of this vector is
            n = nbreak + k - 2 = l + k - 1 = w->ncontrol
        w - bspline workspace

Return: success or error

Notes: The w->knots vector must be initialized prior to calling
       this function (see gsl_bspline_init_*)
*/

int
gsl_bspline_eval_basis (const double x,
                        gsl_vector * B,
                        gsl_bspline_workspace * w)
{
  if (B->size != w->ncontrol)
    {
      GSL_ERROR ("size of B does not match workspace", GSL_EBADLEN);
    }
  else
    {
      const size_t k = w->spline_order;
      size_t i;      /* looping */
      size_t istart; /* first non-zero spline for x */
      int error;     /* error handling */

      /* find all non-zero B_i(x) values */
      error = gsl_bspline_basis (x, w->B, &istart, w);
      if (error)
        return error;

      /* store values in appropriate part of given vector */
      for (i = 0; i < istart; i++)
        gsl_vector_set (B, i, 0.0);

      for (i = istart; i < istart + k; i++)
        gsl_vector_set (B, i, gsl_vector_get (w->B, i - istart));

      for (i = istart + k; i < w->ncontrol; i++)
        gsl_vector_set (B, i, 0.0);

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_basis()
  Evaluate all non-zero B-spline functions at point x.
These are the B_i(x) for i in [istart, istart + k - 1].
Always B_i(x) = 0 for i < istart and for i > istart + k - 1.

Inputs: x      - point at which to evaluate splines
        Bk     - (output) where to store B-spline values (length spline_order)
        istart - (output) B-spline function index of
                 first non-zero basis for given x
        w      - bspline workspace

Return: success or error

Notes: 1) the w->knots vector must be initialized before calling
          this function

       2) On output, B contains

             [B_{istart,k}, B_{istart+1,k}, ..., B_{istart+k-1,k}]

          evaluated at the given x.
*/

int
gsl_bspline_basis (const double x, gsl_vector * Bk, size_t * istart,
                   gsl_bspline_workspace * w)
{
  if (Bk->size != w->spline_order)
    {
      GSL_ERROR ("Bk vector length does not match spline order", GSL_EBADLEN);
    }
  else
    {
      size_t left, j;
      int flag = 0;   /* interval search flag */
      int error = 0;  /* error flag */

      left = gsl_bspline_find_interval (x, &flag, w);
      error = bspline_process_interval_for_eval (x, &left, flag, w);
      if (error)
        return error;

      bspline_pppack_bsplvb (w->knots, w->spline_order, 1, x, left, &j, w->deltal,
                             w->deltar, Bk);

      if (left >= w->ncontrol)
        left = w->ncontrol - 1;

      *istart = left - w->spline_order + 1;

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_eval_deriv_basis()
  Evaluate d^j/dx^j B_i(x) for all i, 0 <= j <= nderiv.
This is a wrapper function for gsl_bspline_basis_deriv()
which formats the output in a nice way.

Inputs: x      - point for evaluation
        nderiv - number of derivatives to compute, inclusive.
        dB     - (output) where to store d^j/dx^j B_i(x)
                 values. the size of this matrix is
                 w->ncontrol by (nderiv + 1)
        w      - bspline derivative workspace

Return: success or error

Notes: 1) The w->knots vector must be initialized prior to calling
          this function (see gsl_bspline_init_*)

       2) based on PPPACK's bsplvd
*/

int
gsl_bspline_eval_deriv_basis (const double x, const size_t nderiv,
                              gsl_matrix * dB, gsl_bspline_workspace * w)
{
  if (dB->size1 != w->ncontrol)
    {
      GSL_ERROR ("dB matrix first dimension does not match workspace", GSL_EBADLEN);
    }
  else if (dB->size2 < nderiv + 1)
    {
      GSL_ERROR ("dB matrix second dimension must be at least length nderiv+1", GSL_EBADLEN);
    }
  else
    {
      const size_t k = w->spline_order;
      size_t i, j;   /* looping */
      size_t istart; /* first non-zero spline for x */
      int error;     /* error handling */

      /* find all non-zero d^j/dx^j B_i(x) values */
      error = gsl_bspline_basis_deriv (x, nderiv, w->dB, &istart, w);
      if (error)
        return error;

      /* store values in appropriate part of given matrix */
      for (j = 0; j <= nderiv; j++)
        {
          for (i = 0; i < istart; i++)
          gsl_matrix_set (dB, i, j, 0.0);

          for (i = istart; i < istart + k; i++)
            gsl_matrix_set (dB, i, j, gsl_matrix_get (w->dB, i - istart, j));

          for (i = istart + k; i < w->ncontrol; i++)
            gsl_matrix_set (dB, i, j, 0.0);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_basis_deriv()
  At point x evaluate all requested, non-zero B-spline function
derivatives and store them in dB.  These are the
d^j/dx^j B_i(x) with i in [istart, istart + k - 1] and j in [0, nderiv].
Always d^j/dx^j B_i(x) = 0 for i < istart and for i > istart + k - 1.

Inputs: x      - point at which to evaluate splines
        nderiv - number of derivatives to request, inclusive
        dB     - (output) where to store dB-spline derivatives
                 (size k by nderiv + 1)
        istart - (output) B-spline function index of
                 first non-zero basis for given x
        w      - bspline derivative workspace

Return: success or error

Notes: 1) the w->knots vector must be initialized before calling
          this function

       2) On output, dB contains

          dB(1,:)   = [ B_{istart,k},       d/dx B_{istart,k},         ..., d^n/dx^n B_{istart,k}   
          dB(2,:)   = [ B_{istart+1,k},     d/dx B_{istart+1,k},       ..., d^n/dx^n B_{istart+1,k} ]
          ...
          dB(k-1,:) = [ B_{istart+k-2,k},   d/dx B_{istart+k-2,k},     ..., d^n/dx^n B_{istart+k-2,k}   ]
          dB(k,:)   = [ B_{istart+k-1,k},   d/dx B_{istart+k-1,k},     ..., d^n/dx^n B_{istart+k-1,k}   ]

          evaluated at x.  B_{istart, k} is stored in dB(0,0).
          Each additional column contains an additional derivative.
          Above, k is the spline order and n is nderiv.

       3) Note that the zero-th column of the result contains the
          0th derivative, which is simply a function evaluation.

       4) based on PPPACK's bsplvd
*/

int
gsl_bspline_basis_deriv (const double x, const size_t nderiv,
                         gsl_matrix * dB, size_t * istart,
                         gsl_bspline_workspace * w)
{
  if (dB->size1 != w->spline_order)
    {
      GSL_ERROR ("dB matrix first dimension does not match workspace", GSL_EBADLEN);
    }
  else if (dB->size2 < nderiv + 1)
    {
      GSL_ERROR ("dB matrix second dimension must be at least length nderiv+1", GSL_EBADLEN);
    }
  else
    {
      size_t left;   /* spline index */
      int flag = 0;  /* interval search flag */
      int error = 0; /* error flag */

      /* check for quick return */
      if (nderiv == 0)
        {
          gsl_vector_view B = gsl_matrix_column(dB, 0);
          return gsl_bspline_basis(x, &B.vector, istart, w);
        }

      left = gsl_bspline_find_interval (x, &flag, w);
      error = bspline_process_interval_for_eval (x, &left, flag, w);
      if (error)
        return error;

      *istart = left - w->spline_order + 1;

      bspline_pppack_bsplvd (w->knots, w->spline_order, x, left,
                             w->deltal, w->deltar, w->A, dB, nderiv);

      /* if they provide columns for derivatives of order k or more, zero them
       * out; the highest order non-zero derivative for a B-spline of order k
       * is nderiv_max = k - 1 */
      if (dB->size2 > w->spline_order)
        {
          gsl_matrix_view m = gsl_matrix_submatrix(dB, 0, w->spline_order,
                                                   w->spline_order,
                                                   dB->size2 - w->spline_order);
          gsl_matrix_set_zero(&m.matrix);
        }

      return GSL_SUCCESS;
    }
}

/*
bspline_process_interval_for_eval()
  Consumes an x location, left knot from gsl_bspline_find_interval, flag
from gsl_bspline_find_interval, and a workspace.  Checks that x lies within
the splines' knots, enforces some endpoint continuity requirements, and
avoids divide by zero errors in the underlying bspline_pppack_* functions.
*/
static int
bspline_process_interval_for_eval (const double x, size_t * i, const int flag,
                                   gsl_bspline_workspace * w)
{
  if (flag == -1)
    {
      GSL_ERROR ("x outside of knot interval", GSL_EINVAL);
    }
  else if (flag == 1)
    {
      if (x <= gsl_vector_get (w->knots, *i) + GSL_DBL_EPSILON)
        {
          *i -= 1;
        }
      else
        {
          GSL_ERROR ("x outside of knot interval", GSL_EINVAL);
        }
    }

  if (gsl_vector_get (w->knots, *i) == gsl_vector_get (w->knots, *i + 1))
    {
      GSL_ERROR ("knot(i) = knot(i+1) will result in division by zero", GSL_EINVAL);
    }

  return GSL_SUCCESS;
}

/********************************************************************
 * PPPACK ROUTINES
 *
 * The routines herein deliberately avoid using the bspline workspace,
 * choosing instead to pass all work areas explicitly.  This allows
 * others to more easily adapt these routines to low memory or
 * parallel scenarios.
 ********************************************************************/

/*
bspline_pppack_bsplvb()
  calculates the value of all possibly nonzero b-splines at x of order
jout = max( jhigh , (j+1)*(index-1) ) with knot sequence t.

Parameters:
   t      - knot sequence, of length left + jout , assumed to be
            nondecreasing.  assumption t(left).lt.t(left + 1).
            division by zero  will result if t(left) = t(left+1)
   jhigh  -
   index  - integers which determine the order jout = max(jhigh,
            (j+1)*(index-1))  of the b-splines whose values at x
            are to be returned.  index  is used to avoid
            recalculations when several columns of the triangular
            array of b-spline values are needed (e.g., in  bsplpp
            or in  bsplvd ).  precisely,

            if  index = 1 ,
               the calculation starts from scratch and the entire
               triangular array of b-spline values of orders
               1,2,...,jhigh  is generated order by order , i.e.,
               column by column .

            if  index = 2 ,
               only the b-spline values of order j+1, j+2, ..., jout
               are generated, the assumption being that biatx, j,
               deltal, deltar are, on entry, as they were on exit
               at the previous call.

            in particular, if jhigh = 0, then jout = j+1, i.e.,
            just the next column of b-spline values is generated.
   x      - the point at which the b-splines are to be evaluated.
   left   - an integer chosen (usually) so that
            t(left) .le. x .le. t(left+1).
   j      - (output) a working scalar for indexing
   deltal - (output) a working area which must be of length at least jout
   deltar - (output) a working area which must be of length at least jout
   biatx  - (output) array of length jout, with  biatx(i)
            containing the value at  x  of the polynomial of order
            jout which agrees with the b-spline b(left-jout+i,jout,t)
            on the interval (t(left), t(left+1)) .

Method:
   the recurrence relation

                      x - t(i)              t(i+j+1) - x
      b(i,j+1)(x) = -----------b(i,j)(x) + ---------------b(i+1,j)(x)
                    t(i+j)-t(i)            t(i+j+1)-t(i+1)

   is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
   ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
   b(left,j)(x), storing the new values in  biatx  over the old. the
   facts that

      b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)

   and that

      b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)

   are used. the particular organization of the calculations follows
   algorithm (8) in chapter x of [1].

Notes:

   (1) This is a direct translation of PPPACK's bsplvb routine with
       j, deltal, deltar rewritten as input parameters and
       utilizing zero-based indexing.

   (2) This routine contains no error checking.  Please use routines
       like gsl_bspline_basis().
*/

static void
bspline_pppack_bsplvb (const gsl_vector * t,
                       const size_t jhigh,
                       const size_t index,
                       const double x,
                       const size_t left,
                       size_t * j,
                       gsl_vector * deltal,
                       gsl_vector * deltar, gsl_vector * biatx)
{
  size_t i;			/* looping */
  double saved;
  double term;

  if (index == 1)
    {
      *j = 0;
      gsl_vector_set (biatx, 0, 1.0);
    }

  for ( /* NOP */ ; *j < jhigh - 1; *j += 1)
    {
      gsl_vector_set (deltar, *j, gsl_vector_get (t, left + *j + 1) - x);
      gsl_vector_set (deltal, *j, x - gsl_vector_get (t, left - *j));

      saved = 0.0;

      for (i = 0; i <= *j; i++)
        {
          term = gsl_vector_get (biatx, i) / (gsl_vector_get (deltar, i)
                 + gsl_vector_get (deltal, *j - i));

          gsl_vector_set (biatx, i, saved + gsl_vector_get (deltar, i) * term);

          saved = gsl_vector_get (deltal, *j - i) * term;
        }

      gsl_vector_set (biatx, *j + 1, saved);
    }

  return;
}

/*
bspline_pppack_bsplvd()
  calculates value and derivs of all b-splines which do not vanish at x

Parameters:
   t      - the knot array, of length left+k (at least)
   k      - the order of the b-splines to be evaluated
   x      - the point at which these values are sought
   left   - an integer indicating the left endpoint of the interval
            of interest. the k b-splines whose support contains the
            interval (t(left), t(left+1)) are to be considered.
            it is assumed that t(left) .lt. t(left+1)
            division by zero will result otherwise (in  bsplvb).
            also, the output is as advertised only if
            t(left) .le. x .le. t(left+1) .
   deltal - a working area which must be of length at least k
   deltar - a working area which must be of length at least k
   a      - an array of order (k,k), to contain b-coeffs of the
            derivatives of a certain order of the k b-splines
            of interest.
   dbiatx - an array of order (k,nderiv). its entry (i,m) contains
            value of (m)th derivative of (left-k+i)-th b-spline
            of order k for knot sequence  t, i=1,...,k, m=0,...,nderiv.
   nderiv - an integer indicating that values of b-splines and
            their derivatives up to AND INCLUDING the nderiv-th
            are asked for. (nderiv is replaced internally by the
            integer mhigh in (1,k) closest to it.)

Method:
   values at x of all the relevant b-splines of order k,k-1,..., k+1-nderiv
   are generated via bsplvb and stored temporarily in dbiatx.  then, the
   b-coeffs of the required derivatives of the b-splines of interest are
   generated by differencing, each from the preceeding one of lower order,
   and combined with the values of b-splines of corresponding order in
   dbiatx  to produce the desired values .

Notes:

   (1) This is a direct translation of PPPACK's bsplvd routine with
       deltal, deltar rewritten as input parameters (to later feed them
       to bspline_pppack_bsplvb) and utilizing zero-based indexing.

   (2) This routine contains no error checking.
*/

static void
bspline_pppack_bsplvd (const gsl_vector * t,
                       const size_t k,
                       const double x,
                       const size_t left,
                       gsl_vector * deltal,
                       gsl_vector * deltar,
                       gsl_matrix * a,
                       gsl_matrix * dbiatx, const size_t nderiv)
{
  int i, ideriv, il, j, jlow, jp1mid, kmm, ldummy, m, mhigh;
  double factor, fkmm, sum;

  size_t bsplvb_j;
  gsl_vector_view dbcol = gsl_matrix_column (dbiatx, 0);

  mhigh = GSL_MIN_INT (nderiv, k - 1);
  bspline_pppack_bsplvb (t, k - mhigh, 1, x, left, &bsplvb_j, deltal, deltar,
                         &dbcol.vector);

  if (mhigh > 0)
    {
      /* the first column of dbiatx always contains the b-spline
         values for the current order. these are stored in column
         k-current order before bsplvb is called to put values
         for the next higher order on top of it.  */
      ideriv = mhigh;
      for (m = 1; m <= mhigh; m++)
        {
          for (j = ideriv, jp1mid = 0; j < (int) k; j++, jp1mid++)
            {
              gsl_matrix_set (dbiatx, j, ideriv, gsl_matrix_get (dbiatx, jp1mid, 0));
            }

          ideriv--;
          bspline_pppack_bsplvb (t, k - ideriv, 2, x, left, &bsplvb_j, deltal,
                                 deltar, &dbcol.vector);
        }

      /* at this point,  b(left-k+i, k+1-j)(x) is in  dbiatx(i,j)
         for i=j,...,k-1 and j=0,...,mhigh. in particular, the
         first column of dbiatx is already in final form. to obtain
         corresponding derivatives of b-splines in subsequent columns,
         generate their b-repr. by differencing, then evaluate at x. */
      jlow = 0;
      for (i = 0; i < (int) k; i++)
        {
          for (j = jlow; j < (int) k; j++)
            {
              gsl_matrix_set (a, j, i, 0.0);
            }
          jlow = i;
          gsl_matrix_set (a, i, i, 1.0);
        }

      /* at this point, a(.,j) contains the b-coeffs for the j-th of the
         k b-splines of interest here. */
      for (m = 1; m <= mhigh; m++)
        {
          kmm = k - m;
          fkmm = (float) kmm;
          il = left;
          i = k - 1;

          /* for j=1,...,k, construct b-coeffs of (m)th  derivative
             of b-splines from those for preceding derivative by
             differencing and store again in  a(.,j) . the fact that
             a(i,j) = 0  for i .lt. j  is used. */
          for (ldummy = 0; ldummy < kmm; ldummy++)
            {
              factor = fkmm / (gsl_vector_get (t, il + kmm) -
                       gsl_vector_get (t, il));

              /* the assumption that t(left).lt.t(left+1) makes
                 denominator in factor nonzero. */
              for (j = 0; j <= i; j++)
                {
                  gsl_matrix_set (a, i, j,
                                  factor * (gsl_matrix_get (a, i, j)
                                  - gsl_matrix_get (a, i - 1, j)));
                }
              il--;
              i--;
            }

          /* for i=1,...,k, combine b-coeffs a(.,i) with b-spline values
             stored in dbiatx(.,m) to get value of (m)th  derivative
             of i-th b-spline (of interest here) at x, and store in
             dbiatx(i,m). storage of this value over the value of a
             b-spline of order m there is safe since the remaining
             b-spline derivatives of the same order do not use this
             value due to the fact that a(j,i) = 0 for j .lt. i . */
          for (i = 0; i < (int) k; i++)
            {
              sum = 0;
              jlow = GSL_MAX_INT (i, m);
              for (j = jlow; j < (int) k; j++)
                {
                  sum += gsl_matrix_get (a, j, i) * gsl_matrix_get (dbiatx, j, m);
                }
              gsl_matrix_set (dbiatx, i, m, sum);
            }
        }
    }

  return;
}

/*
bspline_pppack_bvalue()
  Calculate value of B-spline derivative at a point x with a coefficient
vector c

Inputs: x      - point at which to evaluate B-spline
        c      - coefficient vector, length w->ncontrol
        nderiv - derivative order to compute
        result - (output) d^nderiv/dx^nderiv B(x)
        w      - workspace

Return: success/error
*/

static int
bspline_pppack_bvalue(const double x, const gsl_vector * c, const size_t nderiv,
                      double * result, gsl_bspline_workspace * w)
{
  const int n = (int) w->ncontrol;
  const int k = (int) w->spline_order;
  const int km1 = k - 1;
  gsl_vector_view a = gsl_matrix_row(w->A, 0);
  int flag;
  int jcmin, jcmax, imk, nmi;
  int i, j;

  *result = 0.0;

  /* check for quick return */
  if (nderiv >= w->spline_order)
    return GSL_SUCCESS;

  i = (int) gsl_bspline_find_interval (x, &flag, w);
  if (flag != 0)
    return GSL_SUCCESS;

  if (k == 1)
    {
      *result = gsl_vector_get(c, i);
      return GSL_SUCCESS;
    }

  imk = i - k + 1;
  if (imk < 0)
    {
      jcmin = -imk;

      for (j = 0; j <= i; ++j)
        gsl_vector_set(w->deltal, j, x - gsl_vector_get(w->knots, i - j));

      for (j = i; j < km1; ++j)
        {
          gsl_vector_set(&a.vector, k - j - 2, 0.0);
          gsl_vector_set(w->deltal, j, gsl_vector_get(w->deltal, i));
        }
    }
  else
    {
      jcmin = 0;

      for (j = 0; j < km1; ++j)
        gsl_vector_set(w->deltal, j, x - gsl_vector_get(w->knots, i - j));
    }

  nmi = n - i - 1;
  if (nmi < 0)
    {
      jcmax = k + nmi - 1;

      for (j = 0; j <= jcmax; ++j)
        gsl_vector_set(w->deltar, j, gsl_vector_get(w->knots, i + j + 1) - x);

      for (j = jcmax; j < km1; ++j)
        {
          gsl_vector_set(&a.vector, j + 1, 0.0);
          gsl_vector_set(w->deltar, j, gsl_vector_get(w->deltar, jcmax));
        }
    }
  else
    {
      jcmax = km1;

      for (j = 0; j < km1; ++j)
        gsl_vector_set(w->deltar, j, gsl_vector_get(w->knots, i + j + 1) - x);
    }

  for (j = jcmin; j <= jcmax; ++j)
    gsl_vector_set(&a.vector, j, gsl_vector_get(c, imk + j));

  if (nderiv != 0)
    {
      for (j = 0; j < (int) nderiv; ++j)
        {
          int kmj = k - j - 1;
          int ilo = kmj - 1;
          int jj;

          for (jj = 0; jj < kmj; ++jj)
            {
              double ajj = gsl_vector_get(&a.vector, jj);
              double ajjp1 = gsl_vector_get(&a.vector, jj + 1);
              double dl = gsl_vector_get(w->deltal, ilo);
              double dr = gsl_vector_get(w->deltar, jj);

              gsl_vector_set(&a.vector, jj, ((ajjp1 - ajj) / (dl + dr)) * kmj);
              ilo--;
            }
        }
    }

  if (nderiv != (size_t) km1)
    {
      for (j = (int) nderiv; j < km1; ++j)
        {
          int kmj = k - j - 1;
          int ilo = kmj - 1;
          int jj;

          for (jj = 0; jj < kmj; ++jj)
            {
              double ajj = gsl_vector_get(&a.vector, jj);
              double ajjp1 = gsl_vector_get(&a.vector, jj + 1);
              double dl = gsl_vector_get(w->deltal, ilo);
              double dr = gsl_vector_get(w->deltar, jj);

              gsl_vector_set(&a.vector, jj, (ajjp1 * dl + ajj * dr) / (dl + dr));
              ilo--;
            }
        }
    }

  *result = gsl_vector_get(&a.vector, 0);

  return GSL_SUCCESS;
}
