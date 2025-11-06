/* bspline/old.c
 *
 * Copyright (C) 2006, 2007, 2008, 2009 Patrick Alken
 * Copyright (C) 2008 Rhys Ulerich
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

/* this module contains old routines which will be deprecated in the future */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

/* Return number of coefficients */
size_t
gsl_bspline_ncoeffs (gsl_bspline_workspace * w)
{
  return gsl_bspline_ncontrol(w);
}

/*
gsl_bspline_knots()
  Compute the knots from the given breakpoints:

   knots(1:k) = breakpts(1)
   knots(k+1:k+l-1) = breakpts(i), i = 2 .. l
   knots(n+1:n+k) = breakpts(l + 1)

where l is the number of polynomial pieces (l = nbreak - 1) and
   n = k + l - 1
(using matlab syntax for the arrays)

The repeated knots at the beginning and end of the interval
correspond to the continuity condition there. See pg. 119
of [1].

Inputs: breakpts - breakpoints
        w        - bspline workspace

Return: success or error
*/

int
gsl_bspline_knots (const gsl_vector * breakpts, gsl_bspline_workspace * w)
{
  return gsl_bspline_init_augment(breakpts, w);
}

/*
gsl_bspline_knots_uniform()
  Construct uniformly spaced knots on the interval [a,b] using
the previously specified number of breakpoints. 'a' is the position
of the first breakpoint and 'b' is the position of the last
breakpoint.

Inputs: a - left side of interval
        b - right side of interval
        w - bspline workspace

Return: success or error

Notes: 1) w->knots is modified to contain the uniformly spaced
          knots

       2) The knots vector is set up as follows (using octave syntax):

          knots(1:k) = a
          knots(k+1:k+l-1) = a + i*delta, i = 1 .. l - 1
          knots(n+1:n+k) = b
*/

int
gsl_bspline_knots_uniform (const double a, const double b,
                           gsl_bspline_workspace * w)
{
  return gsl_bspline_init_uniform(a, b, w);
}

/*
gsl_bspline_eval()
  Evaluate the basis functions B_i(x) for all i. This is
a wrapper function for gsl_bspline_eval_nonzero() which
formats the output in a nice way.

Inputs: x - point for evaluation
        B - (output) where to store B_i(x) values
            the length of this vector is
            n = nbreak + k - 2 = l + k - 1 = w->n
        w - bspline workspace

Return: success or error

Notes: The w->knots vector must be initialized prior to calling
       this function (see gsl_bspline_knots())
*/

int
gsl_bspline_eval (const double x, gsl_vector * B, gsl_bspline_workspace * w)
{
  return gsl_bspline_eval_basis(x, B, w);
}

/*
gsl_bspline_eval_nonzero()
  Evaluate all non-zero B-spline functions at point x.
These are the B_i(x) for i in [istart, iend].
Always B_i(x) = 0 for i < istart and for i > iend.

Inputs: x      - point at which to evaluate splines
        Bk     - (output) where to store B-spline values (length k)
        istart - (output) B-spline function index of
                 first non-zero basis for given x
        iend   - (output) B-spline function index of
                 last non-zero basis for given x.
                 This is also the knot index corresponding to x.
        w      - bspline workspace

Return: success or error

Notes: 1) the w->knots vector must be initialized before calling
          this function

       2) On output, B contains

             [B_{istart,k}, B_{istart+1,k},
             ..., B_{iend-1,k}, B_{iend,k}]

          evaluated at the given x.
*/

int
gsl_bspline_eval_nonzero (const double x, gsl_vector * Bk, size_t * istart,
                          size_t * iend, gsl_bspline_workspace * w)
{
  int status;

  status = gsl_bspline_basis(x, Bk, istart, w);
  if (status)
    return status;

  *iend = *istart + w->spline_order - 1;

  return GSL_SUCCESS;
}

/*
gsl_bspline_deriv_eval()
  Evaluate d^j/dx^j B_i(x) for all i, 0 <= j <= nderiv.
This is a wrapper function for gsl_bspline_deriv_eval_nonzero()
which formats the output in a nice way.

Inputs: x      - point for evaluation
        nderiv - number of derivatives to compute, inclusive.
        dB     - (output) where to store d^j/dx^j B_i(x)
                 values. the size of this matrix is
                 (n = nbreak + k - 2 = l + k - 1 = w->n)
                 by (nderiv + 1)
        w      - bspline derivative workspace

Return: success or error

Notes: 1) The w->knots vector must be initialized prior to calling
          this function (see gsl_bspline_knots())

       2) based on PPPACK's bsplvd
*/

int
gsl_bspline_deriv_eval (const double x, const size_t nderiv,
                        gsl_matrix * dB, gsl_bspline_workspace * w)
{
  return gsl_bspline_eval_deriv_basis(x, nderiv, dB, w);
}

/*
gsl_bspline_deriv_eval_nonzero()
  At point x evaluate all requested, non-zero B-spline function
derivatives and store them in dB.  These are the
d^j/dx^j B_i(x) with i in [istart, iend] and j in [0, nderiv].
Always d^j/dx^j B_i(x) = 0 for i < istart and for i > iend.

Inputs: x      - point at which to evaluate splines
        nderiv - number of derivatives to request, inclusive
        dB     - (output) where to store dB-spline derivatives
                 (size k by nderiv + 1)
        istart - (output) B-spline function index of
                 first non-zero basis for given x
        iend   - (output) B-spline function index of
                 last non-zero basis for given x.
                 This is also the knot index corresponding to x.
        w      - bspline derivative workspace

Return: success or error

Notes: 1) the w->knots vector must be initialized before calling
          this function

       2) On output, dB contains

            [[B_{istart,  k}, ..., d^nderiv/dx^nderiv B_{istart  ,k}],
             [B_{istart+1,k}, ..., d^nderiv/dx^nderiv B_{istart+1,k}],
             ...
             [B_{iend-1,  k}, ..., d^nderiv/dx^nderiv B_{iend-1,  k}],
             [B_{iend,    k}, ..., d^nderiv/dx^nderiv B_{iend,    k}]]

          evaluated at x.  B_{istart, k} is stored in dB(0,0).
          Each additional column contains an additional derivative.

       3) Note that the zero-th column of the result contains the
          0th derivative, which is simply a function evaluation.

       4) based on PPPACK's bsplvd
*/

int
gsl_bspline_deriv_eval_nonzero (const double x, const size_t nderiv,
                                gsl_matrix * dB, size_t * istart,
                                size_t * iend, gsl_bspline_workspace * w)
{
  int status;

  status = gsl_bspline_basis_deriv(x, nderiv, dB, istart, w);
  if (status)
    return status;

  *iend = *istart + w->spline_order - 1;

  return GSL_SUCCESS;
}

int
gsl_bspline_knots_greville (const gsl_vector *abscissae,
                            gsl_bspline_workspace *w,
                            double *abserr)
{
  return gsl_bspline_init_greville(abscissae, w, abserr);
}
