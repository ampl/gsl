/* bspline/gsl_bspline.h
 *
 * Copyright (C) 2018, 2019, 2020, 2021 Patrick Alken
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

#ifndef __GSL_BSPLINE_H__
#define __GSL_BSPLINE_H__

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_inline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct
{
    size_t spline_order;  /* spline order */
    size_t nbreak;        /* number of breakpoints */
    size_t ncontrol;      /* number of bspline control points (nbreak + spline_order - 2) */

    gsl_vector *knots;    /* knots vector, length ncontrol + spline_order */
    gsl_vector *deltal;   /* left delta, length spline_order */
    gsl_vector *deltar;   /* right delta, length spline_order */
    gsl_vector *B;        /* temporary spline results, length spline_order */
    gsl_matrix *XTX;      /* stores diagonals of banded normal equations matrix, ncontrol-by-spline_order */
    gsl_matrix *R;        /* R factor for periodic least squares fitting, ncontrol-by-ncontrol */
    gsl_vector *work;     /* workspace, length 3*ncontrol */

    /* bspline derivative parameters */
    gsl_matrix *A;        /* work matrix, spline_order-by-spline_order */
    gsl_matrix *dB;       /* temporary derivative results, spline_order-by-[2*(spline_order+1)] */

    size_t icache;        /* cached index of current interval, in [0,n+k-2] */
} gsl_bspline_workspace;

/* bspline.c */

gsl_bspline_workspace * gsl_bspline_alloc(const size_t k, const size_t nbreak);

gsl_bspline_workspace * gsl_bspline_alloc_ncontrol (const size_t spline_order, const size_t ncontrol);

void gsl_bspline_free(gsl_bspline_workspace *w);

size_t gsl_bspline_ncontrol(const gsl_bspline_workspace * w);
size_t gsl_bspline_order(const gsl_bspline_workspace * w);
size_t gsl_bspline_nbreak(const gsl_bspline_workspace * w);
double gsl_bspline_breakpoint(const size_t i, const gsl_bspline_workspace * w);
double gsl_bspline_greville_abscissa(const size_t i, const gsl_bspline_workspace *w);

int gsl_bspline_init_augment(const gsl_vector *breakpts, gsl_bspline_workspace *w);

int gsl_bspline_init_uniform(const double a, const double b, gsl_bspline_workspace *w);

int gsl_bspline_init_periodic (const double a, const double b, gsl_bspline_workspace * w);

int gsl_bspline_init (const gsl_vector * t, gsl_bspline_workspace * w);

int gsl_bspline_proj_rhs(const gsl_function * F, gsl_vector * y, gsl_bspline_workspace * w);

INLINE_DECL size_t gsl_bspline_find_interval (const double x, int *flag, gsl_bspline_workspace * w);

/* eval.c */

int gsl_bspline_calc(const double x, const gsl_vector * c,
                     double * result, gsl_bspline_workspace * w);

int gsl_bspline_calc_deriv(const double x, const gsl_vector * c, const size_t nderiv,
                           double * result, gsl_bspline_workspace * w);

int gsl_bspline_vector_calc(const double x, const gsl_matrix * c, gsl_vector * result,
                            gsl_bspline_workspace * w);

int gsl_bspline_vector_calc_deriv(const double x, const gsl_matrix * c, const size_t nderiv,
                                  gsl_vector * result, gsl_bspline_workspace * w);

int
gsl_bspline_eval_basis(const double x, gsl_vector *B, gsl_bspline_workspace *w);

int
gsl_bspline_basis(const double x, gsl_vector *Bk, size_t *istart,
                  gsl_bspline_workspace *w);

int
gsl_bspline_eval_deriv_basis(const double x, const size_t nderiv,
                             gsl_matrix *dB, gsl_bspline_workspace *w);

int
gsl_bspline_basis_deriv(const double x, const size_t nderiv,
                        gsl_matrix *dB, size_t *istart,
                        gsl_bspline_workspace *w);

int
gsl_bspline_init_greville(const gsl_vector *abscissae,
                          gsl_bspline_workspace *w,
                          double *abserr);

/* gram.c */

int gsl_bspline_gram(const size_t nderiv, gsl_matrix * G, gsl_bspline_workspace * w);

int gsl_bspline_gram_interval(const double a, const double b, const size_t nderiv,
                              gsl_matrix * G, gsl_bspline_workspace * w);

int gsl_bspline_oprod(const size_t nderiv, const double x, gsl_matrix * A, gsl_bspline_workspace * w);

/* integ.c */

int gsl_bspline_calc_integ(const double a, const double b,
                           const gsl_vector * c, double * result,
                           gsl_bspline_workspace * w);

int gsl_bspline_basis_integ(const double a, const double b,
                            gsl_vector * bint, gsl_bspline_workspace * w);

/* interp.c */

int gsl_bspline_init_interp (const gsl_vector * x, gsl_bspline_workspace * w);

int gsl_bspline_init_hermite(const size_t nderiv, const gsl_vector * x, gsl_bspline_workspace * w);

int gsl_bspline_col_interp(const gsl_vector * tau, gsl_matrix * XB, gsl_bspline_workspace * w);

int gsl_bspline_interp_chermite(const gsl_vector * x, const gsl_vector * y,
                                const gsl_vector * dy, gsl_vector * c,
                                const gsl_bspline_workspace * w);

/* ls.c */

int gsl_bspline_lssolve(const gsl_vector * x, const gsl_vector * y, gsl_vector * c,
                        double * chisq, gsl_bspline_workspace * w);

int gsl_bspline_wlssolve(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts,
                         gsl_vector * c, double * chisq, gsl_bspline_workspace * w);

int gsl_bspline_lsnormal(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts,
                         gsl_vector * XTy, gsl_matrix * XTX, gsl_bspline_workspace * w);

int gsl_bspline_lsnormalm(const gsl_vector * x, const gsl_matrix * Y, const gsl_vector * wts,
                          gsl_matrix * XTY, gsl_matrix * XTX, gsl_bspline_workspace * w);

int gsl_bspline_plssolve(const gsl_vector * x, const gsl_vector * y,
                         gsl_vector * c, double * chisq, gsl_bspline_workspace * w);

int gsl_bspline_pwlssolve(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts,
                          gsl_vector * c, double * chisq, gsl_bspline_workspace * w);

int gsl_bspline_plsqr(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts,
                      gsl_matrix * R, gsl_vector * QTy, double * rnorm,
                      gsl_bspline_workspace * w);

int gsl_bspline_residuals(const gsl_vector * x, const gsl_vector * y, const gsl_vector * c,
                          gsl_vector * r, gsl_bspline_workspace * w);

int gsl_bspline_covariance(const gsl_matrix * XTX, gsl_matrix * cov, gsl_bspline_workspace * w);

int gsl_bspline_rcond(const gsl_matrix * XTX, double * rcond, gsl_bspline_workspace * w);

int gsl_bspline_err(const double x, const size_t nderiv,
                    const gsl_matrix * cov, double * err,
                    gsl_bspline_workspace * w);

/* future to be deprecated functions */

size_t gsl_bspline_ncoeffs (gsl_bspline_workspace * w);
int gsl_bspline_knots (const gsl_vector * breakpts, gsl_bspline_workspace * w);
int gsl_bspline_knots_uniform (const double a, const double b, gsl_bspline_workspace * w);
int gsl_bspline_eval (const double x, gsl_vector * B, gsl_bspline_workspace * w);
int gsl_bspline_deriv_eval (const double x, const size_t nderiv,
                            gsl_matrix * dB, gsl_bspline_workspace * w);
int gsl_bspline_knots_greville (const gsl_vector *abscissae,
                                gsl_bspline_workspace *w,
                                double *abserr);

/* end of future to be deprecated functions */

#ifdef HAVE_INLINE

/*
gsl_bspline_find_interval()
  Find knot interval such that t_i <= x < t_{i + 1}
where the t_i are knot values. The algorithm uses binary search
since the knot array is sorted ascending.

Inputs: x    - x value
        flag - (output) error flag
        w    - bspline workspace

Return: i (index in w->knots corresponding to left limit of interval)

Notes: The error conditions are reported as follows:

       Condition                        Return value        Flag
       ---------                        ------------        ----
       x < t_0                               0               -1
       t_i <= x < t_{i+1}                    i                0
       t_i < x = t_{i+1} = t_{n+k-1}         i                0
       t_i < t_{i+1} = t_{n+k-1} < x         i               +1
*/

INLINE_FUN
size_t
gsl_bspline_find_interval (const double x, int *flag, gsl_bspline_workspace * w)
{
  const double *array = w->knots->data; /* stride is 1 */
  size_t ilo;
  size_t ihi = w->knots->size - 1;

  /* check for quick return */
  if (array[w->icache] <= x && x < array[w->icache + 1])
    {
      *flag = 0;
      return w->icache;
    }
  else if (x < *array)
    {
      *flag = -1;
      return 0;
    }
  else if (x >= array[ihi])
    {
      *flag = (x > array[ihi]);

      while (ihi > 0)
        {
          ihi--;
          if (array[ihi] < array[ihi + 1])
            break;
        }

      return ihi;
    }

  *flag = 0;
  ilo = 0;

  while (ihi > ilo + 1)
    {
      size_t i = (ihi + ilo) >> 1;
      if(array[i] > x)
        ihi = i;
      else
        ilo = i;
    }

  w->icache = ilo;

  return ilo;
}

#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_BSPLINE_H__ */
