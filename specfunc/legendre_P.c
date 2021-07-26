/* specfunc/legendre_P.c
 * 
 * Copyright (C) 2009-2013 Patrick Alken
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

/*
 * The routines in this module compute associated Legendre functions
 * (ALFs) up to order and degree 2700, using the method described
 * in
 *
 * [1] S. A. Holmes and W. E. Featherstone, A unified approach
 *     to the Clenshaw summation and the recursive computation of very
 *     high degree and order normalised associated Legendre functions,
 *     Journal of Geodesy, 76, pg. 279-299, 2002.
 *
 * Further information on ALFs can be found in
 *
 * [2] Abramowitz and Stegun, Handbook of Mathematical Functions,
 *     Chapter 8, 1972.
 */

static void legendre_sqrts(const size_t lmax, double *array);

#define LEGENDRE
#include "legendre_source.c"
#undef LEGENDRE

#define LEGENDRE_DERIV
#include "legendre_source.c"
#undef LEGENDRE_DERIV

#define LEGENDRE_DERIV_ALT
#include "legendre_source.c"
#undef LEGENDRE_DERIV_ALT

#define LEGENDRE_DERIV2
#include "legendre_source.c"
#undef LEGENDRE_DERIV2

#define LEGENDRE_DERIV2_ALT
#include "legendre_source.c"
#undef LEGENDRE_DERIV2_ALT

/* number of P_{lm} functions for a given lmax */
size_t
gsl_sf_legendre_nlm(const size_t lmax)
{
  return ((lmax + 1) * (lmax + 2) / 2);
}

/*
gsl_sf_legendre_array_n()
  This routine returns the minimum result_array[] size needed
for a given lmax
*/

size_t
gsl_sf_legendre_array_n(const size_t lmax)
{
  size_t nlm = gsl_sf_legendre_nlm(lmax);
  size_t nsqrt = 2 * lmax + 2; /* extra room to precompute sqrt factors */

  return (nlm + nsqrt);
}

/*********************************************************
 *                 INTERNAL ROUTINES                     *
 *********************************************************/

/*
legendre_sqrts()
  Precompute square root factors needed for Legendre recurrence.
On output, array[i] = sqrt(i)
*/

static void
legendre_sqrts(const size_t lmax, double *array)
{
  size_t l;
  for (l = 0; l <= 2 * lmax + 1; ++l)
    array[l] = sqrt((double) l);
}
