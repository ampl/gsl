/* specfunc/gsl_sf_alf.h
 * 
 * Copyright (C) 2023 Patrick Alken
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

#ifndef __GSL_SF_ALF_H__
#define __GSL_SF_ALF_H__

#include <stdlib.h>
#include <gsl/gsl_inline.h>
#include <gsl/gsl_sf_result.h>

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

typedef enum
{
  GSL_SF_ALF_SCHMIDT,
  GSL_SF_ALF_SPHARM,
  GSL_SF_ALF_FULL,
  GSL_SF_ALF_FOURPI,
  GSL_SF_ALF_NONE
} gsl_sf_alf_t;

#define GSL_SF_ALF_FLG_CSPHASE       (1 << 0) /* include Condon-Shortley phase */

int gsl_sf_alf_precompute(const gsl_sf_alf_t norm, const size_t lmax,
                          const size_t mmax, const size_t flags, double output_array[]);
size_t gsl_sf_alf_array_size(const size_t lmax, const size_t mmax);
int gsl_sf_alf_array(const size_t lmax, const size_t mmax, const double x, double result_array[]);
int gsl_sf_alf_deriv_array(const size_t lmax, const size_t mmax, const double x,
                           double result_array[], double result_deriv_array[]);
int gsl_sf_alf_vsh_array(const size_t lmax, const size_t mmax, const double x,
                         double result_array[], double result_deriv_array[]);

INLINE_DECL size_t gsl_sf_alf_nlm(const size_t lmax, const size_t mmax);
INLINE_DECL size_t gsl_sf_alf_array_index(const size_t l, const size_t m, const size_t lmax);

#ifdef HAVE_INLINE

/*
gsl_sf_alf_array_index()
This routine computes the index into a result_array[] corresponding
to a given (l,m) using 'M' major indexing:

(l,m) = (0,0) (1,0) (2,0) ... (L,0) (1,1) (2,1) ... (L,1) ... (L,L)

index(l,m,L) = m*L - m(m-1)/2 + l
*/

INLINE_FUN
size_t
gsl_sf_alf_array_index(const size_t l, const size_t m, const size_t lmax)
{
  return (((m * (2 * lmax + 1 - m)) >> 1) + l);
}

/* return number of ALFs for a given lmax = (mmax+1)*(lmax+1) - mmax*(mmax+1)/2 */
INLINE_FUN
size_t
gsl_sf_alf_nlm(const size_t lmax, const size_t mmax)
{
  const size_t M = (mmax > lmax) ? lmax : mmax;
  return ((lmax + 1) * (M + 1) - ((M * (M + 1)) >> 1));
}

#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_SF_ALF_H__ */
