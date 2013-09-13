/* robust_wfun.c
 * 
 * Copyright (C) 2013 Patrick Alken
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

/* default tuning constants */
#define TUNING_BISQUARE       (4.685)
#define TUNING_CAUCHY         (2.385)
#define TUNING_FAIR           (1.4)
#define TUNING_HUBER          (1.345)
#define TUNING_OLS            (1.0)
#define TUNING_WELSCH         (2.985)

static int
bisquare(const gsl_vector *r, gsl_vector *w)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);

      if (fabs(ri) < 1.0)
        gsl_vector_set(w, i, (1.0 - ri*ri)*(1.0 - ri*ri));
      else
        gsl_vector_set(w, i, 0.0);
    }

  return GSL_SUCCESS;
} /* bisquare() */

static int
bisquare_dpsi(const gsl_vector *r, gsl_vector *dpsi)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);

      if (fabs(ri) < 1.0)
        gsl_vector_set(dpsi, i, (1.0 - ri*ri)*(1.0 - 5.0*ri*ri));
      else
        gsl_vector_set(dpsi, i, 0.0);
    }

  return GSL_SUCCESS;
} /* bisquare_dpsi() */

static const gsl_multifit_robust_type bisquare_type = {
  "bisquare",
  &bisquare,
  &bisquare_dpsi,
  TUNING_BISQUARE
};

static int
cauchy(const gsl_vector *r, gsl_vector *w)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);

      gsl_vector_set(w, i, 1.0 / (1.0 + ri*ri));
    }

  return GSL_SUCCESS;
} /* cauchy() */

static int
cauchy_dpsi(const gsl_vector *r, gsl_vector *dpsi)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);
      double rsq = ri * ri;

      gsl_vector_set(dpsi, i, (1 - rsq) / (1.0 + rsq) / (1.0 + rsq));
    }

  return GSL_SUCCESS;
} /* cauchy_dpsi() */

static const gsl_multifit_robust_type cauchy_type = {
  "cauchy",
  &cauchy,
  &cauchy_dpsi,
  TUNING_CAUCHY
};

static int
fair(const gsl_vector *r, gsl_vector *w)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);

      gsl_vector_set(w, i, 1.0 / (1.0 + fabs(ri)));
    }

  return GSL_SUCCESS;
} /* fair() */

static int
fair_dpsi(const gsl_vector *r, gsl_vector *dpsi)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);

      gsl_vector_set(dpsi, i, 1.0 / (1.0 + fabs(ri)) / (1.0 + fabs(ri)));
    }

  return GSL_SUCCESS;
} /* fair_dpsi() */

static const gsl_multifit_robust_type fair_type = {
  "fair",
  &fair,
  &fair_dpsi,
  TUNING_FAIR
};

static int
huber(const gsl_vector *r, gsl_vector *w)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double absri = fabs(gsl_vector_get(r, i));

      if (absri <= 1.0)
        gsl_vector_set(w, i, 1.0);
      else
        gsl_vector_set(w, i, 1.0 / absri);
    }

  return GSL_SUCCESS;
} /* huber() */

static int
huber_dpsi(const gsl_vector *r, gsl_vector *dpsi)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);

      if (fabs(ri) <= 1.0)
        gsl_vector_set(dpsi, i, 1.0);
      else
        gsl_vector_set(dpsi, i, 0.0);
    }

  return GSL_SUCCESS;
} /* huber_dpsi() */

static const gsl_multifit_robust_type huber_type = {
  "huber",
  &huber,
  &huber_dpsi,
  TUNING_HUBER
};

static int
ols(const gsl_vector *r, gsl_vector *w)
{
  gsl_vector_set_all(w, 1.0);

  return GSL_SUCCESS;
}

static int
ols_dpsi(const gsl_vector *r, gsl_vector *dpsi)
{
  gsl_vector_set_all(dpsi, 1.0);

  return GSL_SUCCESS;
}

static const gsl_multifit_robust_type ols_type = {
  "ols",
  &ols,
  &ols_dpsi,
  TUNING_OLS
};

static int
welsch(const gsl_vector *r, gsl_vector *w)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);

      gsl_vector_set(w, i, exp(-ri*ri));
    }

  return GSL_SUCCESS;
} /* welsch() */

static int
welsch_dpsi(const gsl_vector *r, gsl_vector *dpsi)
{
  size_t i;
  size_t n = r->size;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);

      gsl_vector_set(dpsi, i, (1.0 - 2.0*ri*ri) * exp(-ri*ri));
    }

  return GSL_SUCCESS;
} /* welsch_dpsi() */

static const gsl_multifit_robust_type welsch_type = {
  "welsch",
  &welsch,
  &welsch_dpsi,
  TUNING_WELSCH
};

const gsl_multifit_robust_type *gsl_multifit_robust_default = &bisquare_type;
const gsl_multifit_robust_type *gsl_multifit_robust_bisquare = &bisquare_type;
const gsl_multifit_robust_type *gsl_multifit_robust_cauchy = &cauchy_type;
const gsl_multifit_robust_type *gsl_multifit_robust_fair = &fair_type;
const gsl_multifit_robust_type *gsl_multifit_robust_huber = &huber_type;
const gsl_multifit_robust_type *gsl_multifit_robust_ols = &ols_type;
const gsl_multifit_robust_type *gsl_multifit_robust_welsch = &welsch_type;
