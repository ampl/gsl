/* ode-initval/odeiv.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

gsl_odeiv2_step *
gsl_odeiv2_step_alloc (const gsl_odeiv2_step_type * T, size_t dim)
{
  gsl_odeiv2_step *s = (gsl_odeiv2_step *) malloc (sizeof (gsl_odeiv2_step));

  if (s == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for ode struct", GSL_ENOMEM);
    };

  s->type = T;
  s->dimension = dim;

  s->state = s->type->alloc (dim);

  if (s->state == 0)
    {
      free (s);                 /* exception in constructor, avoid memory leak */

      GSL_ERROR_NULL ("failed to allocate space for ode state", GSL_ENOMEM);
    };

  return s;
}

const char *
gsl_odeiv2_step_name (const gsl_odeiv2_step * s)
{
  return s->type->name;
}

unsigned int
gsl_odeiv2_step_order (const gsl_odeiv2_step * s)
{
  return s->type->order (s->state);
}

int
gsl_odeiv2_step_apply (gsl_odeiv2_step * s,
                       double t,
                       double h,
                       double y[],
                       double yerr[],
                       const double dydt_in[],
                       double dydt_out[], const gsl_odeiv2_system * dydt)
{
  return s->type->apply (s->state, s->dimension, t, h, y, yerr, dydt_in,
                         dydt_out, dydt);
}

int
gsl_odeiv2_step_reset (gsl_odeiv2_step * s)
{
  return s->type->reset (s->state, s->dimension);
}

void
gsl_odeiv2_step_free (gsl_odeiv2_step * s)
{
  RETURN_IF_NULL (s);
  s->type->free (s->state);
  free (s);
}

int
gsl_odeiv2_step_set_driver (gsl_odeiv2_step * s, const gsl_odeiv2_driver * d)
{
  if (d != NULL)
    {
      s->type->set_driver (s->state, d);
    }
  else
    {
      GSL_ERROR ("driver pointer is null", GSL_EFAULT);
    }

  return GSL_SUCCESS;
}
