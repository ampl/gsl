/* interpolation/interp2d.c
 * 
 * Copyright 2012 David Zaslavsky
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>

/**
 * Triggers a GSL error if the argument is not equal to GSL_SUCCESS.
 * If the argument is GSL_SUCCESS, this does nothing.
 */
#define DISCARD_STATUS(s) if ((s) != GSL_SUCCESS) { GSL_ERROR_VAL("interpolation error", (s),  GSL_NAN); }

#define IDX2D(i, j, w) ((j) * ((w)->xsize) + (i))

gsl_interp2d *
gsl_interp2d_alloc(const gsl_interp2d_type * T, const size_t xsize,
                   const size_t ysize)
{
  gsl_interp2d * interp;

  if (xsize < T->min_size || ysize < T->min_size)
    {
      GSL_ERROR_NULL ("insufficient number of points for interpolation type",
                      GSL_EINVAL);
    }

  interp = (gsl_interp2d *) calloc(1, sizeof(gsl_interp2d));
  if (interp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for gsl_interp2d struct",
                      GSL_ENOMEM);
    }

  interp->type = T;
  interp->xsize = xsize;
  interp->ysize = ysize;

  if (interp->type->alloc == NULL)
    {
      interp->state = NULL;
      return interp;
    }

  interp->state = interp->type->alloc(xsize, ysize);
  if (interp->state == NULL)
    {
      free(interp);
      GSL_ERROR_NULL ("failed to allocate space for gsl_interp2d state",
                      GSL_ENOMEM);
    }

  return interp;
} /* gsl_interp2d_alloc() */

void
gsl_interp2d_free (gsl_interp2d * interp)
{
  RETURN_IF_NULL(interp);

  if (interp->type->free)
    interp->type->free(interp->state);

  free(interp);
} /* gsl_interp2d_free() */

int
gsl_interp2d_init (gsl_interp2d * interp, const double xarr[], const double yarr[],
                   const double zarr[], const size_t xsize, const size_t ysize)
{
  size_t i;

  if (xsize != interp->xsize || ysize != interp->ysize)
    {
      GSL_ERROR("data must match size of interpolation object", GSL_EINVAL);
    }

  for (i = 1; i < xsize; i++)
    {
      if (xarr[i-1] >= xarr[i])
        {
          GSL_ERROR("x values must be strictly increasing", GSL_EINVAL);
        }
    }

  for (i = 1; i < ysize; i++)
    {
      if (yarr[i-1] >= yarr[i])
        {
          GSL_ERROR("y values must be strictly increasing", GSL_EINVAL);
        }
    }

  interp->xmin = xarr[0];
  interp->xmax = xarr[xsize - 1];
  interp->ymin = yarr[0];
  interp->ymax = yarr[ysize - 1];

  {
    int status = interp->type->init(interp->state, xarr, yarr, zarr,
                                    xsize, ysize);
    return status;
  }
} /* gsl_interp2d_init() */

/*
 * A wrapper function that checks boundary conditions, calls an evaluator
 * which implements the actual calculation of the function value or 
 * derivative etc., and checks the return status.
 */
static int
interp2d_eval(int (*evaluator)(const void *, const double xa[], const double ya[],
                               const double za[], size_t xsize, size_t ysize,
                               double x, double y, gsl_interp_accel *,
                               gsl_interp_accel *, double * z),
                               const gsl_interp2d * interp, const double xarr[],
                               const double yarr[], const double zarr[],
                               const double x, const double y,
                               gsl_interp_accel * xa, gsl_interp_accel * ya,
                               double * result)
{
  if (x < interp->xmin || x > interp->xmax)
    {
      GSL_ERROR ("interpolation x value out of range", GSL_EDOM);
    }
  else if (y < interp->ymin || y > interp->ymax)
    {
      GSL_ERROR ("interpolation y value out of range", GSL_EDOM);
    }

  return evaluator(interp->state, xarr, yarr, zarr,
                   interp->xsize, interp->ysize,
                   x, y, xa, ya, result);
}

/*
 * Another wrapper function that serves as a drop-in replacement for
 * interp2d_eval but does not check the bounds. This can be used
 * for extrapolation.
 */
static int
interp2d_eval_extrap(int (*evaluator)(const void *, const double xa[],
                                      const double ya[], const double za[],
                                      size_t xsize, size_t ysize,
                                      double x, double y,
                                      gsl_interp_accel *,
                                      gsl_interp_accel *, double * z),
                                      const gsl_interp2d * interp, const double xarr[],
                                      const double yarr[], const double zarr[],
                                      const double x, const double y,
                                      gsl_interp_accel * xa, gsl_interp_accel * ya,
                                      double * result)
{
  return evaluator(interp->state, xarr, yarr, zarr,
                   interp->xsize, interp->ysize, x, y, xa, ya, result);
}

double
gsl_interp2d_eval (const gsl_interp2d * interp, const double xarr[],
                   const double yarr[], const double zarr[],
                   const double x, const double y,
                   gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  double z;
  int status = gsl_interp2d_eval_e(interp, xarr, yarr, zarr, x, y, xa, ya, &z);
  DISCARD_STATUS(status)
  return z;
} /* gsl_interp2d_eval() */

double
gsl_interp2d_eval_extrap (const gsl_interp2d * interp,
                          const double xarr[],
                          const double yarr[],
                          const double zarr[],
                          const double x,
                          const double y,
                          gsl_interp_accel * xa,
                          gsl_interp_accel * ya)
{
  double z;
  int status =
    interp2d_eval_extrap(interp->type->eval, interp,
                         xarr, yarr, zarr, x, y, xa, ya, &z);
  DISCARD_STATUS(status)
  return z;
}

int
gsl_interp2d_eval_e (const gsl_interp2d * interp, const double xarr[],
                     const double yarr[], const double zarr[],
                     const double x, const double y,
                     gsl_interp_accel * xa, gsl_interp_accel * ya, double * z)
{
  return interp2d_eval(interp->type->eval, interp,
                       xarr, yarr, zarr, x, y, xa, ya, z);
} /* gsl_interp2d_eval_e() */

#ifndef GSL_DISABLE_DEPRECATED

int
gsl_interp2d_eval_e_extrap (const gsl_interp2d * interp,
                            const double xarr[], const double yarr[],
                            const double zarr[], const double x,
                            const double y, gsl_interp_accel * xa,
                            gsl_interp_accel * ya, double * z)
{
  return interp2d_eval_extrap(interp->type->eval, interp,
                              xarr, yarr, zarr, x, y, xa, ya, z);
}

#endif /* !GSL_DISABLE_DEPRECATED */

int
gsl_interp2d_eval_extrap_e (const gsl_interp2d * interp,
                            const double xarr[], const double yarr[],
                            const double zarr[], const double x,
                            const double y, gsl_interp_accel * xa,
                            gsl_interp_accel * ya, double * z)
{
  return interp2d_eval_extrap(interp->type->eval, interp,
                              xarr, yarr, zarr, x, y, xa, ya, z);
}

double
gsl_interp2d_eval_deriv_x (const gsl_interp2d * interp, const double xarr[],
                           const double yarr[], const double zarr[],
                           const double x, const double y,
                           gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  double z;
  int status = gsl_interp2d_eval_deriv_x_e(interp, xarr, yarr, zarr, x, y, xa, ya, &z);
  DISCARD_STATUS(status)
  return z;
}

int
gsl_interp2d_eval_deriv_x_e (const gsl_interp2d * interp, const double xarr[],
                             const double yarr[], const double zarr[],
                             const double x, const double y,
                             gsl_interp_accel * xa, gsl_interp_accel * ya, double * z)
{
  return interp2d_eval(interp->type->eval_deriv_x, interp,
                       xarr, yarr, zarr, x, y, xa, ya, z);
}

double
gsl_interp2d_eval_deriv_y (const gsl_interp2d * interp, const double xarr[],
                           const double yarr[], const double zarr[],
                           const double x, const double y,
                           gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  double z;
  int status = gsl_interp2d_eval_deriv_y_e(interp, xarr, yarr, zarr, x, y, xa, ya, &z);
  DISCARD_STATUS(status)
  return z;
}

int
gsl_interp2d_eval_deriv_y_e (const gsl_interp2d * interp, const double xarr[],
                             const double yarr[], const double zarr[],
                             const double x, const double y,
                             gsl_interp_accel * xa, gsl_interp_accel * ya, double * z)
{
  return interp2d_eval(interp->type->eval_deriv_y, interp,
                       xarr, yarr, zarr, x, y, xa, ya, z);
}

double
gsl_interp2d_eval_deriv_xx (const gsl_interp2d * interp, const double xarr[],
                            const double yarr[], const double zarr[],
                            const double x, const double y,
                            gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  double z;
  int status = gsl_interp2d_eval_deriv_xx_e(interp, xarr, yarr, zarr, x, y, xa, ya, &z);
  DISCARD_STATUS(status)
  return z;
}

int
gsl_interp2d_eval_deriv_xx_e (const gsl_interp2d * interp, const double xarr[],
                              const double yarr[], const double zarr[],
                              const double x, const double y,
                              gsl_interp_accel * xa, gsl_interp_accel * ya, double * z)
{
  return interp2d_eval(interp->type->eval_deriv_xx, interp,
                       xarr, yarr, zarr, x, y, xa, ya, z);
}

double
gsl_interp2d_eval_deriv_yy (const gsl_interp2d * interp, const double xarr[],
                            const double yarr[], const double zarr[],
                            const double x, const double y,
                            gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  double z;
  int status = gsl_interp2d_eval_deriv_yy_e(interp, xarr, yarr, zarr, x, y, xa, ya, &z);
  DISCARD_STATUS(status)
  return z;
}

int
gsl_interp2d_eval_deriv_yy_e (const gsl_interp2d * interp, const double xarr[],
                              const double yarr[], const double zarr[],
                              const double x, const double y,
                              gsl_interp_accel * xa, gsl_interp_accel * ya, double * z)
{
  return interp2d_eval(interp->type->eval_deriv_yy, interp,
                       xarr, yarr, zarr, x, y, xa, ya, z);
}

double
gsl_interp2d_eval_deriv_xy (const gsl_interp2d * interp, const double xarr[],
                            const double yarr[], const double zarr[],
                            const double x, const double y,
                            gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  double z;
  int status = gsl_interp2d_eval_deriv_xy_e(interp, xarr, yarr, zarr, x, y, xa, ya, &z);
  DISCARD_STATUS(status)
  return z;
}

int
gsl_interp2d_eval_deriv_xy_e (const gsl_interp2d * interp, const double xarr[],
                              const double yarr[], const double zarr[],
                              const double x, const double y,
                              gsl_interp_accel * xa, gsl_interp_accel * ya, double * z)
{
  return interp2d_eval(interp->type->eval_deriv_xy, interp,
                       xarr, yarr, zarr, x, y, xa, ya, z);
}

size_t
gsl_interp2d_type_min_size(const gsl_interp2d_type * T)
{
  return T->min_size;
}

size_t
gsl_interp2d_min_size(const gsl_interp2d * interp)
{
  return interp->type->min_size;
}

const char *
gsl_interp2d_name(const gsl_interp2d * interp)
{
  return interp->type->name;
}

size_t
gsl_interp2d_idx(const gsl_interp2d * interp,
                 const size_t i, const size_t j)
{
  if (i >= interp->xsize)
    {
      GSL_ERROR_VAL ("x index out of range", GSL_ERANGE, 0);
    }
  else if (j >= interp->ysize)
    {
      GSL_ERROR_VAL ("y index out of range", GSL_ERANGE, 0);
    }
  else
    {
      return IDX2D(i, j, interp);
    }
} /* gsl_interp2d_idx() */

int
gsl_interp2d_set(const gsl_interp2d * interp, double zarr[],
                 const size_t i, const size_t j, const double z)
{
  if (i >= interp->xsize)
    {
      GSL_ERROR ("x index out of range", GSL_ERANGE);
    }
  else if (j >= interp->ysize)
    {
      GSL_ERROR ("y index out of range", GSL_ERANGE);
    }
  else
    {
      zarr[IDX2D(i, j, interp)] = z;
      return GSL_SUCCESS;
    }
} /* gsl_interp2d_set() */

double
gsl_interp2d_get(const gsl_interp2d * interp, const double zarr[],
                 const size_t i, const size_t j)
{
  if (i >= interp->xsize)
    {
      GSL_ERROR_VAL ("x index out of range", GSL_ERANGE, 0);
    }
  else if (j >= interp->ysize)
    {
      GSL_ERROR_VAL ("y index out of range", GSL_ERANGE, 0);
    }
  else
    {
      return zarr[IDX2D(i, j, interp)];
    }
} /* gsl_interp2d_get() */

#undef IDX2D
