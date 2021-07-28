/* interpolation/spline2d.c
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
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

gsl_spline2d *
gsl_spline2d_alloc(const gsl_interp2d_type * T, size_t xsize, size_t ysize)
{
  double * array_mem;
  gsl_spline2d * interp;

  if (xsize < T->min_size || ysize < T->min_size)
    {
      GSL_ERROR_NULL("insufficient number of points for interpolation type", GSL_EINVAL);
    }

  interp = calloc(1, sizeof(gsl_spline2d));
  if (interp == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for gsl_spline2d struct", GSL_ENOMEM);
    }

  interp->interp_object.type = T;
  interp->interp_object.xsize = xsize;
  interp->interp_object.ysize = ysize;
  if (interp->interp_object.type->alloc == NULL)
    {
      interp->interp_object.state = NULL;
    }
  else
    {
      interp->interp_object.state = interp->interp_object.type->alloc(xsize, ysize);
      if (interp->interp_object.state == NULL)
        {
          gsl_spline2d_free(interp);
          GSL_ERROR_NULL("failed to allocate space for gsl_spline2d state", GSL_ENOMEM);
        }
    }

  /*
   * Use one contiguous block of memory for all three data arrays.
   * That way the code fails immediately if there isn't sufficient space for everything,
   * rather than allocating one or two and then having to free them.
   */
  array_mem = (double *)calloc(xsize + ysize + xsize * ysize, sizeof(double));
  if (array_mem == NULL)
    {
      gsl_spline2d_free(interp);
      GSL_ERROR_NULL("failed to allocate space for data arrays", GSL_ENOMEM);
    }

  interp->xarr = array_mem;
  interp->yarr = array_mem + xsize;
  interp->zarr = array_mem + xsize + ysize;

  return interp;
} /* gsl_spline2d_alloc() */

int
gsl_spline2d_init(gsl_spline2d * interp, const double xarr[],
                  const double yarr[], const double zarr[],
                  size_t xsize, size_t ysize)
{
  int status = gsl_interp2d_init(&(interp->interp_object), xarr, yarr, zarr, xsize, ysize);

  memcpy(interp->xarr, xarr, xsize * sizeof(double));
  memcpy(interp->yarr, yarr, ysize * sizeof(double));
  memcpy(interp->zarr, zarr, xsize * ysize * sizeof(double));

  return status;
} /* gsl_spline2d_init() */

void
gsl_spline2d_free(gsl_spline2d * interp)
{
  RETURN_IF_NULL(interp);

  if (interp->interp_object.type->free)
    interp->interp_object.type->free(interp->interp_object.state);

  /*
   * interp->xarr points to the beginning of one contiguous block of memory
   * that holds interp->xarr, interp->yarr, and interp->zarr. So it all gets
   * freed with one call. cf. gsl_spline2d_alloc() implementation
   */
  if (interp->xarr)
    free(interp->xarr);

  free(interp);
} /* gsl_spline2d_free() */

double
gsl_spline2d_eval(const gsl_spline2d * interp, const double x,
                  const double y, gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  return gsl_interp2d_eval(&(interp->interp_object), interp->xarr, interp->yarr,
                           interp->zarr, x, y, xa, ya);
}

int
gsl_spline2d_eval_e(const gsl_spline2d * interp, const double x,
                    const double y, gsl_interp_accel * xa, gsl_interp_accel * ya,
                    double * z)
{
  return gsl_interp2d_eval_e(&(interp->interp_object), interp->xarr, interp->yarr,
                             interp->zarr, x, y, xa, ya, z);
}

double
gsl_spline2d_eval_extrap(const gsl_spline2d * interp, const double x,
                         const double y, gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  return gsl_interp2d_eval_extrap(&(interp->interp_object), interp->xarr, interp->yarr,
                                  interp->zarr, x, y, xa, ya);
}

int
gsl_spline2d_eval_extrap_e(const gsl_spline2d * interp, const double x,
                           const double y, gsl_interp_accel * xa, gsl_interp_accel * ya,
                           double * z)
{
  return gsl_interp2d_eval_extrap_e(&(interp->interp_object), interp->xarr, interp->yarr,
                                    interp->zarr, x, y, xa, ya, z);
}

double
gsl_spline2d_eval_deriv_x(const gsl_spline2d * interp, const double x,
                          const double y, gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  return gsl_interp2d_eval_deriv_x(&(interp->interp_object), interp->xarr, interp->yarr,
                                   interp->zarr, x, y, xa, ya);
}

int
gsl_spline2d_eval_deriv_x_e(const gsl_spline2d * interp, const double x,
                            const double y, gsl_interp_accel * xa, gsl_interp_accel * ya,
                            double * z)
{
  return gsl_interp2d_eval_deriv_x_e(&(interp->interp_object), interp->xarr, interp->yarr,
                                     interp->zarr, x, y, xa, ya, z);
}

double
gsl_spline2d_eval_deriv_y(const gsl_spline2d * interp, const double x,
                          const double y, gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  return gsl_interp2d_eval_deriv_y(&(interp->interp_object), interp->xarr, interp->yarr,
                                   interp->zarr, x, y, xa, ya);
}

int
gsl_spline2d_eval_deriv_y_e(const gsl_spline2d * interp, const double x,
                            const double y, gsl_interp_accel * xa, gsl_interp_accel * ya,
                            double * z)
{
  return gsl_interp2d_eval_deriv_y_e(&(interp->interp_object), interp->xarr, interp->yarr,
                                     interp->zarr, x, y, xa, ya, z);
}

double
gsl_spline2d_eval_deriv_xx(const gsl_spline2d * interp, const double x,
                           const double y, gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  return gsl_interp2d_eval_deriv_xx(&(interp->interp_object), interp->xarr, interp->yarr,
                                    interp->zarr, x, y, xa, ya);
}

int
gsl_spline2d_eval_deriv_xx_e(const gsl_spline2d * interp, const double x,
                             const double y, gsl_interp_accel * xa, gsl_interp_accel * ya,
                             double * z)
{
  return gsl_interp2d_eval_deriv_xx_e(&(interp->interp_object), interp->xarr, interp->yarr,
                                      interp->zarr, x, y, xa, ya, z);
}

double
gsl_spline2d_eval_deriv_yy(const gsl_spline2d * interp, const double x,
                           const double y, gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  return gsl_interp2d_eval_deriv_yy(&(interp->interp_object), interp->xarr, interp->yarr,
                                    interp->zarr, x, y, xa, ya);
}

int
gsl_spline2d_eval_deriv_yy_e(const gsl_spline2d * interp, const double x,
                             const double y, gsl_interp_accel * xa, gsl_interp_accel * ya,
                             double * z)
{
  return gsl_interp2d_eval_deriv_yy_e(&(interp->interp_object), interp->xarr, interp->yarr,
                                      interp->zarr, x, y, xa, ya, z);
}

double
gsl_spline2d_eval_deriv_xy(const gsl_spline2d * interp, const double x,
                           const double y, gsl_interp_accel * xa, gsl_interp_accel * ya)
{
  return gsl_interp2d_eval_deriv_xy(&(interp->interp_object), interp->xarr, interp->yarr,
                                    interp->zarr, x, y, xa, ya);
}

int
gsl_spline2d_eval_deriv_xy_e(const gsl_spline2d * interp, const double x,
                             const double y, gsl_interp_accel * xa, gsl_interp_accel * ya,
                             double * z)
{
  return gsl_interp2d_eval_deriv_xy_e(&(interp->interp_object), interp->xarr, interp->yarr,
                                      interp->zarr, x, y, xa, ya, z);
}

size_t
gsl_spline2d_min_size(const gsl_spline2d * interp)
{
  return gsl_interp2d_min_size(&(interp->interp_object));
}

const char *
gsl_spline2d_name(const gsl_spline2d * interp)
{
  return gsl_interp2d_name(&(interp->interp_object));
}

int
gsl_spline2d_set(const gsl_spline2d * interp, double zarr[],
                 const size_t i, const size_t j, const double z)
{
  return gsl_interp2d_set(&(interp->interp_object), zarr, i, j, z);
} /* gsl_spline2d_set() */

double
gsl_spline2d_get(const gsl_spline2d * interp, const double zarr[],
                 const size_t i, const size_t j)
{
  return gsl_interp2d_get(&(interp->interp_object), zarr, i, j);
} /* gsl_spline2d_get() */
