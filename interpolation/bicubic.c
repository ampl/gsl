/* interpolation/bicubic.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>

#define IDX2D(i, j, w) ((j) * ((w)->xsize) + (i))

typedef struct
{
  double * zx;
  double * zy;
  double * zxy;
  size_t xsize;
  size_t ysize;
} bicubic_state_t;

static void bicubic_free (void * vstate);

static void *
bicubic_alloc(size_t xsize, size_t ysize)
{
  bicubic_state_t *state;
  
  state = calloc(1, sizeof (bicubic_state_t));

  if (state == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for state", GSL_ENOMEM);
    }

  state->zx = (double *) malloc (xsize * ysize * sizeof (double));
  if (state->zx == NULL)
    {
      bicubic_free(state);
      GSL_ERROR_NULL("failed to allocate space for zx", GSL_ENOMEM);
    }

  state->zy = (double *) malloc (xsize * ysize * sizeof (double));
  if (state->zy == NULL)
    {
      bicubic_free(state);
      GSL_ERROR_NULL("failed to allocate space for zy", GSL_ENOMEM);
    }

  state->zxy = (double *) malloc (xsize * ysize * sizeof (double));
  if (state->zxy == NULL)
    {
      bicubic_free(state);
      GSL_ERROR_NULL("failed to allocate space for zxy", GSL_ENOMEM);
    }

  state->xsize = xsize;
  state->ysize = ysize;

  return state;
} /* bicubic_alloc() */

static void
bicubic_free (void * vstate)
{
  bicubic_state_t *state = (bicubic_state_t *) vstate;

  RETURN_IF_NULL(state);

  if (state->zx)
    free (state->zx);

  if (state->zy)
    free (state->zy);

  if (state->zxy)
    free (state->zxy);

  free (state);
} /* bicubic_free() */

static int
bicubic_init(void * vstate, const double xa[], const double ya[],
             const double za[], size_t xsize, size_t ysize)
{
  bicubic_state_t *state = (bicubic_state_t *) vstate;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline;
  gsl_vector *x;
  gsl_vector *y;
  size_t i, j;

  x = gsl_vector_alloc(xsize);
  y = gsl_vector_alloc(xsize);
  spline = gsl_spline_alloc(gsl_interp_cspline, xsize);
  for (j = 0; j <= ysize - 1; j++)
    {
      for (i = 0; i <= xsize - 1; i++)
        {
          gsl_vector_set(x, i, xa[i]);
          gsl_vector_set(y, i, za[IDX2D(i, j, state)]);
        }
      gsl_spline_init(spline, x->data, y->data, xsize);
      for (i = 0; i <= xsize - 1; i++)
        {
          state->zx[IDX2D(i, j, state)] = gsl_spline_eval_deriv(spline, xa[i], acc);
        }
    }
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_spline_free(spline);
  gsl_interp_accel_reset(acc);

  x = gsl_vector_alloc(ysize);
  y = gsl_vector_alloc(ysize);
  spline = gsl_spline_alloc(gsl_interp_cspline, ysize);
  for (i = 0; i <= xsize - 1; i++)
    {
      for (j = 0; j <= ysize - 1; j++)
        {
          gsl_vector_set(x, j, ya[j]);
          gsl_vector_set(y, j, za[IDX2D(i, j, state)]);
        }
      gsl_spline_init(spline, x->data, y->data, ysize);
      for (j = 0; j <= ysize - 1; j++)
        {
          state->zy[IDX2D(i, j, state)] = gsl_spline_eval_deriv(spline, ya[j], acc);
        }
    }
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_spline_free(spline);
  gsl_interp_accel_reset(acc);

  x = gsl_vector_alloc(xsize);
  y = gsl_vector_alloc(xsize);
  spline = gsl_spline_alloc(gsl_interp_cspline, xsize);
  for (j = 0; j <= ysize - 1; j++)
    {
      for (i = 0; i <= xsize - 1; i++)
        {
          gsl_vector_set(x, i, xa[i]);
          gsl_vector_set(y, i, state->zy[IDX2D(i, j, state)]);
        }
      gsl_spline_init(spline, x->data, y->data, xsize);
      for (i = 0; i <= xsize - 1; i++)
        {
          state->zxy[IDX2D(i, j, state)] = gsl_spline_eval_deriv(spline, xa[i], acc);
        }
    }
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return GSL_SUCCESS;
} /* bicubic_init() */


static int
bicubic_eval(const void * vstate, const double xarr[], const double yarr[],
             const double zarr[], size_t xsize, size_t ysize,
             double x, double y, gsl_interp_accel * xa,
             gsl_interp_accel * ya, double * z)
{
  bicubic_state_t *state = (bicubic_state_t *) vstate;

  double xmin, xmax, ymin, ymax;
  double zminmin, zminmax, zmaxmin, zmaxmax;
  double zxminmin, zxminmax, zxmaxmin, zxmaxmax;
  double zyminmin, zyminmax, zymaxmin, zymaxmax;
  double zxyminmin, zxyminmax, zxymaxmin, zxymaxmax;

  double dx, dy;   /* size of the grid cell */
  double dt, du;

  /*
   * t and u are the positions within the grid cell at which we are computing
   * the interpolation, in units of grid cell size
   */
  double t, u;
  double t0, t1, t2, t3, u0, u1, u2, u3;
  double v;
  size_t xi, yi;

  /* first compute the indices into the data arrays where we are interpolating */
  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  /* find the minimum and maximum values on the grid cell in each dimension */
  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, state)];
  zminmax = zarr[IDX2D(xi, yi + 1, state)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, state)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, state)];
  /* Get the width and height of the grid cell */
  dx = xmax - xmin;
  dy = ymax - ymin;
  t = (x - xmin)/dx;
  u = (y - ymin)/dy;
  dt = 1./dx; /* partial t / partial x */
  du = 1./dy; /* partial u / partial y */
  zxminmin = state->zx[IDX2D(xi, yi, state)]/dt;
  zxminmax = state->zx[IDX2D(xi, yi + 1, state)]/dt;
  zxmaxmin = state->zx[IDX2D(xi + 1, yi, state)]/dt;
  zxmaxmax = state->zx[IDX2D(xi + 1, yi + 1, state)]/dt;
  zyminmin = state->zy[IDX2D(xi, yi, state)]/du;
  zyminmax = state->zy[IDX2D(xi, yi + 1, state)]/du;
  zymaxmin = state->zy[IDX2D(xi + 1, yi, state)]/du;
  zymaxmax = state->zy[IDX2D(xi + 1, yi + 1, state)]/du;
  zxyminmin = state->zxy[IDX2D(xi, yi, state)]/(dt*du);
  zxyminmax = state->zxy[IDX2D(xi, yi + 1, state)]/(dt*du);
  zxymaxmin = state->zxy[IDX2D(xi + 1, yi, state)]/(dt*du);
  zxymaxmax = state->zxy[IDX2D(xi + 1, yi + 1, state)]/(dt*du);
  t0 = 1;
  t1 = t;
  t2 = t*t;
  t3 = t*t2;
  u0 = 1;
  u1 = u;
  u2 = u*u;
  u3 = u*u2;

  *z = 0;
  v = zminmin;
  *z += v*t0*u0;
  v = zyminmin;
  *z += v*t0*u1;
  v = -3*zminmin + 3*zminmax - 2*zyminmin - zyminmax;
  *z += v*t0*u2;
  v = 2*zminmin - 2*zminmax + zyminmin + zyminmax;
  *z += v*t0*u3;
  v = zxminmin;
  *z += v*t1*u0;
  v = zxyminmin;
  *z += v*t1*u1;
  v = -3*zxminmin + 3*zxminmax - 2*zxyminmin - zxyminmax;
  *z += v*t1*u2;
  v = 2*zxminmin - 2*zxminmax + zxyminmin + zxyminmax;
  *z += v*t1*u3;
  v = -3*zminmin + 3*zmaxmin - 2*zxminmin - zxmaxmin;
  *z += v*t2*u0;
  v = -3*zyminmin + 3*zymaxmin - 2*zxyminmin - zxymaxmin;
  *z += v*t2*u1;
  v = 9*zminmin - 9*zmaxmin + 9*zmaxmax - 9*zminmax + 6*zxminmin + 3*zxmaxmin - 3*zxmaxmax - 6*zxminmax + 6*zyminmin - 6*zymaxmin - 3*zymaxmax + 3*zyminmax + 4*zxyminmin + 2*zxymaxmin + zxymaxmax + 2*zxyminmax;
  *z += v*t2*u2;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 4*zxminmin - 2*zxmaxmin + 2*zxmaxmax + 4*zxminmax - 3*zyminmin + 3*zymaxmin + 3*zymaxmax - 3*zyminmax - 2*zxyminmin - zxymaxmin - zxymaxmax - 2*zxyminmax;
  *z += v*t2*u3;
  v = 2*zminmin - 2*zmaxmin + zxminmin + zxmaxmin;
  *z += v*t3*u0;
  v = 2*zyminmin - 2*zymaxmin + zxyminmin + zxymaxmin;
  *z += v*t3*u1;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 3*zxminmin - 3*zxmaxmin + 3*zxmaxmax + 3*zxminmax - 4*zyminmin + 4*zymaxmin + 2*zymaxmax - 2*zyminmax - 2*zxyminmin - 2*zxymaxmin - zxymaxmax - zxyminmax;
  *z += v*t3*u2;
  v = 4*zminmin - 4*zmaxmin + 4*zmaxmax - 4*zminmax + 2*zxminmin + 2*zxmaxmin - 2*zxmaxmax - 2*zxminmax + 2*zyminmin - 2*zymaxmin - 2*zymaxmax + 2*zyminmax + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
  *z += v*t3*u3;

  return GSL_SUCCESS;
} /* bicubic_eval() */

static int
bicubic_deriv_x(const void * vstate, const double xarr[], const double yarr[],
                const double zarr[], size_t xsize, size_t ysize,
                double x, double y,
                gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_p)
{
  bicubic_state_t *state = (bicubic_state_t *) vstate;

  double xmin, xmax, ymin, ymax;
  double zminmin, zminmax, zmaxmin, zmaxmax;
  double zxminmin, zxminmax, zxmaxmin, zxmaxmax;
  double zyminmin, zyminmax, zymaxmin, zymaxmax;
  double zxyminmin, zxyminmax, zxymaxmin, zxymaxmax;
  double dx, dy; /* size of the grid cell */
  double dt, du;

  /*
   * t and u are the positions within the grid cell at which we are computing
   * the interpolation, in units of grid cell size
   */
  double t, u;
  double t0, t1, t2, u0, u1, u2, u3;
  double v;
  size_t xi, yi;

  /* first compute the indices into the data arrays where we are interpolating */
  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  /* find the minimum and maximum values on the grid cell in each dimension */
  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, state)];
  zminmax = zarr[IDX2D(xi, yi + 1, state)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, state)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, state)];

  /* get the width and height of the grid cell */
  dx = xmax - xmin;
  dy = ymax - ymin;
  t = (x - xmin)/dx;
  u = (y - ymin)/dy;
  dt = 1./dx; /* partial t / partial x */
  du = 1./dy; /* partial u / partial y */

  zxminmin = state->zx[IDX2D(xi, yi, state)]/dt;
  zxminmax = state->zx[IDX2D(xi, yi + 1, state)]/dt;
  zxmaxmin = state->zx[IDX2D(xi + 1, yi, state)]/dt;
  zxmaxmax = state->zx[IDX2D(xi + 1, yi + 1, state)]/dt;
  zyminmin = state->zy[IDX2D(xi, yi, state)]/du;
  zyminmax = state->zy[IDX2D(xi, yi + 1, state)]/du;
  zymaxmin = state->zy[IDX2D(xi + 1, yi, state)]/du;
  zymaxmax = state->zy[IDX2D(xi + 1, yi + 1, state)]/du;
  zxyminmin = state->zxy[IDX2D(xi, yi, state)]/(dt*du);
  zxyminmax = state->zxy[IDX2D(xi, yi + 1, state)]/(dt*du);
  zxymaxmin = state->zxy[IDX2D(xi + 1, yi, state)]/(dt*du);
  zxymaxmax = state->zxy[IDX2D(xi + 1, yi + 1, state)]/(dt*du);

  t0 = 1;
  t1 = t;
  t2 = t*t;
  u0 = 1;
  u1 = u;
  u2 = u*u;
  u3 = u*u2;

  *z_p = 0;
  v = zxminmin;
  *z_p += v*t0*u0;
  v = zxyminmin;
  *z_p += v*t0*u1;
  v = -3*zxminmin + 3*zxminmax - 2*zxyminmin - zxyminmax;
  *z_p += v*t0*u2;
  v = 2*zxminmin - 2*zxminmax + zxyminmin + zxyminmax;
  *z_p += v*t0*u3;
  v = -3*zminmin + 3*zmaxmin - 2*zxminmin - zxmaxmin;
  *z_p += 2*v*t1*u0;
  v = -3*zyminmin + 3*zymaxmin - 2*zxyminmin - zxymaxmin;
  *z_p += 2*v*t1*u1;
  v = 9*zminmin - 9*zmaxmin + 9*zmaxmax - 9*zminmax + 6*zxminmin + 3*zxmaxmin - 3*zxmaxmax - 6*zxminmax + 6*zyminmin - 6*zymaxmin - 3*zymaxmax + 3*zyminmax + 4*zxyminmin + 2*zxymaxmin + zxymaxmax + 2*zxyminmax;
  *z_p += 2*v*t1*u2;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 4*zxminmin - 2*zxmaxmin + 2*zxmaxmax + 4*zxminmax - 3*zyminmin + 3*zymaxmin + 3*zymaxmax - 3*zyminmax - 2*zxyminmin - zxymaxmin - zxymaxmax - 2*zxyminmax;
  *z_p += 2*v*t1*u3;
  v = 2*zminmin - 2*zmaxmin + zxminmin + zxmaxmin;
  *z_p += 3*v*t2*u0;
  v = 2*zyminmin - 2*zymaxmin + zxyminmin + zxymaxmin;
  *z_p += 3*v*t2*u1;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 3*zxminmin - 3*zxmaxmin + 3*zxmaxmax + 3*zxminmax - 4*zyminmin + 4*zymaxmin + 2*zymaxmax - 2*zyminmax - 2*zxyminmin - 2*zxymaxmin - zxymaxmax - zxyminmax;
  *z_p += 3*v*t2*u2;
  v = 4*zminmin - 4*zmaxmin + 4*zmaxmax - 4*zminmax + 2*zxminmin + 2*zxmaxmin - 2*zxmaxmax - 2*zxminmax + 2*zyminmin - 2*zymaxmin - 2*zymaxmax + 2*zyminmax + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
  *z_p += 3*v*t2*u3;
  *z_p *= dt;

  return GSL_SUCCESS;
} /* bicubic_deriv_x() */

static int
bicubic_deriv_y(const void * vstate, const double xarr[], const double yarr[],
                const double zarr[], size_t xsize, size_t ysize,
                double x, double y,
                gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_p)
{
  bicubic_state_t *state = (bicubic_state_t *) vstate;

  double xmin, xmax, ymin, ymax;
  double zminmin, zminmax, zmaxmin, zmaxmax;
  double zxminmin, zxminmax, zxmaxmin, zxmaxmax;
  double zyminmin, zyminmax, zymaxmin, zymaxmax;
  double zxyminmin, zxyminmax, zxymaxmin, zxymaxmax;
  /* dx and dy are the size of the grid cell */
  double dx, dy;
  double dt, du;
  /* t and u are the positions within the grid cell at which we are
   * computing the interpolation, in units of grid cell size */
  double t, u;
  double t0, t1, t2, t3, u0, u1, u2;
  double v;
  size_t xi, yi;

  /* first compute the indices into the data arrays where we are interpolating */
  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  /* find the minimum and maximum values on the grid cell in each dimension */
  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, state)];
  zminmax = zarr[IDX2D(xi, yi + 1, state)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, state)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, state)];

  /* get the width and height of the grid cell */
  dx = xmax - xmin;
  dy = ymax - ymin;
  t = (x - xmin)/dx;
  u = (y - ymin)/dy;
  dt = 1./dx; /* partial t / partial x */
  du = 1./dy; /* partial u / partial y */

  zxminmin = state->zx[IDX2D(xi, yi, state)]/dt;
  zxminmax = state->zx[IDX2D(xi, yi + 1, state)]/dt;
  zxmaxmin = state->zx[IDX2D(xi + 1, yi, state)]/dt;
  zxmaxmax = state->zx[IDX2D(xi + 1, yi + 1, state)]/dt;
  zyminmin = state->zy[IDX2D(xi, yi, state)]/du;
  zyminmax = state->zy[IDX2D(xi, yi + 1, state)]/du;
  zymaxmin = state->zy[IDX2D(xi + 1, yi, state)]/du;
  zymaxmax = state->zy[IDX2D(xi + 1, yi + 1, state)]/du;
  zxyminmin = state->zxy[IDX2D(xi, yi, state)]/(dt*du);
  zxyminmax = state->zxy[IDX2D(xi, yi + 1, state)]/(dt*du);
  zxymaxmin = state->zxy[IDX2D(xi + 1, yi, state)]/(dt*du);
  zxymaxmax = state->zxy[IDX2D(xi + 1, yi + 1, state)]/(dt*du);

  t0 = 1;
  t1 = t;
  t2 = t*t;
  t3 = t*t2;
  u0 = 1;
  u1 = u;
  u2 = u*u;

  *z_p = 0;
  v = zyminmin;
  *z_p += v*t0*u0;
  v = -3*zminmin + 3*zminmax - 2*zyminmin - zyminmax;
  *z_p += 2*v*t0*u1;
  v = 2*zminmin-2*zminmax + zyminmin + zyminmax;
  *z_p += 3*v*t0*u2;
  v = zxyminmin;
  *z_p += v*t1*u0;
  v = -3*zxminmin + 3*zxminmax - 2*zxyminmin - zxyminmax;
  *z_p += 2*v*t1*u1;
  v = 2*zxminmin - 2*zxminmax + zxyminmin + zxyminmax;
  *z_p += 3*v*t1*u2;
  v = -3*zyminmin + 3*zymaxmin - 2*zxyminmin - zxymaxmin;
  *z_p += v*t2*u0;
  v = 9*zminmin - 9*zmaxmin + 9*zmaxmax - 9*zminmax + 6*zxminmin + 3*zxmaxmin - 3*zxmaxmax - 6*zxminmax + 6*zyminmin - 6*zymaxmin - 3*zymaxmax + 3*zyminmax + 4*zxyminmin + 2*zxymaxmin + zxymaxmax + 2*zxyminmax;
  *z_p += 2*v*t2*u1;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 4*zxminmin - 2*zxmaxmin + 2*zxmaxmax + 4*zxminmax - 3*zyminmin + 3*zymaxmin + 3*zymaxmax - 3*zyminmax - 2*zxyminmin - zxymaxmin - zxymaxmax - 2*zxyminmax;
  *z_p += 3*v*t2*u2;
  v = 2*zyminmin - 2*zymaxmin + zxyminmin + zxymaxmin;
  *z_p += v*t3*u0;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 3*zxminmin - 3*zxmaxmin + 3*zxmaxmax + 3*zxminmax - 4*zyminmin + 4*zymaxmin + 2*zymaxmax - 2*zyminmax - 2*zxyminmin - 2*zxymaxmin - zxymaxmax - zxyminmax;
  *z_p += 2*v*t3*u1;
  v = 4*zminmin - 4*zmaxmin + 4*zmaxmax - 4*zminmax + 2*zxminmin + 2*zxmaxmin - 2*zxmaxmax - 2*zxminmax + 2*zyminmin - 2*zymaxmin - 2*zymaxmax + 2*zyminmax + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
  *z_p += 3*v*t3*u2;
  *z_p *= du;

  return GSL_SUCCESS;
}

static int
bicubic_deriv_xx(const void * vstate, const double xarr[], const double yarr[],
                 const double zarr[], size_t xsize, size_t ysize,
                 double x, double y,
                 gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_pp)
{
  bicubic_state_t *state = (bicubic_state_t *) vstate;

  double xmin, xmax, ymin, ymax;
  double zminmin, zminmax, zmaxmin, zmaxmax;
  double zxminmin, zxminmax, zxmaxmin, zxmaxmax;
  double zyminmin, zyminmax, zymaxmin, zymaxmax;
  double zxyminmin, zxyminmax, zxymaxmin, zxymaxmax;

  double dx, dy; /* size of the grid cell */
  double dt, du;

  /*
   * t and u are the positions within the grid cell at which we are computing
   * the interpolation, in units of grid cell size
   */
  double t, u;
  double t0, t1, u0, u1, u2, u3;
  double v;
  size_t xi, yi;

  /* first compute the indices into the data arrays where we are interpolating */
  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  /* find the minimum and maximum values on the grid cell in each dimension */
  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, state)];
  zminmax = zarr[IDX2D(xi, yi + 1, state)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, state)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, state)];

  /* get the width and height of the grid cell */
  dx = xmax - xmin;
  dy = ymax - ymin;
  t = (x - xmin)/dx;
  u = (y - ymin)/dy;
  dt = 1./dx; /* partial t / partial x */
  du = 1./dy; /* partial u / partial y */

  zxminmin = state->zx[IDX2D(xi, yi, state)]/dt;
  zxminmax = state->zx[IDX2D(xi, yi + 1, state)]/dt;
  zxmaxmin = state->zx[IDX2D(xi + 1, yi, state)]/dt;
  zxmaxmax = state->zx[IDX2D(xi + 1, yi + 1, state)]/dt;
  zyminmin = state->zy[IDX2D(xi, yi, state)]/du;
  zyminmax = state->zy[IDX2D(xi, yi + 1, state)]/du;
  zymaxmin = state->zy[IDX2D(xi + 1, yi, state)]/du;
  zymaxmax = state->zy[IDX2D(xi + 1, yi + 1, state)]/du;
  zxyminmin = state->zxy[IDX2D(xi, yi, state)]/(dt*du);
  zxyminmax = state->zxy[IDX2D(xi, yi + 1, state)]/(dt*du);
  zxymaxmin = state->zxy[IDX2D(xi + 1, yi, state)]/(dt*du);
  zxymaxmax = state->zxy[IDX2D(xi + 1, yi + 1, state)]/(dt*du);

  t0 = 1;
  t1 = t;
  u0 = 1;
  u1 = u;
  u2 = u*u;
  u3 = u*u2;

  *z_pp = 0;
  v = -3*zminmin + 3*zmaxmin - 2*zxminmin - zxmaxmin;
  *z_pp += 2*v*t0*u0;
  v = -3*zyminmin + 3*zymaxmin - 2*zxyminmin - zxymaxmin;
  *z_pp += 2*v*t0*u1;
  v = 9*zminmin - 9*zmaxmin + 9*zmaxmax - 9*zminmax + 6*zxminmin + 3*zxmaxmin - 3*zxmaxmax - 6*zxminmax + 6*zyminmin - 6*zymaxmin - 3*zymaxmax + 3*zyminmax + 4*zxyminmin + 2*zxymaxmin + zxymaxmax + 2*zxyminmax;
  *z_pp += 2*v*t0*u2;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 4*zxminmin - 2*zxmaxmin + 2*zxmaxmax + 4*zxminmax - 3*zyminmin + 3*zymaxmin + 3*zymaxmax - 3*zyminmax - 2*zxyminmin - zxymaxmin - zxymaxmax - 2*zxyminmax;
  *z_pp += 2*v*t0*u3;
  v = 2*zminmin - 2*zmaxmin + zxminmin + zxmaxmin;
  *z_pp += 6*v*t1*u0;
  v = 2*zyminmin - 2*zymaxmin + zxyminmin + zxymaxmin;
  *z_pp += 6*v*t1*u1;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 3*zxminmin - 3*zxmaxmin + 3*zxmaxmax + 3*zxminmax - 4*zyminmin + 4*zymaxmin + 2*zymaxmax - 2*zyminmax - 2*zxyminmin - 2*zxymaxmin - zxymaxmax - zxyminmax;
  *z_pp += 6*v*t1*u2;
  v = 4*zminmin - 4*zmaxmin + 4*zmaxmax - 4*zminmax + 2*zxminmin + 2*zxmaxmin - 2*zxmaxmax - 2*zxminmax + 2*zyminmin - 2*zymaxmin - 2*zymaxmax + 2*zyminmax + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
  *z_pp += 6*v*t1*u3;
  *z_pp *= dt*dt;

  return GSL_SUCCESS;
}

static int
bicubic_deriv_xy(const void * vstate, const double xarr[], const double yarr[],
                 const double zarr[], size_t xsize, size_t ysize,
                 double x, double y,
                 gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_pp)
{
  bicubic_state_t *state = (bicubic_state_t *) vstate;

  double xmin, xmax, ymin, ymax;
  double zminmin, zminmax, zmaxmin, zmaxmax;
  double zxminmin, zxminmax, zxmaxmin, zxmaxmax;
  double zyminmin, zyminmax, zymaxmin, zymaxmax;
  double zxyminmin, zxyminmax, zxymaxmin, zxymaxmax;

  double dx, dy; /* size of the grid cell */
  double dt, du;

  /*
   * t and u are the positions within the grid cell at which we are computing
   * the interpolation, in units of grid cell size
   */
  double t, u;
  double t0, t1, t2, u0, u1, u2;
  double v;
  size_t xi, yi;

  /* first compute the indices into the data arrays where we are interpolating */
  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  /* find the minimum and maximum values on the grid cell in each dimension */
  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, state)];
  zminmax = zarr[IDX2D(xi, yi + 1, state)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, state)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, state)];

  /* get the width and height of the grid cell */
  dx = xmax - xmin;
  dy = ymax - ymin;
  t = (x - xmin)/dx;
  u = (y - ymin)/dy;
  dt = 1./dx; /* partial t / partial x */
  du = 1./dy; /* partial u / partial y */

  zxminmin = state->zx[IDX2D(xi, yi, state)]/dt;
  zxminmax = state->zx[IDX2D(xi, yi + 1, state)]/dt;
  zxmaxmin = state->zx[IDX2D(xi + 1, yi, state)]/dt;
  zxmaxmax = state->zx[IDX2D(xi + 1, yi + 1, state)]/dt;
  zyminmin = state->zy[IDX2D(xi, yi, state)]/du;
  zyminmax = state->zy[IDX2D(xi, yi + 1, state)]/du;
  zymaxmin = state->zy[IDX2D(xi + 1, yi, state)]/du;
  zymaxmax = state->zy[IDX2D(xi + 1, yi + 1, state)]/du;
  zxyminmin = state->zxy[IDX2D(xi, yi, state)]/(dt*du);
  zxyminmax = state->zxy[IDX2D(xi, yi + 1, state)]/(dt*du);
  zxymaxmin = state->zxy[IDX2D(xi + 1, yi, state)]/(dt*du);
  zxymaxmax = state->zxy[IDX2D(xi + 1, yi + 1, state)]/(dt*du);

  t0 = 1;
  t1 = t;
  t2 = t*t;
  u0 = 1;
  u1 = u;
  u2 = u*u;

  *z_pp = 0;
  v = zxyminmin;
  *z_pp += v*t0*u0;
  v = -3*zxminmin + 3*zxminmax - 2*zxyminmin - zxyminmax;
  *z_pp += 2*v*t0*u1;
  v = 2*zxminmin - 2*zxminmax + zxyminmin + zxyminmax;
  *z_pp += 3*v*t0*u2;
  v = -3*zyminmin + 3*zymaxmin - 2*zxyminmin - zxymaxmin;
  *z_pp += 2*v*t1*u0;
  v = 9*zminmin - 9*zmaxmin + 9*zmaxmax - 9*zminmax + 6*zxminmin + 3*zxmaxmin - 3*zxmaxmax - 6*zxminmax + 6*zyminmin - 6*zymaxmin - 3*zymaxmax + 3*zyminmax + 4*zxyminmin + 2*zxymaxmin + zxymaxmax + 2*zxyminmax;
  *z_pp += 4*v*t1*u1;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 4*zxminmin - 2*zxmaxmin + 2*zxmaxmax + 4*zxminmax - 3*zyminmin + 3*zymaxmin + 3*zymaxmax - 3*zyminmax - 2*zxyminmin - zxymaxmin - zxymaxmax - 2*zxyminmax;
  *z_pp += 6*v*t1*u2;
  v = 2*zyminmin - 2*zymaxmin + zxyminmin + zxymaxmin;
  *z_pp += 3*v*t2*u0;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 3*zxminmin - 3*zxmaxmin + 3*zxmaxmax + 3*zxminmax - 4*zyminmin + 4*zymaxmin + 2*zymaxmax - 2*zyminmax - 2*zxyminmin - 2*zxymaxmin - zxymaxmax - zxyminmax;
  *z_pp += 6*v*t2*u1;
  v = 4*zminmin - 4*zmaxmin + 4*zmaxmax - 4*zminmax + 2*zxminmin + 2*zxmaxmin - 2*zxmaxmax - 2*zxminmax + 2*zyminmin - 2*zymaxmin - 2*zymaxmax + 2*zyminmax + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
  *z_pp += 9*v*t2*u2;
  *z_pp *= dt*du;

  return GSL_SUCCESS;
}

static int
bicubic_deriv_yy(const void * vstate, const double xarr[], const double yarr[],
                 const double zarr[], size_t xsize, size_t ysize,
                 double x, double y,
                 gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_pp)
{
  bicubic_state_t *state = (bicubic_state_t *) vstate;

  double xmin, xmax, ymin, ymax;
  double zminmin, zminmax, zmaxmin, zmaxmax;
  double zxminmin, zxminmax, zxmaxmin, zxmaxmax;
  double zyminmin, zyminmax, zymaxmin, zymaxmax;
  double zxyminmin, zxyminmax, zxymaxmin, zxymaxmax;

  double dx, dy; /* size of the grid cell */
  double dt, du;

  /*
   * t and u are the positions within the grid cell at which we are computing
   * the interpolation, in units of grid cell size
   */
  double t, u;
  double t0, t1, t2, t3, u0, u1;
  double v;
  size_t xi, yi;

  /* first compute the indices into the data arrays where we are interpolating */
  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  /* find the minimum and maximum values on the grid cell in each dimension */
  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, state)];
  zminmax = zarr[IDX2D(xi, yi + 1, state)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, state)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, state)];

  /* get the width and height of the grid cell */
  dx = xmax - xmin;
  dy = ymax - ymin;
  t = (x - xmin)/dx;
  u = (y - ymin)/dy;
  dt = 1./dx; /* partial t / partial x */
  du = 1./dy; /* partial u / partial y */

  zxminmin = state->zx[IDX2D(xi, yi, state)]/dt;
  zxminmax = state->zx[IDX2D(xi, yi + 1, state)]/dt;
  zxmaxmin = state->zx[IDX2D(xi + 1, yi, state)]/dt;
  zxmaxmax = state->zx[IDX2D(xi + 1, yi + 1, state)]/dt;
  zyminmin = state->zy[IDX2D(xi, yi, state)]/du;
  zyminmax = state->zy[IDX2D(xi, yi + 1, state)]/du;
  zymaxmin = state->zy[IDX2D(xi + 1, yi, state)]/du;
  zymaxmax = state->zy[IDX2D(xi + 1, yi + 1, state)]/du;
  zxyminmin = state->zxy[IDX2D(xi, yi, state)]/(dt*du);
  zxyminmax = state->zxy[IDX2D(xi, yi + 1, state)]/(dt*du);
  zxymaxmin = state->zxy[IDX2D(xi + 1, yi, state)]/(dt*du);
  zxymaxmax = state->zxy[IDX2D(xi + 1, yi + 1, state)]/(dt*du);

  t0 = 1;
  t1 = t;
  t2 = t*t;
  t3 = t*t2;
  u0 = 1;
  u1 = u;

  *z_pp = 0;
  v = -3*zminmin + 3*zminmax - 2*zyminmin - zyminmax;
  *z_pp += 2*v*t0*u0;
  v = 2*zminmin-2*zminmax + zyminmin + zyminmax;
  *z_pp += 6*v*t0*u1;
  v = -3*zxminmin + 3*zxminmax - 2*zxyminmin - zxyminmax;
  *z_pp += 2*v*t1*u0;
  v = 2*zxminmin - 2*zxminmax + zxyminmin + zxyminmax;
  *z_pp += 6*v*t1*u1;
  v = 9*zminmin - 9*zmaxmin + 9*zmaxmax - 9*zminmax + 6*zxminmin + 3*zxmaxmin - 3*zxmaxmax - 6*zxminmax + 6*zyminmin - 6*zymaxmin - 3*zymaxmax + 3*zyminmax + 4*zxyminmin + 2*zxymaxmin + zxymaxmax + 2*zxyminmax;
  *z_pp += 2*v*t2*u0;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 4*zxminmin - 2*zxmaxmin + 2*zxmaxmax + 4*zxminmax - 3*zyminmin + 3*zymaxmin + 3*zymaxmax - 3*zyminmax - 2*zxyminmin - zxymaxmin - zxymaxmax - 2*zxyminmax;
  *z_pp += 6*v*t2*u1;
  v = -6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 3*zxminmin - 3*zxmaxmin + 3*zxmaxmax + 3*zxminmax - 4*zyminmin + 4*zymaxmin + 2*zymaxmax - 2*zyminmax - 2*zxyminmin - 2*zxymaxmin - zxymaxmax - zxyminmax;
  *z_pp += 2*v*t3*u0;
  v = 4*zminmin - 4*zmaxmin + 4*zmaxmax - 4*zminmax + 2*zxminmin + 2*zxmaxmin - 2*zxmaxmax - 2*zxminmax + 2*zyminmin - 2*zymaxmin - 2*zymaxmax + 2*zyminmax + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
  *z_pp += 6*v*t3*u1;
  *z_pp *= du*du;

  return GSL_SUCCESS;
}

static const gsl_interp2d_type bicubic_type = {
  "bicubic",
  4,
  &bicubic_alloc,
  &bicubic_init,
  &bicubic_eval,
  &bicubic_deriv_x,
  &bicubic_deriv_y,
  &bicubic_deriv_xx,
  &bicubic_deriv_xy,
  &bicubic_deriv_yy,
  &bicubic_free
};

const gsl_interp2d_type * gsl_interp2d_bicubic = &bicubic_type;

#undef IDX2D
