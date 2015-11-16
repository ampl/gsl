/* interpolation/bilinear.c
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
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>

#define IDX2D(i, j, xsize, ysize) ((j) * (xsize) + (i))

static int
bilinear_init(void * state, const double xa[], const double ya[],
              const double za[], size_t xsize, size_t ysize)
{
  return GSL_SUCCESS;
}

static int
bilinear_eval(const void * state, const double xarr[], const double yarr[],
              const double zarr[], size_t xsize, size_t ysize,
              double x, double y, gsl_interp_accel * xa,
              gsl_interp_accel * ya, double * z)
{
  double xmin, xmax, ymin, ymax, zminmin, zminmax, zmaxmin, zmaxmax;
  double dx, dy;
  double t, u;
  size_t xi, yi;

  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, xsize, ysize)];
  zminmax = zarr[IDX2D(xi, yi + 1, xsize, ysize)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, xsize, ysize)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, xsize, ysize)];
  dx = xmax - xmin;
  dy = ymax - ymin;
  t = (x - xmin)/dx;
  u = (y - ymin)/dy;
  *z = (1.-t)*(1.-u)*zminmin + t*(1.-u)*zmaxmin + (1.-t)*u*zminmax + t*u*zmaxmax;

  return GSL_SUCCESS;
}

static int
bilinear_deriv_x(const void * state, const double xarr[],
                 const double yarr[], const double zarr[],
                 size_t xsize, size_t ysize, double x, double y,
                 gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_p)
{
  double xmin, xmax, ymin, ymax, zminmin, zminmax, zmaxmin, zmaxmax;
  double dx, dy;
  double dt, u;
  size_t xi, yi;

  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, xsize, ysize)];
  zminmax = zarr[IDX2D(xi, yi + 1, xsize, ysize)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, xsize, ysize)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, xsize, ysize)];
  dx = xmax - xmin;
  dy = ymax - ymin;
  dt = 1./dx; /* partial t / partial x */
  u = (y - ymin)/dy;
  *z_p = dt*(-(1.-u)*zminmin + (1.-u)*zmaxmin - u*zminmax + u*zmaxmax);

  return GSL_SUCCESS;
}

static int
bilinear_deriv_y(const void * state, const double xarr[],
                 const double yarr[], const double zarr[],
                 size_t xsize, size_t ysize, double x, double y,
                 gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_p)
{
  double xmin, xmax, ymin, ymax, zminmin, zminmax, zmaxmin, zmaxmax;
  double dx, dy;
  double t, du;
  size_t xi, yi;

  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, xsize, ysize)];
  zminmax = zarr[IDX2D(xi, yi + 1, xsize, ysize)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, xsize, ysize)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, xsize, ysize)];
  dx = xmax - xmin;
  dy = ymax - ymin;
  t = (x - xmin)/dx;
  du = 1./dy; /* partial u / partial y */
  *z_p = du*(-(1.-t)*zminmin - t*zmaxmin + (1.-t)*zminmax + t*zmaxmax);

  return GSL_SUCCESS;
}

static int
bilinear_deriv2(const void * state, const double xarr[],
                const double yarr[], const double zarr[],
                size_t xsize, size_t ysize, double x, double y,
                gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_pp)
{
  *z_pp = 0.0;
  return GSL_SUCCESS;
}

static int
bilinear_derivxy(const void * state, const double xarr[],
                 const double yarr[], const double zarr[],
                 size_t xsize, size_t ysize, double x, double y,
                 gsl_interp_accel * xa, gsl_interp_accel * ya, double * z_pp)
{
  double xmin, xmax, ymin, ymax, zminmin, zminmax, zmaxmin, zmaxmax;
  double dx, dy;
  double dt, du;
  size_t xi, yi;

  if (xa != NULL)
    xi = gsl_interp_accel_find(xa, xarr, xsize, x);
  else
    xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);

  if (ya != NULL)
    yi = gsl_interp_accel_find(ya, yarr, ysize, y);
  else
    yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);

  xmin = xarr[xi];
  xmax = xarr[xi + 1];
  ymin = yarr[yi];
  ymax = yarr[yi + 1];
  zminmin = zarr[IDX2D(xi, yi, xsize, ysize)];
  zminmax = zarr[IDX2D(xi, yi + 1, xsize, ysize)];
  zmaxmin = zarr[IDX2D(xi + 1, yi, xsize, ysize)];
  zmaxmax = zarr[IDX2D(xi + 1, yi + 1, xsize, ysize)];
  dx = xmax - xmin;
  dy = ymax - ymin;
  dt = 1./dx; /* partial t / partial x */
  du = 1./dy; /* partial u / partial y */
  *z_pp = dt*du*(zminmin-zmaxmin-zminmax+zmaxmax);

  return GSL_SUCCESS;
}

static const gsl_interp2d_type bilinear_type = {
  "bilinear",
  2,
  NULL,
  &bilinear_init,
  &bilinear_eval,
  &bilinear_deriv_x,
  &bilinear_deriv_y,
  &bilinear_deriv2,
  &bilinear_derivxy,
  &bilinear_deriv2,
  NULL
};

const gsl_interp2d_type * gsl_interp2d_bilinear = &bilinear_type;

#undef IDX2D
