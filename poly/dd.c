/* interpolation/interp_poly.c
 * 
 * Copyright (C) 2001 DAN, HO-JIN
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

/* Modified for standalone use in polynomial directory, B.Gough 2001 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

int
gsl_poly_dd_init (double dd[], const double xa[], const double ya[],
                  size_t size)
{
  size_t i, j;

  /* Newton's divided differences */

  dd[0] = ya[0];

  for (j = size - 1; j >= 1; j--)
    {
      dd[j] = (ya[j] - ya[j - 1]) / (xa[j] - xa[j - 1]);
    }

  for (i = 2; i < size; i++)
    {
      for (j = size - 1; j >= i; j--)
        {
          dd[j] = (dd[j] - dd[j - 1]) / (xa[j] - xa[j - i]);
        }
    }

  return GSL_SUCCESS;
}

int
gsl_poly_dd_taylor (double c[], double xp, 
                    const double dd[], const double xa[], size_t size,
                    double w[])
{
  size_t i, j;

  for (i = 0; i < size; i++)
    {
      c[i] = 0.0;
      w[i] = 0.0;
    }

  w[size - 1] = 1.0;

  c[0] = dd[0];

  for (i = size - 1; i-- > 0;)
    {
      w[i] = -w[i + 1] * (xa[size - 2 - i] - xp);

      for (j = i + 1; j < size - 1; j++)
        {
          w[j] = w[j] - w[j + 1] * (xa[size - 2 - i] - xp);
        }

      for (j = i; j < size; j++)
        {
          c[j - i] += w[j] * dd[size - i - 1];
        }
    }

  return GSL_SUCCESS;
}

/*
gsl_poly_dd_hermite_init()
  Compute divided difference representation of data
for Hermite polynomial interpolation

Inputs: dd   - (output) array of size 2*size containing
               divided differences, dd[k] = f[z_0,z_1,...,z_k]
        za   - (output) array of size 2*size containing
               z values
        xa   - x data
        ya   - y data
        dya  - dy/dx data
        size - size of xa,ya,dya arrays

Return: success
*/

int
gsl_poly_dd_hermite_init (double dd[], double za[], const double xa[], const double ya[],
                          const double dya[], const size_t size)
{
  const size_t N = 2 * size;
  size_t i, j;

  /* Hermite divided differences */

  dd[0] = ya[0];

  /* compute: dd[j] = f[z_{j-1},z_j] for j \in [1,N-1] */
  for (j = 0; j < size; ++j)
    {
      za[2*j] = xa[j];
      za[2*j + 1] = xa[j];

      if (j != 0)
        {
          dd[2*j] = (ya[j] - ya[j - 1]) / (xa[j] - xa[j - 1]);
          dd[2*j - 1] = dya[j - 1];
        }
    }

  dd[N - 1] = dya[size - 1];

  for (i = 2; i < N; i++)
    {
      for (j = N - 1; j >= i; j--)
        {
          dd[j] = (dd[j] - dd[j - 1]) / (za[j] - za[j - i]);
        }
    }

  return GSL_SUCCESS;
} /* gsl_poly_dd_hermite_init() */
