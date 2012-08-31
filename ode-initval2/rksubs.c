/* ode-initval2/rksubs.c
 * 
 * Copyright (C) 2008, 2009, 2010 Tuomo Keskitalo
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

static int
rksubs (double y[], const double h, const double y0[], const double fY[],
        const double b[], const size_t stage, const size_t dim)
{
  /* The final substitution step in Runge-Kutta equation:
     Calculates new values y by substituting the value of step size
     (h), current initial values of y (y0), function f values at
     Runge-Kutta points (fY), Runge-Kutta b-coefficients (b) and
     method stage (stage) into the equation

     y = y0 + h * sum j=1..stage (b_j * fY_j)

     dim is the number of ODEs.
   */

  size_t i, j;

  for (i = 0; i < dim; i++)
    {
      y[i] = 0.0;

      for (j = 0; j < stage; j++)
        y[i] += b[j] * fY[j * dim + i];
    }

  for (i = 0; i < dim; i++)
    {
      y[i] *= h;
      y[i] += y0[i];
    }

  return GSL_SUCCESS;
}
