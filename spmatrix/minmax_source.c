/* spmatrix/minmax_source.c
 * 
 * Copyright (C) 2019 Patrick Alken
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

int
FUNCTION (gsl_spmatrix, minmax) (const TYPE (gsl_spmatrix) * m,
                                 ATOMIC * min_out, ATOMIC * max_out)
{
  ATOMIC min, max;
  size_t n;

  if (m->nz == 0)
    {
      GSL_ERROR("matrix is empty", GSL_EINVAL);
    }

  min = m->data[0];
  max = m->data[0];

  for (n = 1; n < m->nz; ++n)
    {
      ATOMIC x = m->data[n];

      if (x < min)
        min = x;

      if (x > max)
        max = x;
    }

  *min_out = min;
  *max_out = max;

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_spmatrix, min_index) (const TYPE (gsl_spmatrix) * m,
                                    size_t * imin_out, size_t * jmin_out)
{
  if (m->nz == 0)
    {
      GSL_ERROR("matrix is empty", GSL_EINVAL);
    }
  else
    {
      ATOMIC * Ad = m->data;
      ATOMIC min = Ad[0];
      int imin = 0;
      int jmin = 0;
      size_t n;

      if (GSL_SPMATRIX_ISCOO(m))
        {
          imin = m->i[0];
          jmin = m->p[0];

          for (n = 1; n < m->nz; ++n)
            {
              ATOMIC x = m->data[n];

              if (x < min)
                {
                  min = x;
                  imin = m->i[n];
                  jmin = m->p[n];
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSC(m))
        {
          const int *Ap = m->p;
          int p;
          size_t j;

          for (j = 0; j < m->size2; ++j)
            {
              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  if (Ad[p] < min)
                    {
                      min = Ad[p];
                      imin = m->i[p];
                      jmin = (int) j;
                    }
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSR(m))
        {
          const int *Ap = m->p;
          int p;
          size_t i;

          for (i = 0; i < m->size1; ++i)
            {
              for (p = Ap[i]; p < Ap[i + 1]; ++p)
                {
                  if (Ad[p] < min)
                    {
                      min = Ad[p];
                      imin = (int) i;
                      jmin = m->i[p];
                    }
                }
            }
        }
      else
        {
          GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
        }

      *imin_out = (size_t) imin;
      *jmin_out = (size_t) jmin;

      return GSL_SUCCESS;
    }
}
