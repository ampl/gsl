/* vector/copy_complex_source.c
 * 
 * Copyright (C) 2021 Patrick Alken
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
FUNCTION (gsl_vector, conj_memcpy) (TYPE (gsl_vector) * dest,
                                    const TYPE (gsl_vector) * src)
{
  const size_t src_size = src->size;
  const size_t dest_size = dest->size;

  if (src_size != dest_size)
    {
      GSL_ERROR ("vector lengths are not equal", GSL_EBADLEN);
    }

  {
    const size_t src_stride = src->stride;
    const size_t dest_stride = dest->stride;
    size_t j;

    for (j = 0; j < src_size; j++)
      {
        size_t e1 = 2 * j * dest_stride;
        size_t e2 = 2 * j * src_stride;

        dest->data[e1] = src->data[e2];
        dest->data[e1 + 1] = -src->data[e2 + 1];
      }
  }

  return GSL_SUCCESS;
}
