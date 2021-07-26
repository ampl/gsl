/* matrix/swap_complex_source.c
 * 
 * Copyright (C) 2020 Patrick Alken
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
FUNCTION (gsl_matrix, conjtrans_memcpy) (TYPE (gsl_matrix) * dest, 
                                         const TYPE (gsl_matrix) * src)
{
  const size_t src_size1 = src->size1;
  const size_t src_size2 = src->size2;
  const size_t dest_size1 = dest->size1;
  const size_t dest_size2 = dest->size2;
  size_t i;

  if (dest_size2 != src_size1 || dest_size1 != src_size2)
    {
      GSL_ERROR ("dimensions of dest matrix must be transpose of src matrix", 
                 GSL_EBADLEN);
    }

  for (i = 0; i < dest_size1; i++)
    {
      size_t j;

      for (j = 0 ; j < dest_size2; j++) 
        {
          size_t e1 = (i *  dest->tda + j) * 2;
          size_t e2 = (j *  src->tda + i) * 2;

          dest->data[e1] = src->data[e2];
          dest->data[e1 + 1] = -src->data[e2 + 1];
        }
    }

  return GSL_SUCCESS;
}
