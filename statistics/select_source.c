/* statistics/select_source.c
 * 
 * Copyright (C) 2018 Patrick Alken
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

#define SWAP(a,b) do { tmp = b ; b = a ; a = tmp ; } while(0)

/*
gsl_stats_select()
  Select k-th largest element from an unsorted array using
quickselect algorithm

Inputs: data   - unsorted array containing the observations
        stride - stride
        n      - length of 'data'
        k      - desired element in [0,n-1]

Return: k-th largest element of data[]
*/

BASE
FUNCTION(gsl_stats,select) (BASE data[],
                            const size_t stride,
                            const size_t n,
                            const size_t k)
{
  if (n == 0)
    {
      GSL_ERROR_VAL("array size must be positive", GSL_EBADLEN, 0.0);
    }
  else
    {
      size_t left = 0;
      size_t right = n - 1;
      size_t mid, i, j;
      BASE pivot, tmp;

      while (1)
        {
          if (right <= left + 1)
            {
              if (right == left + 1 && data[right * stride] < data[left * stride])
                {
                  SWAP(data[left * stride], data[right * stride]);
                }

              return data[k * stride];
            }
          else
            {
              mid = (left + right) >> 1;
              SWAP(data[mid * stride], data[(left + 1) * stride]);

              if (data[left * stride] > data[right * stride])
                {
                  SWAP(data[left * stride], data[right * stride]);
                }

              if (data[(left + 1) * stride] > data[right * stride])
                {
                  SWAP(data[(left + 1) * stride], data[right * stride]);
                }

              if (data[left * stride] > data[(left + 1) * stride])
                {
                  SWAP(data[left * stride], data[(left + 1) * stride]);
                }

              i = left + 1;
              j = right;
              pivot = data[(left + 1) * stride];

              while (1)
                {
                  do i++; while (data[i * stride] < pivot);
                  do j--; while (data[j * stride] > pivot);

                  if (j < i)
                    break;

                  SWAP(data[i * stride], data[j * stride]);
                }

              data[(left + 1) * stride] = data[j * stride];
              data[j * stride] = pivot;

              if (j >= k)
                right = j - 1;

              if (j <= k)
                left = i;
            }
        }

      /* will never get here */
      GSL_ERROR_VAL("select error", GSL_FAILURE, 0.0);
    }
}
