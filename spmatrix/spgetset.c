/* spgetset.c
 * 
 * Copyright (C) 2012 Patrick Alken
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
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>

#include "avl.c"

static void *tree_find(const gsl_spmatrix *m, const size_t i, const size_t j);

double
gsl_spmatrix_get(const gsl_spmatrix *m, const size_t i, const size_t j)
{
  if (i >= m->size1)
    {
      GSL_ERROR_VAL("first index out of range", GSL_EINVAL, 0.0);
    }
  else if (j >= m->size2)
    {
      GSL_ERROR_VAL("second index out of range", GSL_EINVAL, 0.0);
    }
  else if (m->nz == 0)
    {
      /* no non-zero elements added to matrix */
      return (0.0);
    }
  else
    {
      if (GSL_SPMATRIX_ISTRIPLET(m))
        {
          /* traverse binary tree to search for (i,j) element */
          void *ptr = tree_find(m, i, j);
          double x = ptr ? *(double *) ptr : 0.0;

          return x;
        }
      else if (GSL_SPMATRIX_ISCCS(m))
        {
          const size_t *mi = m->i;
          const size_t *mp = m->p;
          size_t p;

          /* loop over column j and search for row index i */
          for (p = mp[j]; p < mp[j + 1]; ++p)
            {
              if (mi[p] == i)
                return m->data[p];
            }
        }
      else if (GSL_SPMATRIX_ISCRS(m))
        {
          const size_t *mj = m->i;
          const size_t *mp = m->p;
          size_t p;

          /* loop over row i and search for column index j */
          for (p = mp[i]; p < mp[i + 1]; ++p)
            {
              if (mj[p] == j)
                return m->data[p];
            }
        }
      else
        {
          GSL_ERROR_VAL("unknown sparse matrix type", GSL_EINVAL, 0.0);
        }

      /* element not found; return 0 */
      return 0.0;
    }
} /* gsl_spmatrix_get() */

/*
gsl_spmatrix_set()
  Add an element to a matrix in triplet form

Inputs: m - spmatrix
        i - row index
        j - column index
        x - matrix value
*/

int
gsl_spmatrix_set(gsl_spmatrix *m, const size_t i, const size_t j,
                 const double x)
{
  if (!GSL_SPMATRIX_ISTRIPLET(m))
    {
      GSL_ERROR("matrix not in triplet representation", GSL_EINVAL);
    }
  else if (x == 0.0)
    {
      /* traverse binary tree to search for (i,j) element */
      void *ptr = tree_find(m, i, j);

      /*
       * just set the data element to 0; it would be easy to
       * delete the node from the tree with avl_delete(), but
       * we'd also have to delete it from the data arrays which
       * is less simple
       */
      if (ptr != NULL)
        *(double *) ptr = 0.0;

      return GSL_SUCCESS;
    }
  else
    {
      int s = GSL_SUCCESS;
      void *ptr;

      /* check if matrix needs to be realloced */
      if (m->nz >= m->nzmax)
        {
          s = gsl_spmatrix_realloc(2 * m->nzmax, m);
          if (s)
            return s;
        }

      /* store the triplet (i, j, x) */
      m->i[m->nz] = i;
      m->p[m->nz] = j;
      m->data[m->nz] = x;

      ptr = avl_insert(m->tree_data->tree, &m->data[m->nz]);
      if (ptr != NULL)
        {
          /* found duplicate entry (i,j), replace with new x */
          *((double *) ptr) = x;
        }
      else
        {
          /* no duplicate (i,j) found, update indices as needed */

          /* increase matrix dimensions if needed */
          m->size1 = GSL_MAX(m->size1, i + 1);
          m->size2 = GSL_MAX(m->size2, j + 1);

          ++(m->nz);
        }

      return s;
    }
} /* gsl_spmatrix_set() */

double *
gsl_spmatrix_ptr(gsl_spmatrix *m, const size_t i, const size_t j)
{
  if (i >= m->size1)
    {
      GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
    }
  else if (j >= m->size2)
    {
      GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
    }
  else
    {
      if (GSL_SPMATRIX_ISTRIPLET(m))
        {
          /* traverse binary tree to search for (i,j) element */
          void *ptr = tree_find(m, i, j);
          return (double *) ptr;
        }
      else if (GSL_SPMATRIX_ISCCS(m))
        {
          const size_t *mi = m->i;
          const size_t *mp = m->p;
          size_t p;

          /* loop over column j and search for row index i */
          for (p = mp[j]; p < mp[j + 1]; ++p)
            {
              if (mi[p] == i)
                return &(m->data[p]);
            }
        }
      else if (GSL_SPMATRIX_ISCRS(m))
        {
          const size_t *mj = m->i;
          const size_t *mp = m->p;
          size_t p;

          /* loop over row i and search for column index j */
          for (p = mp[i]; p < mp[i + 1]; ++p)
            {
              if (mj[p] == j)
                return &(m->data[p]);
            }
        }
      else
        {
          GSL_ERROR_NULL("unknown sparse matrix type", GSL_EINVAL);
        }

      /* element not found; return 0 */
      return NULL;
    }
} /* gsl_spmatrix_ptr() */

/*
tree_find()
  Find node in tree corresponding to matrix entry (i,j). Adapted
from avl_find()

Inputs: m - spmatrix
        i - row index
        j - column index

Return: pointer to tree node data if found, NULL if not found
*/

static void *
tree_find(const gsl_spmatrix *m, const size_t i, const size_t j)
{
  const struct avl_table *tree = (struct avl_table *) m->tree_data->tree;
  const struct avl_node *p;

  for (p = tree->avl_root; p != NULL; )
    {
      size_t n = (double *) p->avl_data - m->data;
      size_t pi = m->i[n];
      size_t pj = m->p[n];
      int cmp = gsl_spmatrix_compare_idx(i, j, pi, pj);

      if (cmp < 0)
        p = p->avl_link[0];
      else if (cmp > 0)
        p = p->avl_link[1];
      else /* |cmp == 0| */
        return p->avl_data;
    }

  return NULL;
} /* tree_find() */
