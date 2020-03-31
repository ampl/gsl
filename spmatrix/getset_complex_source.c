/* spmatrix/getset_complex_source.c
 * 
 * Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Patrick Alken
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

static void * FUNCTION (tree, find) (const TYPE (gsl_spmatrix) * m, const size_t i, const size_t j);

BASE
FUNCTION (gsl_spmatrix, get) (const TYPE (gsl_spmatrix) * m, const size_t i, const size_t j)
{
  BASE zero = ZERO;

  if (i >= m->size1)
    {
      GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
    }
  else if (j >= m->size2)
    {
      GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
    }
  else if (m->nz == 0)
    {
      /* no non-zero elements added to matrix */
      return zero;
    }
  else
    {
      if (GSL_SPMATRIX_ISCOO(m))
        {
          /* traverse binary tree to search for (i,j) element */
          void *ptr = FUNCTION (tree, find) (m, i, j);
          BASE x = ptr ? *(BASE *) ptr : zero;
          return x;
        }
      else if (GSL_SPMATRIX_ISCSC(m))
        {
          const int *mi = m->i;
          const int *mp = m->p;
          int p;

          /* loop over column j and search for row index i */
          for (p = mp[j]; p < mp[j + 1]; ++p)
            {
              if (mi[p] == (int) i)
                return *(BASE *) &m->data[2 * p];
            }
        }
      else if (GSL_SPMATRIX_ISCSR(m))
        {
          const int *mj = m->i;
          const int *mp = m->p;
          int p;

          /* loop over row i and search for column index j */
          for (p = mp[i]; p < mp[i + 1]; ++p)
            {
              if (mj[p] == (int) j)
                return *(BASE *) &m->data[2 * p];
            }
        }
      else
        {
          GSL_ERROR_VAL("unknown sparse matrix type", GSL_EINVAL, zero);
        }

      /* element not found; return 0 */
      return zero;
    }
}

int
FUNCTION (gsl_spmatrix, set) (TYPE (gsl_spmatrix) * m, const size_t i,
                              const size_t j, const BASE x)
{
  if (!GSL_SPMATRIX_ISCOO(m))
    {
      GSL_ERROR("matrix not in COO representation", GSL_EINVAL);
    }
  else if (!(m->spflags & GSL_SPMATRIX_FLG_GROW) && (i >= m->size1 || j >= m->size2))
    {
      GSL_ERROR ("indices out of range", GSL_EINVAL);
    }
  else if (m->spflags & GSL_SPMATRIX_FLG_FIXED)
    {
      /*
       * sparsity pattern is fixed - no elements may be inserted, but existing
       * elements can be replaced
       */
      void * ptr = FUNCTION (tree, find) (m, i, j);

      if (ptr == NULL)
        {
          GSL_ERROR("attempt to add new matrix element to fixed sparsity pattern", GSL_EINVAL);
        }
      else
        {
          *(BASE *) ptr = x;
        }

      return GSL_SUCCESS;
    }
  else
    {
      int status = GSL_SUCCESS;
      void *ptr;

      /* check if matrix needs to be reallocated */
      if (m->nz >= m->nzmax)
        {
          status = FUNCTION (gsl_spmatrix, realloc) (2 * m->nzmax, m);
          if (status)
            return status;
        }

      /* store the triplet (i, j, x) */
      m->i[m->nz] = i;
      m->p[m->nz] = j;
      m->data[2 * m->nz] = GSL_REAL (x);
      m->data[2 * m->nz + 1] = GSL_IMAG (x);

      ptr = gsl_bst_insert(&m->data[2 * m->nz], m->tree);
      if (ptr != NULL)
        {
          /* found duplicate entry (i,j), replace with new x */
          *((BASE *) ptr) = x;
        }
      else
        {
          /* no duplicate (i,j) found */

          /* increase matrix dimensions if needed */
          if (m->spflags & GSL_SPMATRIX_FLG_GROW)
            {
              m->size1 = GSL_MAX(m->size1, i + 1);
              m->size2 = GSL_MAX(m->size2, j + 1);
            }

          ++(m->nz);
        }

      return status;
    }
}

BASE *
FUNCTION (gsl_spmatrix, ptr) (const TYPE (gsl_spmatrix) * m, const size_t i, const size_t j)
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
      if (GSL_SPMATRIX_ISCOO(m))
        {
          /* traverse binary tree to search for (i,j) element */
          void *ptr = FUNCTION (tree, find) (m, i, j);
          return (BASE *) ptr;
        }
      else if (GSL_SPMATRIX_ISCSC(m))
        {
          const int *mi = m->i;
          const int *mp = m->p;
          int p;

          /* loop over column j and search for row index i */
          for (p = mp[j]; p < mp[j + 1]; ++p)
            {
              if (mi[p] == (int) i)
                return (BASE *) &(m->data[2 * p]);
            }
        }
      else if (GSL_SPMATRIX_ISCSR(m))
        {
          const int *mj = m->i;
          const int *mp = m->p;
          int p;

          /* loop over row i and search for column index j */
          for (p = mp[i]; p < mp[i + 1]; ++p)
            {
              if (mj[p] == (int) j)
                return (BASE *) &(m->data[2 * p]);
            }
        }
      else
        {
          GSL_ERROR_NULL("unknown sparse matrix type", GSL_EINVAL);
        }

      /* element not found; return NULL */
      return NULL;
    }
}

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
FUNCTION (tree, find) (const TYPE (gsl_spmatrix) * m, const size_t i, const size_t j)
{
  const gsl_bst_avl_table *table = (const gsl_bst_avl_table *) &(m->tree->table.avl_table);
  struct gsl_bst_avl_node *p;

  for (p = table->avl_root; p != NULL; )
    {
      size_t n = ((ATOMIC *) p->avl_data - m->data) >> 1;
      int pi = m->i[n];
      int pj = m->p[n];
      int cmp = GSL_SPMATRIX_COMPARE_ROWCOL(m, (int) i, (int) j, pi, pj);

      if (cmp < 0)
        p = p->avl_link[0];
      else if (cmp > 0)
        p = p->avl_link[1];
      else /* |cmp == 0| */
        return p->avl_data;
    }

  return NULL;
}
