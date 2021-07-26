/* spmatrix/init_source.c
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

static int FUNCTION(compare, func) (const void * pa, const void * pb, void * param);
static int FUNCTION (spmatrix, pool_init) (TYPE (gsl_spmatrix) * m);
static int FUNCTION (spmatrix, pool_free) (TYPE (gsl_spmatrix) * m);
static void * FUNCTION (spmatrix, malloc) (size_t size, void * params);
static void FUNCTION (spmatrix, free) (void * block, void * params);

static const gsl_bst_allocator FUNCTION(spmatrix, allocator) =
{
  FUNCTION (spmatrix, malloc),
  FUNCTION (spmatrix, free)
};

TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, alloc) (const size_t n1, const size_t n2)
{
  const double density = 0.1; /* estimate */
  size_t nzmax = (size_t) floor(n1 * n2 * density);

  if (nzmax == 0)
    nzmax = 10;

  return FUNCTION(gsl_spmatrix, alloc_nzmax) (n1, n2, nzmax, GSL_SPMATRIX_COO);
}

TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, alloc_nzmax) (const size_t n1, const size_t n2,
                                      const size_t nzmax, const int sptype)
{
  TYPE(gsl_spmatrix) * m;

  if (n1 == 0)
    {
      GSL_ERROR_NULL ("matrix dimension n1 must be positive integer",
                      GSL_EINVAL);
    }
  else if (n2 == 0)
    {
      GSL_ERROR_NULL ("matrix dimension n2 must be positive integer",
                      GSL_EINVAL);
    }

  m = calloc(1, sizeof(TYPE(gsl_spmatrix)));
  if (!m)
    {
      GSL_ERROR_NULL("failed to allocate space for spmatrix struct",
                     GSL_ENOMEM);
    }

  m->size1 = n1;
  m->size2 = n2;
  m->nz = 0;
  m->nzmax = GSL_MAX(nzmax, 1);
  m->sptype = sptype;

  if (n1 == 1 && n2 == 1)
    m->spflags = GSL_SPMATRIX_FLG_GROW; /* allow matrix size to grow */
  else
    m->spflags = 0;

  m->i = malloc(m->nzmax * sizeof(int));
  if (!m->i)
    {
      FUNCTION(gsl_spmatrix, free) (m);
      GSL_ERROR_NULL("failed to allocate space for row indices",
                     GSL_ENOMEM);
    }

  m->work.work_void = malloc(GSL_MAX(n1, n2) * MULTIPLICITY *
                             GSL_MAX(sizeof(int), sizeof(ATOMIC)));
  if (!m->work.work_void)
    {
      FUNCTION(gsl_spmatrix, free) (m);
      GSL_ERROR_NULL("failed to allocate space for work", GSL_ENOMEM);
    }

  if (sptype == GSL_SPMATRIX_COO)
    {
      m->tree = gsl_bst_alloc(gsl_bst_avl, &FUNCTION(spmatrix, allocator), FUNCTION (compare, func), (void *) m);
      if (!m->tree)
        {
          FUNCTION(gsl_spmatrix, free) (m);
          GSL_ERROR_NULL("failed to allocate space for binary tree",
                         GSL_ENOMEM);
        }

      m->node_size = gsl_bst_node_size(m->tree);

      /* initialize memory pool */
      FUNCTION (spmatrix, pool_init) (m);

      m->p = malloc(m->nzmax * sizeof(int));
      if (!m->p)
        {
          FUNCTION(gsl_spmatrix, free) (m);
          GSL_ERROR_NULL("failed to allocate space for column indices",
                         GSL_ENOMEM);
        }
    }
  else if (sptype == GSL_SPMATRIX_CSC)
    {
      m->p = malloc((n2 + 1) * sizeof(int));
      if (!m->p)
        {
          FUNCTION(gsl_spmatrix, free) (m);
          GSL_ERROR_NULL("failed to allocate space for column pointers",
                         GSL_ENOMEM);
        }
    }
  else if (sptype == GSL_SPMATRIX_CSR)
    {
      m->p = malloc((n1 + 1) * sizeof(int));
      if (!m->p)
        {
          FUNCTION(gsl_spmatrix, free) (m);
          GSL_ERROR_NULL("failed to allocate space for row pointers",
                         GSL_ENOMEM);
        }
    }

  m->data = malloc(m->nzmax * MULTIPLICITY * sizeof (ATOMIC));
  if (!m->data)
    {
      FUNCTION(gsl_spmatrix, free) (m);
      GSL_ERROR_NULL("failed to allocate space for data",
                     GSL_ENOMEM);
    }

  return m;
}

void
FUNCTION (gsl_spmatrix, free) (TYPE (gsl_spmatrix) * m)
{
  if (m->i)
    free(m->i);

  if (m->p)
    free(m->p);

  if (m->data)
    free(m->data);

  if (m->work.work_void)
    free(m->work.work_void);

  /* binary tree should be freed before pool */
  if (m->tree)
    gsl_bst_free(m->tree);

  FUNCTION (spmatrix, pool_free) (m);

  free(m);
}

int
FUNCTION (gsl_spmatrix, realloc) (const size_t nzmax, TYPE (gsl_spmatrix) * m)
{
  int status = GSL_SUCCESS;
  void *ptr;
  ATOMIC * ptr_atomic;

  if (nzmax < m->nz)
    {
      GSL_ERROR("new nzmax is less than current nz", GSL_EINVAL);
    }

  ptr = realloc(m->i, nzmax * sizeof(int));
  if (!ptr)
    {
      GSL_ERROR("failed to allocate space for row indices", GSL_ENOMEM);
    }

  m->i = (int *) ptr;

  if (GSL_SPMATRIX_ISCOO(m))
    {
      ptr = realloc(m->p, nzmax * sizeof(int));
      if (!ptr)
        {
          GSL_ERROR("failed to allocate space for column indices", GSL_ENOMEM);
        }

      m->p = (int *) ptr;
    }

  ptr_atomic = realloc(m->data, nzmax * MULTIPLICITY * sizeof (ATOMIC));
  if (!ptr_atomic)
    {
      GSL_ERROR("failed to allocate space for data", GSL_ENOMEM);
    }

  if (GSL_SPMATRIX_ISCOO(m))
    {
      const size_t nnew = nzmax - m->nz; /* number of new nodes to allocate in memory pool */
      gsl_spmatrix_pool * node;

      /*
       * if the address of m->data has moved, traverse binary tree and change
       * all node data to point to the new array 'ptr_atomic'
       */
      if (ptr_atomic != m->data)
        {
          gsl_bst_trav trav;

          ptr = gsl_bst_trav_first (&trav, m->tree);
          while (ptr != NULL)
            {
              ATOMIC * q = (ATOMIC *) ptr;
              size_t idx = q - m->data;

              gsl_bst_trav_replace (&trav, &ptr_atomic[idx]);
              ptr = gsl_bst_trav_next (&trav);
            }
        }

      /* allocate a new block in the memory pool to accomodate the additional nodes */
      node = malloc(sizeof(gsl_spmatrix_pool));
      if (!node)
        {
          GSL_ERROR ("failed to allocate space for memory pool node", GSL_ENOMEM);
        }

      node->block_ptr = malloc(nnew * m->node_size);
      if (!node->block_ptr)
        {
          GSL_ERROR ("failed to allocate space for memory pool block", GSL_ENOMEM);
        }

      node->free_slot = (unsigned char *) node->block_ptr;

      /* insert node at beginning of pool linked list */
      node->next = m->pool;
      m->pool = node;
    }

  m->data = ptr_atomic;
  m->nzmax = nzmax;

  return status;
}

size_t
FUNCTION (gsl_spmatrix, nnz) (const TYPE (gsl_spmatrix) * m)
{
  return m->nz;
}

const char *
FUNCTION (gsl_spmatrix, type) (const TYPE (gsl_spmatrix) * m)
{
  if (GSL_SPMATRIX_ISCOO(m))
    return "COO";
  else if (GSL_SPMATRIX_ISCSR(m))
    return "CSR";
  else if (GSL_SPMATRIX_ISCSC(m))
    return "CSC";
  else
    return "unknown";
}

int
FUNCTION (gsl_spmatrix, set_zero) (TYPE (gsl_spmatrix) * m)
{
  m->nz = 0;

  if (m->tree != NULL)
    {
      /* reset tree to empty state */
      gsl_bst_empty(m->tree);

      /* reset memory pool */
      FUNCTION (spmatrix, pool_free) (m);
      FUNCTION (spmatrix, pool_init) (m);
    }

  return GSL_SUCCESS;
}

/*
gsl_spmatrix_tree_rebuild()
  When reading a triplet matrix from disk, or when
copying a triplet matrix, it is necessary to rebuild the
binary tree for element searches.

Inputs: m - triplet matrix
*/

int
FUNCTION (gsl_spmatrix, tree_rebuild) (TYPE (gsl_spmatrix) * m)
{
  if (!GSL_SPMATRIX_ISCOO(m))
    {
      GSL_ERROR("matrix must be in COO format", GSL_EINVAL);
    }
  else
    {
      size_t n;

      /* empty the tree */
      gsl_bst_empty(m->tree);

      /* reset memory pool */
      FUNCTION (spmatrix, pool_free) (m);
      FUNCTION (spmatrix, pool_init) (m);

      /* re-insert all tree elements */
      for (n = 0; n < m->nz; ++n)
        {
          void *ptr = gsl_bst_insert(&m->data[MULTIPLICITY * n], m->tree);
          if (ptr != NULL)
            {
              GSL_ERROR("detected duplicate entry", GSL_EINVAL);
            }
        }

      return GSL_SUCCESS;
    }
}

/*
compare_func()
  Comparison function for searching binary tree in triplet
representation.

To detect duplicate elements in the tree, we want to determine
if there already exists an entry for (i,j) in the tree. Since
the actual tree node stores only the data elements data[n],
we will do pointer arithmetic to get from the given data[n]
to the row/column indices i[n] and j[n].

This compare function will sort the tree first by row i,
and for equal rows, it will then sort by column j

Inputs: pa    - element 1 for comparison (ATOMIC *) 
        pb    - element 2 for comparison (ATOMIC *)
        param - parameter (gsl_spmatrix)

Return:
  -1 if pa < pb: (ia,ja) < (ib,jb)
  +1 if pa > pb: (ia,ja) > (ib,jb)
   0 if pa = pb: (ia,ja) == (ib,jb)
*/

static int
FUNCTION(compare, func) (const void * pa, const void * pb, void * param)
{
  TYPE (gsl_spmatrix) * m = (TYPE (gsl_spmatrix) *) param;

  /* pointer arithmetic to find indices in data array */
#if (MULTIPLICITY == 1)
  const size_t idxa = (const ATOMIC *) pa - m->data;
  const size_t idxb = (const ATOMIC *) pb - m->data;
#else
  const size_t idxa = ((const ATOMIC *) pa - m->data) >> 1;
  const size_t idxb = ((const ATOMIC *) pb - m->data) >> 1;
#endif

  return GSL_SPMATRIX_COMPARE_ROWCOL(m, m->i[idxa], m->p[idxa], m->i[idxb], m->p[idxb]);
}

/* initialize m->pool by preallocating nzmax tree nodes */
static int
FUNCTION (spmatrix, pool_init) (TYPE (gsl_spmatrix) * m)
{
  m->pool = malloc(sizeof(gsl_spmatrix_pool));
  if (m->pool == NULL)
    {
      GSL_ERROR ("failed to allocate space for memory pool", GSL_ENOMEM);
    }

  m->pool->block_ptr = malloc(m->nzmax * m->node_size);
  if (m->pool->block_ptr == NULL)
    {
      GSL_ERROR ("failed to allocate space for memory block", GSL_ENOMEM);
    }

  m->pool->next = NULL;
  m->pool->free_slot = (unsigned char *) m->pool->block_ptr;

  return GSL_SUCCESS;
}

/* free memory pool */
static int
FUNCTION (spmatrix, pool_free) (TYPE (gsl_spmatrix) * m)
{
  while (m->pool != NULL)
    {
      gsl_spmatrix_pool * next = m->pool->next;
      free(m->pool->block_ptr);
      free(m->pool);
      m->pool = next;
    }

  return GSL_SUCCESS;
}

#if 0

static void *
FUNCTION (spmatrix, malloc) (size_t size, void * params)
{
  (void) params;
  return malloc (size);
}

static void
FUNCTION (spmatrix, free) (void * block, void * params)
{
  (void) params;
  free (block);
}

#else

static void *
FUNCTION (spmatrix, malloc) (size_t size, void * params)
{
  TYPE (gsl_spmatrix) * m = (TYPE (gsl_spmatrix) *) params;
  gsl_spmatrix_pool * pool = m->pool;
  void * ptr = (void *) pool->free_slot;
  pool->free_slot += size;
  return ptr;
}

static void
FUNCTION (spmatrix, free) (void * block, void * params)
{
  (void) block;
  (void) params;
}

#endif
