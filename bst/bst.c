/* bst.c
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

#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bst.h>
#include <gsl/gsl_errno.h>

static void * bst_malloc(size_t size, void * params);
static void bst_free(void * block, void * params);

static const gsl_bst_allocator bst_default_allocator =
{
  bst_malloc,
  bst_free
};

/*
gsl_bst_alloc()
  Allocate binary search tree

Inputs: T         - tree type
        allocator - memory allocator
        compare   - comparison function
        params    - parameters to pass to allocator and compare
*/

gsl_bst_workspace *
gsl_bst_alloc(const gsl_bst_type * T, const gsl_bst_allocator * allocator,
              gsl_bst_cmp_function * compare, void * params)
{
  int status;
  gsl_bst_workspace *w;

  w = calloc(1, sizeof(gsl_bst_workspace));
  if (w == NULL)
    {
      GSL_ERROR_NULL("failed to allocate bst workspace", GSL_ENOMEM);
    }

  w->type = T;

  status = (w->type->init)(allocator != NULL ? allocator : &bst_default_allocator, compare, params, (void *) &w->table);
  if (status)
    {
      gsl_bst_free(w);
      GSL_ERROR_NULL("failed to initialize bst", GSL_EFAILED);
    }

  return w;
}

void
gsl_bst_free(gsl_bst_workspace * w)
{
  /* free tree nodes */
  gsl_bst_empty(w);

  free(w);
}

/* delete all nodes from tree */
int
gsl_bst_empty(gsl_bst_workspace * w)
{
  return (w->type->empty)((void *) &w->table);
}

/*
gsl_bst_insert()
  Inserts |item| into tree
  
AVL: if duplicate found, returns pointer to item without inserting
AVLmult: if duplicate found, increase multiplicity for that node

If no duplicate found, insert item and return pointer to item.
Returns NULL if a memory allocation error occurred.
*/

void *
gsl_bst_insert(void * item, gsl_bst_workspace * w)
{
  return (w->type->insert)(item, (void *) &w->table);
}

void *
gsl_bst_find(const void * item, const gsl_bst_workspace * w)
{
  return (w->type->find)(item, (const void *) &w->table);
}

void *
gsl_bst_remove(const void * item, gsl_bst_workspace * w)
{
  return (w->type->remove)(item, (void *) &w->table);
}

/* return number of nodes in tree */
size_t
gsl_bst_nodes(const gsl_bst_workspace * w)
{
  return (w->type->nodes)((const void *) &w->table);
}

/* return size (in bytes) of each node in tree */
size_t
gsl_bst_node_size(const gsl_bst_workspace * w)
{
  return w->type->node_size;
}

const char *
gsl_bst_name(const gsl_bst_workspace * w)
{
  return w->type->name;
}

/**********************************************
 * INTERNAL ROUTINES                          *
 **********************************************/

static void *
bst_malloc(size_t size, void * params)
{
  (void) params; /* avoid unused parameter warning */
  return malloc(size);
}

static void
bst_free(void * block, void * params)
{
  (void) params; /* avoid unused parameter warning */
  free(block);
}
