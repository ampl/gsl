/* avl.c
 * 
 * Copyright (C) 1998-2002, 2004 Free Software Foundation, Inc.
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

/* This code is originally from GNU libavl, with some modifications */

#include <config.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bst.h>
#include <gsl/gsl_errno.h>

typedef struct gsl_bst_avl_node avl_node;
typedef gsl_bst_avl_table avl_table;
typedef gsl_bst_avl_traverser avl_traverser;

#ifndef AVL_MAX_HEIGHT
#define AVL_MAX_HEIGHT GSL_BST_AVL_MAX_HEIGHT
#endif

/* Function types. */
typedef void avl_item_func (void *avl_item, void *avl_param);
typedef void *avl_copy_func (void *avl_item, void *avl_param);

/* tree functions */
static int avl_init(const gsl_bst_allocator * allocator,
                    gsl_bst_cmp_function * compare, void * params, void * vtable);
static size_t avl_nodes (const void * vtable);
static int avl_empty (void * vtable);
static void ** avl_probe (void * item, avl_table * table);
static void * avl_insert (void * item, void * vtable);
static void * avl_find (const void *item, const void * vtable);
static void * avl_remove (const void *item, void * vtable);

/* traverser functions */
static int avl_t_init (void * vtrav, const void * vtable);
static void * avl_t_first (void * vtrav, const void * vtable);
static void * avl_t_last (void * vtrav, const void * vtable);
static void * avl_t_find (const void * item, void * vtrav, const void * vtable);
static void * avl_t_insert (void * item, void * vtrav, void * vtable);
static void * avl_t_copy (void * vtrav, const void * vsrc);
static void * avl_t_next (void * vtrav);
static void * avl_t_prev (void * vtrav);
static void * avl_t_cur (const void * vtrav);
static void * avl_t_replace (void * vtrav, void * new_item);
static void avl_trav_refresh (avl_traverser * trav);

#if 0
static avl_table * avl_copy (const avl_table *, avl_copy_func *,
                             avl_item_func *);
static void *avl_replace (avl_table *, void *);
#endif

static int
avl_init(const gsl_bst_allocator * allocator,
         gsl_bst_cmp_function * compare, void * params, void * vtable)
{
  avl_table * table = (avl_table *) vtable;

  table->avl_alloc = allocator;
  table->avl_compare = compare;
  table->avl_param = params;

  table->avl_root = NULL;
  table->avl_count = 0;
  table->avl_generation = 0;

  return GSL_SUCCESS;
}

static size_t
avl_nodes (const void * vtable)
{
  const avl_table * table = (const avl_table *) vtable;
  return table->avl_count;
}

/* empty tree (delete all nodes) but do not free the tree itself */
static int
avl_empty (void * vtable)
{
  avl_table * table = (avl_table *) vtable;
  avl_node *p, *q;

  for (p = table->avl_root; p != NULL; p = q)
    {
      if (p->avl_link[0] == NULL)
        {
          q = p->avl_link[1];
          table->avl_alloc->free (p, table->avl_param);
        }
      else
        {
          q = p->avl_link[0];
          p->avl_link[0] = q->avl_link[1];
          q->avl_link[1] = p;
        }
    }

  table->avl_root = NULL;
  table->avl_count = 0;
  table->avl_generation = 0;

  return GSL_SUCCESS;
}

/*
avl_probe()
  Inserts |item| into |tree| and returns a pointer to |item|'s address.
If a duplicate item is found in the tree, returns a pointer to the existing
item without inserting |item|.
Returns |NULL| in case of memory allocation failure.
*/

static void **
avl_probe (void * item, avl_table * table)
{
  avl_node *y, *z; /* top node to update balance factor, and parent */
  avl_node *p, *q; /* iterator, and parent */
  avl_node *n;     /* newly inserted node */
  avl_node *w;     /* new root of rebalanced subtree */
  int dir;         /* direction to descend */

  unsigned char da[AVL_MAX_HEIGHT]; /* cached comparison results */
  int k = 0;                        /* number of cached results */

  z = (avl_node *) &table->avl_root;
  y = table->avl_root;
  dir = 0;
  for (q = z, p = y; p != NULL; q = p, p = p->avl_link[dir])
    {
      int cmp = table->avl_compare (item, p->avl_data, table->avl_param);

      if (cmp == 0)
        return &p->avl_data;

      if (p->avl_balance != 0)
        z = q, y = p, k = 0;

      da[k++] = dir = cmp > 0;
    }

  /* allocate a new node */
  n = q->avl_link[dir] = table->avl_alloc->alloc (sizeof *n, table->avl_param);
  if (n == NULL)
    return NULL;

  table->avl_count++;

  n->avl_data = item;
  n->avl_link[0] = n->avl_link[1] = NULL;
  n->avl_balance = 0;

  if (y == NULL)
    return &n->avl_data;

  for (p = y, k = 0; p != n; p = p->avl_link[da[k]], k++)
    if (da[k] == 0)
      p->avl_balance--;
    else
      p->avl_balance++;

  if (y->avl_balance == -2)
    {
      avl_node *x = y->avl_link[0];
      if (x->avl_balance == -1)
        {
          w = x;
          y->avl_link[0] = x->avl_link[1];
          x->avl_link[1] = y;
          x->avl_balance = y->avl_balance = 0;
        }
      else
        {
          w = x->avl_link[1];
          x->avl_link[1] = w->avl_link[0];
          w->avl_link[0] = x;
          y->avl_link[0] = w->avl_link[1];
          w->avl_link[1] = y;
          if (w->avl_balance == -1)
            x->avl_balance = 0, y->avl_balance = +1;
          else if (w->avl_balance == 0)
            x->avl_balance = y->avl_balance = 0;
          else /* |w->avl_balance == +1| */
            x->avl_balance = -1, y->avl_balance = 0;
          w->avl_balance = 0;
        }
    }
  else if (y->avl_balance == +2)
    {
      avl_node *x = y->avl_link[1];
      if (x->avl_balance == +1)
        {
          w = x;
          y->avl_link[1] = x->avl_link[0];
          x->avl_link[0] = y;
          x->avl_balance = y->avl_balance = 0;
        }
      else
        {
          w = x->avl_link[0];
          x->avl_link[0] = w->avl_link[1];
          w->avl_link[1] = x;
          y->avl_link[1] = w->avl_link[0];
          w->avl_link[0] = y;
          if (w->avl_balance == +1)
            x->avl_balance = 0, y->avl_balance = -1;
          else if (w->avl_balance == 0)
            x->avl_balance = y->avl_balance = 0;
          else /* |w->avl_balance == -1| */
            x->avl_balance = +1, y->avl_balance = 0;
          w->avl_balance = 0;
        }
    }
  else
    return &n->avl_data;

  z->avl_link[y != z->avl_link[0]] = w;
  table->avl_generation++;

  return &n->avl_data;
}

/*
avl_insert()
  Inserts |item| into |table|. Returns NULL if item successfully inserted,
or if a memory allocation error occurred. Otherwise, return duplicate
item.
*/

static void *
avl_insert (void * item, void * vtable)
{
  void **p = avl_probe (item, vtable);
  return p == NULL || *p == item ? NULL : *p;
}

/*
avl_find()
  Search for |item| in |table| and return a pointer to item if found.
Return NULL if not found.
*/

static void *
avl_find (const void * item, const void * vtable)
{
  const avl_table * table = (const avl_table *) vtable;
  avl_node *p;

  for (p = table->avl_root; p != NULL; )
    {
      int cmp = table->avl_compare (item, p->avl_data, table->avl_param);

      if (cmp < 0)
        p = p->avl_link[0];
      else if (cmp > 0)
        p = p->avl_link[1];
      else /* |cmp == 0| */
        return p->avl_data;
    }

  return NULL;
}

/*
avl_remove()
  Deletes from |table| and returns an item matching |item|.
Returns a null pointer if no matching item found.
*/

static void *
avl_remove (const void * item, void * vtable)
{
  avl_table * table = (avl_table *) vtable;

  /* stack of nodes */
  avl_node *pa[AVL_MAX_HEIGHT];      /* nodes */
  unsigned char da[AVL_MAX_HEIGHT];  /* |link[]| indexes */
  int k;                             /* stack pointer */

  avl_node *p;                       /* traverses tree to find node to delete */
  int cmp;                           /* result of comparison between |item| and |p| */

  k = 0;
  p = (avl_node *) &table->avl_root;
  for (cmp = -1; cmp != 0;
       cmp = table->avl_compare (item, p->avl_data, table->avl_param))
    {
      int dir = cmp > 0;

      pa[k] = p;
      da[k++] = dir;

      p = p->avl_link[dir];
      if (p == NULL)
        return NULL;
    }

  item = p->avl_data;

  if (p->avl_link[1] == NULL)
    pa[k - 1]->avl_link[da[k - 1]] = p->avl_link[0];
  else
    {
      avl_node *r = p->avl_link[1];
      if (r->avl_link[0] == NULL)
        {
          r->avl_link[0] = p->avl_link[0];
          r->avl_balance = p->avl_balance;
          pa[k - 1]->avl_link[da[k - 1]] = r;
          da[k] = 1;
          pa[k++] = r;
        }
      else
        {
          avl_node *s;
          int j = k++;

          for (;;)
            {
              da[k] = 0;
              pa[k++] = r;
              s = r->avl_link[0];
              if (s->avl_link[0] == NULL)
                break;

              r = s;
            }

          s->avl_link[0] = p->avl_link[0];
          r->avl_link[0] = s->avl_link[1];
          s->avl_link[1] = p->avl_link[1];
          s->avl_balance = p->avl_balance;

          pa[j - 1]->avl_link[da[j - 1]] = s;
          da[j] = 1;
          pa[j] = s;
        }
    }

  table->avl_alloc->free (p, table->avl_param);

  while (--k > 0)
    {
      avl_node *y = pa[k];

      if (da[k] == 0)
        {
          y->avl_balance++;
          if (y->avl_balance == +1)
            break;
          else if (y->avl_balance == +2)
            {
              avl_node *x = y->avl_link[1];
              if (x->avl_balance == -1)
                {
                  avl_node *w;
                  w = x->avl_link[0];
                  x->avl_link[0] = w->avl_link[1];
                  w->avl_link[1] = x;
                  y->avl_link[1] = w->avl_link[0];
                  w->avl_link[0] = y;
                  if (w->avl_balance == +1)
                    x->avl_balance = 0, y->avl_balance = -1;
                  else if (w->avl_balance == 0)
                    x->avl_balance = y->avl_balance = 0;
                  else /* |w->avl_balance == -1| */
                    x->avl_balance = +1, y->avl_balance = 0;
                  w->avl_balance = 0;
                  pa[k - 1]->avl_link[da[k - 1]] = w;
                }
              else
                {
                  y->avl_link[1] = x->avl_link[0];
                  x->avl_link[0] = y;
                  pa[k - 1]->avl_link[da[k - 1]] = x;
                  if (x->avl_balance == 0)
                    {
                      x->avl_balance = -1;
                      y->avl_balance = +1;
                      break;
                    }
                  else
                    x->avl_balance = y->avl_balance = 0;
                }
            }
        }
      else
        {
          y->avl_balance--;
          if (y->avl_balance == -1)
            break;
          else if (y->avl_balance == -2)
            {
              avl_node *x = y->avl_link[0];
              if (x->avl_balance == +1)
                {
                  avl_node *w;
                  w = x->avl_link[1];
                  x->avl_link[1] = w->avl_link[0];
                  w->avl_link[0] = x;
                  y->avl_link[0] = w->avl_link[1];
                  w->avl_link[1] = y;
                  if (w->avl_balance == -1)
                    x->avl_balance = 0, y->avl_balance = +1;
                  else if (w->avl_balance == 0)
                    x->avl_balance = y->avl_balance = 0;
                  else /* |w->avl_balance == +1| */
                    x->avl_balance = -1, y->avl_balance = 0;
                  w->avl_balance = 0;
                  pa[k - 1]->avl_link[da[k - 1]] = w;
                }
              else
                {
                  y->avl_link[0] = x->avl_link[1];
                  x->avl_link[1] = y;
                  pa[k - 1]->avl_link[da[k - 1]] = x;
                  if (x->avl_balance == 0)
                    {
                      x->avl_balance = +1;
                      y->avl_balance = -1;
                      break;
                    }
                  else
                    x->avl_balance = y->avl_balance = 0;
                }
            }
        }
    }

  table->avl_count--;
  table->avl_generation++;

  return (void *) item;
}

#if 0
/* Inserts |item| into |table|, replacing any duplicate item.
   Returns |NULL| if |item| was inserted without replacing a duplicate,
   or if a memory allocation error occurred.
   Otherwise, returns the item that was replaced. */
static void *
avl_replace (avl_table *table, void *item)
{
  void **p = avl_probe (table, item);
  if (p == NULL || *p == item)
    return NULL;
  else
    {
      void *r = *p;
      *p = item;
      return r;
    }
}


/* Destroys |new| with |avl_destroy (new, destroy)|,
   first setting right links of nodes in |stack| within |new|
   to null pointers to avoid touching uninitialized data. */
static void
copy_error_recovery (avl_node **stack, int height,
                     avl_table *new, avl_item_func *destroy)
{
  for (; height > 2; height -= 2)
    stack[height - 1]->avl_link[1] = NULL;
  avl_destroy (new, destroy);
}

/* Copies |org| to a newly created tree, which is returned.
   If |copy != NULL|, each data item in |org| is first passed to |copy|,
   and the return values are inserted into the tree,
   with |NULL| return values taken as indications of failure.
   On failure, destroys the partially created new tree,
   applying |destroy|, if non-null, to each item in the new tree so far,
   and returns |NULL|.
   If |allocator != NULL|, it is used for allocation in the new tree.
   Otherwise, the same allocator used for |org| is used. */
static avl_table *
avl_copy (const avl_table *org, avl_copy_func *copy,
          avl_item_func *destroy, struct libavl_allocator *allocator)
{
  avl_node *stack[2 * (AVL_MAX_HEIGHT + 1)];
  int height = 0;

  avl_table *new;
  const avl_node *x;
  avl_node *y;

  new = avl_alloc (org->avl_compare, org->avl_param,
                   allocator != NULL ? allocator : org->avl_alloc);
  if (new == NULL)
    return NULL;
  new->avl_count = org->avl_count;
  if (new->avl_count == 0)
    return new;

  x = (const avl_node *) &org->avl_root;
  y = (avl_node *) &new->avl_root;
  for (;;)
    {
      while (x->avl_link[0] != NULL)
        {
          y->avl_link[0] =
            new->allocator->alloc (sizeof *y->avl_link[0],
                                   new->avl_param);
          if (y->avl_link[0] == NULL)
            {
              if (y != (avl_node *) &new->avl_root)
                {
                  y->avl_data = NULL;
                  y->avl_link[1] = NULL;
                }

              copy_error_recovery (stack, height, new, destroy);
              return NULL;
            }

          stack[height++] = (avl_node *) x;
          stack[height++] = y;
          x = x->avl_link[0];
          y = y->avl_link[0];
        }
      y->avl_link[0] = NULL;

      for (;;)
        {
          y->avl_balance = x->avl_balance;
          if (copy == NULL)
            y->avl_data = x->avl_data;
          else
            {
              y->avl_data = copy (x->avl_data, org->avl_param);
              if (y->avl_data == NULL)
                {
                  y->avl_link[1] = NULL;
                  copy_error_recovery (stack, height, new, destroy);
                  return NULL;
                }
            }

          if (x->avl_link[1] != NULL)
            {
              y->avl_link[1] =
                new->allocator->alloc (sizeof *y->avl_link[1],
                                       new->avl_param);
              if (y->avl_link[1] == NULL)
                {
                  copy_error_recovery (stack, height, new, destroy);
                  return NULL;
                }

              x = x->avl_link[1];
              y = y->avl_link[1];
              break;
            }
          else
            y->avl_link[1] = NULL;

          if (height <= 2)
            return new;

          y = stack[--height];
          x = stack[--height];
        }
    }
}

#endif

/*
avl_t_init()
  Initializes |trav| for use with |tree| and selects the null node.
*/

static int
avl_t_init (void * vtrav, const void * vtable)
{
  avl_traverser * trav = (avl_traverser *) vtrav;
  const avl_table * table = (const avl_table *) vtable;

  trav->avl_table = table;
  trav->avl_node = NULL;
  trav->avl_height = 0;
  trav->avl_generation = table->avl_generation;

  return GSL_SUCCESS;
}

/*
avl_t_first()
  Initializes |trav| for |tree| and selects and returns a pointer to its least-valued item.
  Returns |NULL| if |tree| contains no nodes.
*/

static void *
avl_t_first (void * vtrav, const void * vtable)
{
  const avl_table * table = (const avl_table *) vtable;
  avl_traverser * trav = (avl_traverser *) vtrav;
  avl_node *x;

  trav->avl_table = table;
  trav->avl_height = 0;
  trav->avl_generation = table->avl_generation;

  x = table->avl_root;
  if (x != NULL)
    {
      while (x->avl_link[0] != NULL)
        {
          if (trav->avl_height >= AVL_MAX_HEIGHT)
            {
              GSL_ERROR_NULL("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->avl_stack[trav->avl_height++] = x;
          x = x->avl_link[0];
        }
    }

  trav->avl_node = x;

  return x != NULL ? x->avl_data : NULL;
}

/*
avl_t_last()
  Initializes |trav| for |tree| and selects and returns a pointer to its greatest-valued item.
  Returns |NULL| if |tree| contains no nodes.
*/

static void *
avl_t_last (void * vtrav, const void * vtable)
{
  const avl_table * table = (const avl_table *) vtable;
  avl_traverser * trav = (avl_traverser *) vtrav;
  avl_node *x;

  trav->avl_table = table;
  trav->avl_height = 0;
  trav->avl_generation = table->avl_generation;

  x = table->avl_root;
  if (x != NULL)
    {
      while (x->avl_link[1] != NULL)
        {
          if (trav->avl_height >= AVL_MAX_HEIGHT)
            {
              GSL_ERROR_NULL("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->avl_stack[trav->avl_height++] = x;
          x = x->avl_link[1];
        }
    }

  trav->avl_node = x;

  return x != NULL ? x->avl_data : NULL;
}

/*
avl_t_find()
  Searches for |item| in |table|. If found, initializes |trav| to
the item found and returns the item as well.  If there is no matching
item, initializes |trav| to the null item and returns |NULL|.
*/

static void *
avl_t_find (const void * item, void * vtrav, const void * vtable)
{
  const avl_table * table = (const avl_table *) vtable;
  avl_traverser * trav = (avl_traverser *) vtrav;
  avl_node *p, *q;

  trav->avl_table = table;
  trav->avl_height = 0;
  trav->avl_generation = table->avl_generation;

  for (p = table->avl_root; p != NULL; p = q)
    {
      int cmp = table->avl_compare (item, p->avl_data, table->avl_param);

      if (cmp < 0)
        q = p->avl_link[0];
      else if (cmp > 0)
        q = p->avl_link[1];
      else /* |cmp == 0| */
        {
          trav->avl_node = p;
          return p->avl_data;
        }

      if (trav->avl_height >= AVL_MAX_HEIGHT)
        {
          GSL_ERROR_NULL("traverser height exceeds maximum", GSL_ETABLE);
        }

      trav->avl_stack[trav->avl_height++] = p;
    }

  trav->avl_height = 0;
  trav->avl_node = NULL;

  return NULL;
}

/*
avl_t_insert()
  Attempts to insert |item| into |table|.  If |item| is inserted
successfully, it is returned and |trav| is initialized to its location.
If a duplicate is found, it is returned and |trav| is initialized to
its location.  No replacement of the item occurs.
If a memory allocation failure occurs, |NULL| is returned and |trav|
is initialized to the null item.
*/

static void *
avl_t_insert (void * item, void * vtrav, void * vtable)
{
  avl_table * table = (avl_table *) vtable;
  avl_traverser * trav = (avl_traverser *) vtrav;
  void **p;

  p = avl_probe (item, table);
  if (p != NULL)
    {
      trav->avl_table = table;
      trav->avl_node = ((avl_node *) ((char *) p - offsetof (avl_node, avl_data)));
      trav->avl_generation = table->avl_generation - 1;
      return *p;
    }
  else
    {
      avl_t_init (vtrav, vtable);
      return NULL;
    }
}

/*
avl_t_copy()
  Initializes |trav| to have the same current node as |src|.
*/

static void *
avl_t_copy (void * vtrav, const void * vsrc)
{
  const avl_traverser * src = (const avl_traverser *) vsrc;
  avl_traverser * trav = (avl_traverser *) vtrav;

  if (trav != src)
    {
      trav->avl_table = src->avl_table;
      trav->avl_node = src->avl_node;
      trav->avl_generation = src->avl_generation;
      if (trav->avl_generation == trav->avl_table->avl_generation)
        {
          trav->avl_height = src->avl_height;
          memcpy (trav->avl_stack, (const void *) src->avl_stack,
                  sizeof *trav->avl_stack * trav->avl_height);
        }
    }

  return trav->avl_node != NULL ? trav->avl_node->avl_data : NULL;
}

/*
avl_t_next()
  Returns the next data item in in-order within the tree
being traversed with |trav|, or if there are no more data items returns NULL.
*/

static void *
avl_t_next (void * vtrav)
{
  avl_traverser * trav = (avl_traverser *) vtrav;
  avl_node *x;

  if (trav->avl_generation != trav->avl_table->avl_generation)
    avl_trav_refresh (trav);

  x = trav->avl_node;
  if (x == NULL)
    {
      return avl_t_first (vtrav, trav->avl_table);
    }
  else if (x->avl_link[1] != NULL)
    {
      if (trav->avl_height >= AVL_MAX_HEIGHT)
        {
          GSL_ERROR_NULL("traverser height exceeds maximum", GSL_ETABLE);
        }

      trav->avl_stack[trav->avl_height++] = x;
      x = x->avl_link[1];

      while (x->avl_link[0] != NULL)
        {
          if (trav->avl_height >= AVL_MAX_HEIGHT)
            {
              GSL_ERROR_NULL("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->avl_stack[trav->avl_height++] = x;
          x = x->avl_link[0];
        }
    }
  else
    {
      avl_node *y;

      do
        {
          if (trav->avl_height == 0)
            {
              trav->avl_node = NULL;
              return NULL;
            }

          y = x;
          x = trav->avl_stack[--trav->avl_height];
        }
      while (y == x->avl_link[1]);
    }

  trav->avl_node = x;

  return x->avl_data;
}

/*
avl_t_prev()
  Returns the previous data item in inorder within the tree being
traversed with |trav|, or if there are no more data items returns NULL.
*/

static void *
avl_t_prev (void * vtrav)
{
  avl_traverser * trav = (avl_traverser *) vtrav;
  avl_node *x;

  if (trav->avl_generation != trav->avl_table->avl_generation)
    avl_trav_refresh (trav);

  x = trav->avl_node;
  if (x == NULL)
    {
      return avl_t_last (vtrav, trav->avl_table);
    }
  else if (x->avl_link[0] != NULL)
    {
      if (trav->avl_height >= AVL_MAX_HEIGHT)
        {
          GSL_ERROR_NULL("traverser height exceeds maximum", GSL_ETABLE);
        }

      trav->avl_stack[trav->avl_height++] = x;
      x = x->avl_link[0];

      while (x->avl_link[1] != NULL)
        {
          if (trav->avl_height >= AVL_MAX_HEIGHT)
            {
              GSL_ERROR_NULL("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->avl_stack[trav->avl_height++] = x;
          x = x->avl_link[1];
        }
    }
  else
    {
      avl_node *y;

      do
        {
          if (trav->avl_height == 0)
            {
              trav->avl_node = NULL;
              return NULL;
            }

          y = x;
          x = trav->avl_stack[--trav->avl_height];
        }
      while (y == x->avl_link[0]);
    }

  trav->avl_node = x;

  return x->avl_data;
}

static void *
avl_t_cur (const void * vtrav)
{
  const avl_traverser * trav = (const avl_traverser *) vtrav;
  return trav->avl_node != NULL ? trav->avl_node->avl_data : NULL;
}

/*
avl_t_replace()
  Replaces the current item in |trav| by |new| and returns the item replaced.
|trav| must not have the null item selected. The new item must not upset the
ordering of the tree.
*/

static void *
avl_t_replace (void * vtrav, void * new_item)
{
  avl_traverser * trav = (avl_traverser *) vtrav;
  void *old;

  old = trav->avl_node->avl_data;
  trav->avl_node->avl_data = new_item;

  return old;
}

/*
avl_trav_refresh()
  Refreshes the stack of parent pointers in |trav| and updates its generation number
*/

static void
avl_trav_refresh (avl_traverser * trav)
{
  trav->avl_generation = trav->avl_table->avl_generation;

  if (trav->avl_node != NULL)
    {
      gsl_bst_cmp_function *cmp = trav->avl_table->avl_compare;
      void *param = trav->avl_table->avl_param;
      avl_node *node = trav->avl_node;
      avl_node *i;

      trav->avl_height = 0;
      for (i = trav->avl_table->avl_root; i != node; )
        {
          if (trav->avl_height >= AVL_MAX_HEIGHT)
            {
              GSL_ERROR_VOID("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->avl_stack[trav->avl_height++] = i;
          i = i->avl_link[cmp (node->avl_data, i->avl_data, param) > 0];
        }
    }
}

static const gsl_bst_type avl_tree_type =
{
  "AVL",
  sizeof(avl_node),
  avl_init,
  avl_nodes,
  avl_insert,
  avl_find,
  avl_remove,
  avl_empty,

  avl_t_init,
  avl_t_first,
  avl_t_last,
  avl_t_find,
  avl_t_insert,
  avl_t_copy,
  avl_t_next,
  avl_t_prev,
  avl_t_cur,
  avl_t_replace
};

const gsl_bst_type * gsl_bst_avl = &avl_tree_type;
