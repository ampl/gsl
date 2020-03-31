/* rb.c
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

typedef struct gsl_bst_rb_node rb_node;
typedef gsl_bst_rb_table rb_table;
typedef gsl_bst_rb_traverser rb_traverser;

enum rb_color
{
  RB_BLACK, /* black */
  RB_RED    /* red */
};

#ifndef RB_MAX_HEIGHT
#define RB_MAX_HEIGHT GSL_BST_RB_MAX_HEIGHT
#endif

/* tree functions */
static int rb_init(const gsl_bst_allocator * allocator, gsl_bst_cmp_function * compare,
                   void * params, void * vtable);
static size_t rb_nodes (const void * vtable);
static int rb_empty (void * vtable);
static void ** rb_probe (void * item, rb_table * table);
static void * rb_insert (void * item, void * vtable);
static void * rb_find (const void * item, const void * vtable);
static void * rb_remove (const void * item, void * vtable);

/* traverser functions */
static int rb_t_init (void * vtrav, const void * vtable);
static void * rb_t_first (void * vtrav, const void * vtable);
static void * rb_t_last (void * vtrav, const void * vtable);
static void * rb_t_find (const void * item, void * vtrav, const void * vtable);
static void * rb_t_insert (void * item, void * vtrav, void * vtable);
static void * rb_t_copy (void * vtrav, const void * vsrc);
static void * rb_t_next (void * vtrav);
static void * rb_t_prev (void * vtrav);
static void * rb_t_cur (const void * vtrav);
static void * rb_t_replace (void * vtrav, void * new_item);
static void rb_trav_refresh (rb_traverser *trav);

static int
rb_init(const gsl_bst_allocator * allocator, gsl_bst_cmp_function * compare,
        void * params, void * vtable)
{
  rb_table * table = (rb_table *) vtable;

  table->rb_alloc = allocator;
  table->rb_compare = compare;
  table->rb_param = params;

  table->rb_root = NULL;
  table->rb_count = 0;
  table->rb_generation = 0;

  return GSL_SUCCESS;
}

static size_t
rb_nodes (const void * vtable)
{
  const rb_table * table = (const rb_table *) vtable;
  return table->rb_count;
}

/* empty tree (delete all nodes) but do not free the tree itself */
static int
rb_empty (void * vtable)
{
  rb_table * table = (rb_table *) vtable;
  rb_node *p, *q;

  for (p = table->rb_root; p != NULL; p = q)
    {
      if (p->rb_link[0] == NULL)
        {
          q = p->rb_link[1];
          table->rb_alloc->free (p, table->rb_param);
        }
      else
        {
          q = p->rb_link[0];
          p->rb_link[0] = q->rb_link[1];
          q->rb_link[1] = p;
        }
    }

  table->rb_root = NULL;
  table->rb_count = 0;
  table->rb_generation = 0;

  return GSL_SUCCESS;
}

/* Inserts |item| into |table| and returns a pointer to |item|'s address.
   If a duplicate item is found in the tree,
   returns a pointer to the duplicate without inserting |item|.
   Returns |NULL| in case of memory allocation failure. */

static void **
rb_probe (void * item, rb_table * table)
{
  rb_node *pa[RB_MAX_HEIGHT];      /* nodes on stack */
  unsigned char da[RB_MAX_HEIGHT]; /* directions moved from stack nodes */
  int k;                           /* stack height */

  rb_node *p; /* traverses tree looking for insertion point */
  rb_node *n; /* newly inserted node */

  pa[0] = (rb_node *) &table->rb_root;
  da[0] = 0;
  k = 1;
  for (p = table->rb_root; p != NULL; p = p->rb_link[da[k - 1]])
    {
      int cmp = table->rb_compare (item, p->rb_data, table->rb_param);
      if (cmp == 0)
        return &p->rb_data;

      pa[k] = p;
      da[k++] = cmp > 0;
    }

  n = pa[k - 1]->rb_link[da[k - 1]] =
    table->rb_alloc->alloc (sizeof *n, table->rb_param);
  if (n == NULL)
    return NULL;

  n->rb_data = item;
  n->rb_link[0] = n->rb_link[1] = NULL;
  n->rb_color = RB_RED;
  table->rb_count++;
  table->rb_generation++;

  while (k >= 3 && pa[k - 1]->rb_color == RB_RED)
    {
      if (da[k - 2] == 0)
        {
          rb_node *y = pa[k - 2]->rb_link[1];
          if (y != NULL && y->rb_color == RB_RED)
            {
              pa[k - 1]->rb_color = y->rb_color = RB_BLACK;
              pa[k - 2]->rb_color = RB_RED;
              k -= 2;
            }
          else
            {
              rb_node *x;

              if (da[k - 1] == 0)
                y = pa[k - 1];
              else
                {
                  x = pa[k - 1];
                  y = x->rb_link[1];
                  x->rb_link[1] = y->rb_link[0];
                  y->rb_link[0] = x;
                  pa[k - 2]->rb_link[0] = y;
                }

              x = pa[k - 2];
              x->rb_color = RB_RED;
              y->rb_color = RB_BLACK;

              x->rb_link[0] = y->rb_link[1];
              y->rb_link[1] = x;
              pa[k - 3]->rb_link[da[k - 3]] = y;
              break;
            }
        }
      else
        {
          rb_node *y = pa[k - 2]->rb_link[0];
          if (y != NULL && y->rb_color == RB_RED)
            {
              pa[k - 1]->rb_color = y->rb_color = RB_BLACK;
              pa[k - 2]->rb_color = RB_RED;
              k -= 2;
            }
          else
            {
              rb_node *x;

              if (da[k - 1] == 1)
                y = pa[k - 1];
              else
                {
                  x = pa[k - 1];
                  y = x->rb_link[0];
                  x->rb_link[0] = y->rb_link[1];
                  y->rb_link[1] = x;
                  pa[k - 2]->rb_link[1] = y;
                }

              x = pa[k - 2];
              x->rb_color = RB_RED;
              y->rb_color = RB_BLACK;

              x->rb_link[1] = y->rb_link[0];
              y->rb_link[0] = x;
              pa[k - 3]->rb_link[da[k - 3]] = y;
              break;
            }
        }
    }

  table->rb_root->rb_color = RB_BLACK;

  return &n->rb_data;
}

/* Inserts |item| into |table|.
   Returns |NULL| if |item| was successfully inserted
   or if a memory allocation error occurred.
   Otherwise, returns the duplicate item. */

static void *
rb_insert (void * item, void * vtable)
{
  void **p = rb_probe (item, vtable);
  return p == NULL || *p == item ? NULL : *p;
}

/* Search |table| for an item matching |item|, and return it if found.
   Otherwise return |NULL|. */
static void *
rb_find (const void * item, const void * vtable)
{
  const rb_table * table = (const rb_table *) vtable;
  const rb_node *p;

  for (p = table->rb_root; p != NULL; )
    {
      int cmp = table->rb_compare (item, p->rb_data, table->rb_param);

      if (cmp < 0)
        p = p->rb_link[0];
      else if (cmp > 0)
        p = p->rb_link[1];
      else /* |cmp == 0| */
        return p->rb_data;
    }

  return NULL;
}

#if 0 /*XXX*/

/* Inserts |item| into |table|, replacing any duplicate item.
   Returns |NULL| if |item| was inserted without replacing a duplicate,
   or if a memory allocation error occurred.
   Otherwise, returns the item that was replaced. */
void *
rb_replace (struct rb_table *table, void *item)
{
  void **p = rb_probe (table, item);
  if (p == NULL || *p == item)
    return NULL;
  else
    {
      void *r = *p;
      *p = item;
      return r;
    }
}

#endif

/* Deletes from |table| and returns an item matching |item|.
   Returns a null pointer if no matching item found. */

static void *
rb_remove (const void * item, void * vtable)
{
  rb_table * table = (rb_table *) vtable;
  rb_node *pa[RB_MAX_HEIGHT];      /* nodes on stack */
  unsigned char da[RB_MAX_HEIGHT]; /* directions moved from stack nodes */
  int k;                           /* stack height */

  rb_node *p;    /* the node to delete, or a node part way to it */
  int cmp;       /* result of comparison between |item| and |p| */

  k = 0;
  p = (rb_node *) &table->rb_root;
  for (cmp = -1; cmp != 0;
       cmp = table->rb_compare (item, p->rb_data, table->rb_param))
    {
      int dir = cmp > 0;

      pa[k] = p;
      da[k++] = dir;

      p = p->rb_link[dir];
      if (p == NULL)
        return NULL;
    }
  item = p->rb_data;

  if (p->rb_link[1] == NULL)
    pa[k - 1]->rb_link[da[k - 1]] = p->rb_link[0];
  else
    {
      enum rb_color t;
      rb_node *r = p->rb_link[1];

      if (r->rb_link[0] == NULL)
        {
          r->rb_link[0] = p->rb_link[0];
          t = r->rb_color;
          r->rb_color = p->rb_color;
          p->rb_color = t;
          pa[k - 1]->rb_link[da[k - 1]] = r;
          da[k] = 1;
          pa[k++] = r;
        }
      else
        {
          rb_node *s;
          int j = k++;

          for (;;)
            {
              da[k] = 0;
              pa[k++] = r;
              s = r->rb_link[0];
              if (s->rb_link[0] == NULL)
                break;

              r = s;
            }

          da[j] = 1;
          pa[j] = s;
          pa[j - 1]->rb_link[da[j - 1]] = s;

          s->rb_link[0] = p->rb_link[0];
          r->rb_link[0] = s->rb_link[1];
          s->rb_link[1] = p->rb_link[1];

          t = s->rb_color;
          s->rb_color = p->rb_color;
          p->rb_color = t;
        }
    }

  if (p->rb_color == RB_BLACK)
    {
      for (;;)
        {
          rb_node *x = pa[k - 1]->rb_link[da[k - 1]];
          if (x != NULL && x->rb_color == RB_RED)
            {
              x->rb_color = RB_BLACK;
              break;
            }
          if (k < 2)
            break;

          if (da[k - 1] == 0)
            {
              rb_node *w = pa[k - 1]->rb_link[1];

              if (w->rb_color == RB_RED)
                {
                  w->rb_color = RB_BLACK;
                  pa[k - 1]->rb_color = RB_RED;

                  pa[k - 1]->rb_link[1] = w->rb_link[0];
                  w->rb_link[0] = pa[k - 1];
                  pa[k - 2]->rb_link[da[k - 2]] = w;

                  pa[k] = pa[k - 1];
                  da[k] = 0;
                  pa[k - 1] = w;
                  k++;

                  w = pa[k - 1]->rb_link[1];
                }

              if ((w->rb_link[0] == NULL
                   || w->rb_link[0]->rb_color == RB_BLACK)
                  && (w->rb_link[1] == NULL
                      || w->rb_link[1]->rb_color == RB_BLACK))
                w->rb_color = RB_RED;
              else
                {
                  if (w->rb_link[1] == NULL
                      || w->rb_link[1]->rb_color == RB_BLACK)
                    {
                      rb_node *y = w->rb_link[0];
                      y->rb_color = RB_BLACK;
                      w->rb_color = RB_RED;
                      w->rb_link[0] = y->rb_link[1];
                      y->rb_link[1] = w;
                      w = pa[k - 1]->rb_link[1] = y;
                    }

                  w->rb_color = pa[k - 1]->rb_color;
                  pa[k - 1]->rb_color = RB_BLACK;
                  w->rb_link[1]->rb_color = RB_BLACK;

                  pa[k - 1]->rb_link[1] = w->rb_link[0];
                  w->rb_link[0] = pa[k - 1];
                  pa[k - 2]->rb_link[da[k - 2]] = w;
                  break;
                }
            }
          else
            {
              rb_node *w = pa[k - 1]->rb_link[0];

              if (w->rb_color == RB_RED)
                {
                  w->rb_color = RB_BLACK;
                  pa[k - 1]->rb_color = RB_RED;

                  pa[k - 1]->rb_link[0] = w->rb_link[1];
                  w->rb_link[1] = pa[k - 1];
                  pa[k - 2]->rb_link[da[k - 2]] = w;

                  pa[k] = pa[k - 1];
                  da[k] = 1;
                  pa[k - 1] = w;
                  k++;

                  w = pa[k - 1]->rb_link[0];
                }

              if ((w->rb_link[0] == NULL
                   || w->rb_link[0]->rb_color == RB_BLACK)
                  && (w->rb_link[1] == NULL
                      || w->rb_link[1]->rb_color == RB_BLACK))
                w->rb_color = RB_RED;
              else
                {
                  if (w->rb_link[0] == NULL
                      || w->rb_link[0]->rb_color == RB_BLACK)
                    {
                      rb_node *y = w->rb_link[1];
                      y->rb_color = RB_BLACK;
                      w->rb_color = RB_RED;
                      w->rb_link[1] = y->rb_link[0];
                      y->rb_link[0] = w;
                      w = pa[k - 1]->rb_link[0] = y;
                    }

                  w->rb_color = pa[k - 1]->rb_color;
                  pa[k - 1]->rb_color = RB_BLACK;
                  w->rb_link[0]->rb_color = RB_BLACK;

                  pa[k - 1]->rb_link[0] = w->rb_link[1];
                  w->rb_link[1] = pa[k - 1];
                  pa[k - 2]->rb_link[da[k - 2]] = w;
                  break;
                }
            }

          k--;
        }

    }

  table->rb_alloc->free (p, table->rb_param);
  table->rb_count--;
  table->rb_generation++;

  return (void *) item;
}

/* Initializes |trav| for use with |tree| and selects the null node. */
static int
rb_t_init (void * vtrav, const void * vtable)
{
  rb_traverser * trav = (rb_traverser *) vtrav;
  const rb_table * table = (const rb_table *) vtable;

  trav->rb_table = table;
  trav->rb_node = NULL;
  trav->rb_height = 0;
  trav->rb_generation = table->rb_generation;

  return GSL_SUCCESS;
}

/* Initializes |trav| for |table|
   and selects and returns a pointer to its least-valued item.
   Returns |NULL| if |table| contains no nodes. */
static void *
rb_t_first (void * vtrav, const void * vtable)
{
  const rb_table * table = (const rb_table *) vtable;
  rb_traverser * trav = (rb_traverser *) vtrav;
  rb_node *x;

  trav->rb_table = table;
  trav->rb_height = 0;
  trav->rb_generation = table->rb_generation;

  x = table->rb_root;
  if (x != NULL)
    {
      while (x->rb_link[0] != NULL)
        {
          if (trav->rb_height >= RB_MAX_HEIGHT)
            {
              GSL_ERROR_NULL ("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->rb_stack[trav->rb_height++] = x;
          x = x->rb_link[0];
        }
    }

  trav->rb_node = x;

  return x != NULL ? x->rb_data : NULL;
}

/* Initializes |trav| for |table|
   and selects and returns a pointer to its greatest-valued item.
   Returns |NULL| if |table| contains no nodes. */
static void *
rb_t_last (void * vtrav, const void * vtable)
{
  const rb_table * table = (const rb_table *) vtable;
  rb_traverser * trav = (rb_traverser *) vtrav;
  rb_node *x;

  trav->rb_table = table;
  trav->rb_height = 0;
  trav->rb_generation = table->rb_generation;

  x = table->rb_root;
  if (x != NULL)
    {
      while (x->rb_link[1] != NULL)
        {
          if (trav->rb_height >= RB_MAX_HEIGHT)
            {
              GSL_ERROR_NULL ("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->rb_stack[trav->rb_height++] = x;
          x = x->rb_link[1];
        }
    }

  trav->rb_node = x;

  return x != NULL ? x->rb_data : NULL;
}

/* Searches for |item| in |table|.
   If found, initializes |trav| to the item found and returns the item
   as well.
   If there is no matching item, initializes |trav| to the null item
   and returns |NULL|. */
static void *
rb_t_find (const void * item, void * vtrav, const void * vtable)
{
  const rb_table * table = (const rb_table *) vtable;
  rb_traverser * trav = (rb_traverser *) vtrav;
  rb_node *p, *q;

  trav->rb_table = table;
  trav->rb_height = 0;
  trav->rb_generation = table->rb_generation;
  for (p = table->rb_root; p != NULL; p = q)
    {
      int cmp = table->rb_compare (item, p->rb_data, table->rb_param);

      if (cmp < 0)
        q = p->rb_link[0];
      else if (cmp > 0)
        q = p->rb_link[1];
      else /* |cmp == 0| */
        {
          trav->rb_node = p;
          return p->rb_data;
        }

      if (trav->rb_height >= RB_MAX_HEIGHT)
        {
          GSL_ERROR_NULL ("traverser height exceeds maximum", GSL_ETABLE);
        }

      trav->rb_stack[trav->rb_height++] = p;
    }

  trav->rb_height = 0;
  trav->rb_node = NULL;

  return NULL;
}

/* Attempts to insert |item| into |table|.
   If |item| is inserted successfully, it is returned and |trav| is
   initialized to its location.
   If a duplicate is found, it is returned and |trav| is initialized to
   its location.  No replacement of the item occurs.
   If a memory allocation failure occurs, |NULL| is returned and |trav|
   is initialized to the null item. */
static void *
rb_t_insert (void * item, void * vtrav, void * vtable)
{
  rb_table * table = (rb_table *) vtable;
  rb_traverser * trav = (rb_traverser *) vtrav;
  void **p;

  p = rb_probe (item, table);
  if (p != NULL)
    {
      trav->rb_table = table;
      trav->rb_node = ((rb_node *) ((char *) p - offsetof (rb_node, rb_data)));
      trav->rb_generation = table->rb_generation - 1;
      return *p;
    }
  else
    {
      rb_t_init (vtrav, vtable);
      return NULL;
    }
}

/* Initializes |trav| to have the same current node as |src|. */
static void *
rb_t_copy (void * vtrav, const void * vsrc)
{
  const rb_traverser * src = (const rb_traverser *) vsrc;
  rb_traverser * trav = (rb_traverser *) vtrav;

  if (trav != src)
    {
      trav->rb_table = src->rb_table;
      trav->rb_node = src->rb_node;
      trav->rb_generation = src->rb_generation;
      if (trav->rb_generation == trav->rb_table->rb_generation)
        {
          trav->rb_height = src->rb_height;
          memcpy (trav->rb_stack, (const void *) src->rb_stack,
                  sizeof *trav->rb_stack * trav->rb_height);
        }
    }

  return trav->rb_node != NULL ? trav->rb_node->rb_data : NULL;
}

/* Returns the next data item in inorder
   within the tree being traversed with |trav|,
   or if there are no more data items returns |NULL|. */
static void *
rb_t_next (void * vtrav)
{
  rb_traverser * trav = (rb_traverser *) vtrav;
  rb_node *x;

  if (trav->rb_generation != trav->rb_table->rb_generation)
    rb_trav_refresh (trav);

  x = trav->rb_node;
  if (x == NULL)
    {
      return rb_t_first (vtrav, trav->rb_table);
    }
  else if (x->rb_link[1] != NULL)
    {
      if (trav->rb_height >= RB_MAX_HEIGHT)
        {
          GSL_ERROR_NULL ("traverser height exceeds maximum", GSL_ETABLE);
        }

      trav->rb_stack[trav->rb_height++] = x;
      x = x->rb_link[1];

      while (x->rb_link[0] != NULL)
        {
          if (trav->rb_height >= RB_MAX_HEIGHT)
            {
              GSL_ERROR_NULL ("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->rb_stack[trav->rb_height++] = x;
          x = x->rb_link[0];
        }
    }
  else
    {
      rb_node *y;

      do
        {
          if (trav->rb_height == 0)
            {
              trav->rb_node = NULL;
              return NULL;
            }

          y = x;
          x = trav->rb_stack[--trav->rb_height];
        }
      while (y == x->rb_link[1]);
    }

  trav->rb_node = x;

  return x->rb_data;
}

/* Returns the previous data item in inorder
   within the tree being traversed with |trav|,
   or if there are no more data items returns |NULL|. */
static void *
rb_t_prev (void * vtrav)
{
  rb_traverser * trav = (rb_traverser *) vtrav;
  rb_node *x;

  if (trav->rb_generation != trav->rb_table->rb_generation)
    rb_trav_refresh (trav);

  x = trav->rb_node;
  if (x == NULL)
    {
      return rb_t_last (vtrav, trav->rb_table);
    }
  else if (x->rb_link[0] != NULL)
    {
      if (trav->rb_height >= RB_MAX_HEIGHT)
        {
          GSL_ERROR_NULL ("traverser height exceeds maximum", GSL_ETABLE);
        }

      trav->rb_stack[trav->rb_height++] = x;
      x = x->rb_link[0];

      while (x->rb_link[1] != NULL)
        {
          if (trav->rb_height >= RB_MAX_HEIGHT)
            {
              GSL_ERROR_NULL ("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->rb_stack[trav->rb_height++] = x;
          x = x->rb_link[1];
        }
    }
  else
    {
      rb_node *y;

      do
        {
          if (trav->rb_height == 0)
            {
              trav->rb_node = NULL;
              return NULL;
            }

          y = x;
          x = trav->rb_stack[--trav->rb_height];
        }
      while (y == x->rb_link[0]);
    }

  trav->rb_node = x;

  return x->rb_data;
}

/* Returns |trav|'s current item. */
static void *
rb_t_cur (const void * vtrav)
{
  const rb_traverser * trav = (const rb_traverser *) vtrav;
  return trav->rb_node != NULL ? trav->rb_node->rb_data : NULL;
}

/* Replaces the current item in |trav| by |new| and returns the item replaced.
   |trav| must not have the null item selected.
   The new item must not upset the ordering of the tree. */
static void *
rb_t_replace (void * vtrav, void * new_item)
{
  rb_traverser * trav = (rb_traverser *) vtrav;
  void *old;

  old = trav->rb_node->rb_data;
  trav->rb_node->rb_data = new_item;

  return old;
}

#if 0 /*XXX*/
/* Destroys |new| with |rb_destroy (new, destroy)|,
   first setting right links of nodes in |stack| within |new|
   to null pointers to avoid touching uninitialized data. */
static void
copy_error_recovery (rb_node **stack, int height,
                     struct rb_table *new, rb_item_func *destroy)
{
  assert (stack != NULL && height >= 0 && new != NULL);

  for (; height > 2; height -= 2)
    stack[height - 1]->rb_link[1] = NULL;
  rb_destroy (new, destroy);
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
struct rb_table *
rb_copy (const struct rb_table *org, rb_copy_func *copy,
          rb_item_func *destroy, struct libavl_allocator *allocator)
{
  rb_node *stack[2 * (RB_MAX_HEIGHT + 1)];
  int height = 0;

  struct rb_table *new;
  const rb_node *x;
  rb_node *y;

  assert (org != NULL);
  new = rb_create (org->rb_compare, org->rb_param,
                    allocator != NULL ? allocator : org->rb_alloc);
  if (new == NULL)
    return NULL;
  new->rb_count = org->rb_count;
  if (new->rb_count == 0)
    return new;

  x = (const rb_node *) &org->rb_root;
  y = (rb_node *) &new->rb_root;
  for (;;)
    {
      while (x->rb_link[0] != NULL)
        {
          assert (height < 2 * (RB_MAX_HEIGHT + 1));

          y->rb_link[0] =
            new->rb_alloc->libavl_malloc (new->rb_alloc,
                                           sizeof *y->rb_link[0]);
          if (y->rb_link[0] == NULL)
            {
              if (y != (rb_node *) &new->rb_root)
                {
                  y->rb_data = NULL;
                  y->rb_link[1] = NULL;
                }

              copy_error_recovery (stack, height, new, destroy);
              return NULL;
            }

          stack[height++] = (rb_node *) x;
          stack[height++] = y;
          x = x->rb_link[0];
          y = y->rb_link[0];
        }
      y->rb_link[0] = NULL;

      for (;;)
        {
          y->rb_color = x->rb_color;
          if (copy == NULL)
            y->rb_data = x->rb_data;
          else
            {
              y->rb_data = copy (x->rb_data, org->rb_param);
              if (y->rb_data == NULL)
                {
                  y->rb_link[1] = NULL;
                  copy_error_recovery (stack, height, new, destroy);
                  return NULL;
                }
            }

          if (x->rb_link[1] != NULL)
            {
              y->rb_link[1] =
                new->rb_alloc->libavl_malloc (new->rb_alloc,
                                               sizeof *y->rb_link[1]);
              if (y->rb_link[1] == NULL)
                {
                  copy_error_recovery (stack, height, new, destroy);
                  return NULL;
                }

              x = x->rb_link[1];
              y = y->rb_link[1];
              break;
            }
          else
            y->rb_link[1] = NULL;

          if (height <= 2)
            return new;

          y = stack[--height];
          x = stack[--height];
        }
    }
}

#endif

/* Refreshes the stack of parent pointers in |trav|
   and updates its generation number. */
static void
rb_trav_refresh (rb_traverser *trav)
{
  trav->rb_generation = trav->rb_table->rb_generation;

  if (trav->rb_node != NULL)
    {
      gsl_bst_cmp_function *cmp = trav->rb_table->rb_compare;
      void *param = trav->rb_table->rb_param;
      rb_node *node = trav->rb_node;
      rb_node *i;

      trav->rb_height = 0;
      for (i = trav->rb_table->rb_root; i != node; )
        {
          if (trav->rb_height >= RB_MAX_HEIGHT)
            {
              GSL_ERROR_VOID ("traverser height exceeds maximum", GSL_ETABLE);
            }

          trav->rb_stack[trav->rb_height++] = i;
          i = i->rb_link[cmp (node->rb_data, i->rb_data, param) > 0];
        }
    }
}

static const gsl_bst_type rb_tree_type =
{
  "red-black",
  sizeof(rb_node),
  rb_init,
  rb_nodes,
  rb_insert,
  rb_find,
  rb_remove,
  rb_empty,

  rb_t_init,
  rb_t_first,
  rb_t_last,
  rb_t_find,
  rb_t_insert,
  rb_t_copy,
  rb_t_next,
  rb_t_prev,
  rb_t_cur,
  rb_t_replace
};

const gsl_bst_type * gsl_bst_rb = &rb_tree_type;
