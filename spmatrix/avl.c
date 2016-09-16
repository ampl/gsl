/* avl.c
 * 
 * Copyright (C) 1998, 1999, 2000, 2001, 2002, 2004 Free Software
 * Foundation, Inc.
 * Copyright (C) 2014 Patrick Alken
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

/*
 * This code is originally from GNU libavl. The memory management
 * was slightly modified for use with preallocating GSL sparse matrices
 *
 * The allocator->libavl_malloc function is called only for creating
 * a new avl_node (tree node). This allows GSL to preallocate some number
 * of avl_node structs and then return pointers to them while the tree
 * is being assembled, avoiding multiple malloc calls
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

/* Function types. */
typedef int avl_comparison_func (const void *avl_a, const void *avl_b,
                                 void *avl_param);
typedef void avl_item_func (void *avl_item, void *avl_param);
typedef void *avl_copy_func (void *avl_item, void *avl_param);

/* Memory allocator. */
struct libavl_allocator
  {
    void *(*libavl_malloc) (size_t libavl_size, void *param);
    void (*libavl_free) (void *libavl_block, void *param);
  };

/* Default memory allocator. */
static struct libavl_allocator avl_allocator_default;
static void *avl_malloc (size_t, void *param);
static void avl_free (void *, void *param);

/* Maximum AVL tree height. */
#ifndef AVL_MAX_HEIGHT
#define AVL_MAX_HEIGHT 92
#endif

/* An AVL tree node. */
struct avl_node
  {
    struct avl_node *avl_link[2];  /* Subtrees. */
    void *avl_data;                /* Pointer to data. */
    signed char avl_balance;       /* Balance factor. */
  };

/* Tree data structure. */
struct avl_table
  {
    struct avl_node *avl_root;          /* Tree's root. */
    avl_comparison_func *avl_compare;   /* Comparison function. */
    void *avl_param;                    /* Extra argument to |avl_compare|. */
    struct libavl_allocator *avl_alloc; /* Memory allocator. */
    size_t avl_count;                   /* Number of items in tree. */
    unsigned long avl_generation;       /* Generation number. */
  };

/* Table functions. */
static struct avl_table *avl_create (avl_comparison_func *, void *,
                                     struct libavl_allocator *);
static struct avl_table *avl_copy (const struct avl_table *, avl_copy_func *,
                                   avl_item_func *, struct libavl_allocator *);
static void avl_empty (struct avl_table *, avl_item_func *);
static void avl_destroy (struct avl_table *, avl_item_func *);
static void **avl_probe (struct avl_table *, void *);
static void *avl_insert (struct avl_table *, void *);
static void *avl_replace (struct avl_table *, void *);
static void *avl_delete (struct avl_table *, const void *);
static void *avl_find (const struct avl_table *, const void *);

/* Creates and returns a new table
   with comparison function |compare| using parameter |param|
   and memory allocator |allocator|.
   Returns |NULL| if memory allocation failed. */
static struct avl_table *
avl_create (avl_comparison_func *compare, void *param,
            struct libavl_allocator *allocator)
{
  struct avl_table *tree;

  if (allocator == NULL)
    allocator = &avl_allocator_default;

  /*tree = allocator->libavl_malloc (allocator, sizeof *tree);*/
  tree = malloc(sizeof *tree);
  if (tree == NULL)
    return NULL;

  tree->avl_root = NULL;
  tree->avl_compare = compare;
  tree->avl_param = param;
  tree->avl_alloc = allocator;
  tree->avl_count = 0;
  tree->avl_generation = 0;

  return tree;
}

/* Search |tree| for an item matching |item|, and return it if found.
   Otherwise return |NULL|. */
static void *
avl_find (const struct avl_table *tree, const void *item)
{
  const struct avl_node *p;

  for (p = tree->avl_root; p != NULL; )
    {
      int cmp = tree->avl_compare (item, p->avl_data, tree->avl_param);

      if (cmp < 0)
        p = p->avl_link[0];
      else if (cmp > 0)
        p = p->avl_link[1];
      else /* |cmp == 0| */
        return p->avl_data;
    }

  return NULL;
}

/* Inserts |item| into |tree| and returns a pointer to |item|'s address.
   If a duplicate item is found in the tree,
   returns a pointer to the duplicate without inserting |item|.
   Returns |NULL| in case of memory allocation failure. */
static void **
avl_probe (struct avl_table *tree, void *item)
{
  struct avl_node *y, *z; /* Top node to update balance factor, and parent. */
  struct avl_node *p, *q; /* Iterator, and parent. */
  struct avl_node *n;     /* Newly inserted node. */
  struct avl_node *w;     /* New root of rebalanced subtree. */
  int dir;                /* Direction to descend. */

  unsigned char da[AVL_MAX_HEIGHT]; /* Cached comparison results. */
  int k = 0;              /* Number of cached results. */

  z = (struct avl_node *) &tree->avl_root;
  y = tree->avl_root;
  dir = 0;
  for (q = z, p = y; p != NULL; q = p, p = p->avl_link[dir])
    {
      int cmp = tree->avl_compare (item, p->avl_data, tree->avl_param);
      if (cmp == 0)
        return &p->avl_data;

      if (p->avl_balance != 0)
        z = q, y = p, k = 0;
      da[k++] = dir = cmp > 0;
    }

  n = q->avl_link[dir] =
    tree->avl_alloc->libavl_malloc (sizeof *n, tree->avl_param);
  if (n == NULL)
    return NULL;

  tree->avl_count++;
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
      struct avl_node *x = y->avl_link[0];
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
      struct avl_node *x = y->avl_link[1];
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

  tree->avl_generation++;
  return &n->avl_data;
}

/* Inserts |item| into |table|.
   Returns |NULL| if |item| was successfully inserted
   or if a memory allocation error occurred.
   Otherwise, returns the duplicate item. */
static void *
avl_insert (struct avl_table *table, void *item)
{
  void **p = avl_probe (table, item);
  return p == NULL || *p == item ? NULL : *p;
}

/* Inserts |item| into |table|, replacing any duplicate item.
   Returns |NULL| if |item| was inserted without replacing a duplicate,
   or if a memory allocation error occurred.
   Otherwise, returns the item that was replaced. */
static void *
avl_replace (struct avl_table *table, void *item)
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

/* Deletes from |tree| and returns an item matching |item|.
   Returns a null pointer if no matching item found. */
static void *
avl_delete (struct avl_table *tree, const void *item)
{
  /* Stack of nodes. */
  struct avl_node *pa[AVL_MAX_HEIGHT]; /* Nodes. */
  unsigned char da[AVL_MAX_HEIGHT];    /* |avl_link[]| indexes. */
  int k;                               /* Stack pointer. */

  struct avl_node *p;   /* Traverses tree to find node to delete. */
  int cmp;              /* Result of comparison between |item| and |p|. */

  k = 0;
  p = (struct avl_node *) &tree->avl_root;
  for (cmp = -1; cmp != 0;
       cmp = tree->avl_compare (item, p->avl_data, tree->avl_param))
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
      struct avl_node *r = p->avl_link[1];
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
          struct avl_node *s;
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

  tree->avl_alloc->libavl_free (p, tree->avl_param);

  while (--k > 0)
    {
      struct avl_node *y = pa[k];

      if (da[k] == 0)
        {
          y->avl_balance++;
          if (y->avl_balance == +1)
            break;
          else if (y->avl_balance == +2)
            {
              struct avl_node *x = y->avl_link[1];
              if (x->avl_balance == -1)
                {
                  struct avl_node *w;
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
              struct avl_node *x = y->avl_link[0];
              if (x->avl_balance == +1)
                {
                  struct avl_node *w;
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

  tree->avl_count--;
  tree->avl_generation++;
  return (void *) item;
}

/* Destroys |new| with |avl_destroy (new, destroy)|,
   first setting right links of nodes in |stack| within |new|
   to null pointers to avoid touching uninitialized data. */
static void
copy_error_recovery (struct avl_node **stack, int height,
                     struct avl_table *new, avl_item_func *destroy)
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
static struct avl_table *
avl_copy (const struct avl_table *org, avl_copy_func *copy,
          avl_item_func *destroy, struct libavl_allocator *allocator)
{
  struct avl_node *stack[2 * (AVL_MAX_HEIGHT + 1)];
  int height = 0;

  struct avl_table *new;
  const struct avl_node *x;
  struct avl_node *y;

  new = avl_create (org->avl_compare, org->avl_param,
                    allocator != NULL ? allocator : org->avl_alloc);
  if (new == NULL)
    return NULL;
  new->avl_count = org->avl_count;
  if (new->avl_count == 0)
    return new;

  x = (const struct avl_node *) &org->avl_root;
  y = (struct avl_node *) &new->avl_root;
  for (;;)
    {
      while (x->avl_link[0] != NULL)
        {
          y->avl_link[0] =
            new->avl_alloc->libavl_malloc (sizeof *y->avl_link[0],
                                           new->avl_param);
          if (y->avl_link[0] == NULL)
            {
              if (y != (struct avl_node *) &new->avl_root)
                {
                  y->avl_data = NULL;
                  y->avl_link[1] = NULL;
                }

              copy_error_recovery (stack, height, new, destroy);
              return NULL;
            }

          stack[height++] = (struct avl_node *) x;
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
                new->avl_alloc->libavl_malloc (sizeof *y->avl_link[1],
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

/* empty tree (delete all nodes) but do not free the tree itself */
static void
avl_empty (struct avl_table *tree, avl_item_func *destroy)
{
  struct avl_node *p, *q;

  for (p = tree->avl_root; p != NULL; p = q)
    if (p->avl_link[0] == NULL)
      {
        q = p->avl_link[1];
        if (destroy != NULL && p->avl_data != NULL)
          destroy (p->avl_data, tree->avl_param);
        tree->avl_alloc->libavl_free (p, tree->avl_param);
      }
    else
      {
        q = p->avl_link[0];
        p->avl_link[0] = q->avl_link[1];
        q->avl_link[1] = p;
      }

  tree->avl_root = NULL;
  tree->avl_count = 0;
  tree->avl_generation = 0;
}

/* Frees storage allocated for |tree|.
   If |destroy != NULL|, applies it to each data item in inorder. */
static void
avl_destroy (struct avl_table *tree, avl_item_func *destroy)
{
  avl_empty(tree, destroy);
  free(tree);
}

/* Allocates |size| bytes of space using |malloc()|.
   Returns a null pointer if allocation fails. */
static void *
avl_malloc (size_t size, void *param)
{
  (void)param; /* avoid unused parameter warning */
  return malloc (size);
}

/* Frees |block|. */
static void
avl_free (void *block, void *param)
{
  (void)param; /* avoid unused parameter warning */
  free (block);
}

/* Default memory allocator that uses |malloc()| and |free()|. */
static struct libavl_allocator avl_allocator_default =
  {
    avl_malloc,
    avl_free
  };
