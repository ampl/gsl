/* test.c
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
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bst.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_test.h>

enum array_order
{
  ORD_RANDOM = 0,               /* random order */
  ORD_ASCENDING,                /* ascending order */
  ORD_DESCENDING,               /* descending order */
  ORD_BALANCED,                 /* balanced tree order */
  ORD_ZIGZAG,                   /* zig-zag order */
  ORD_ASCENDING_SHIFTED,        /* ascending from middle, then beginning */
  ORD_END_NODUP,                /* end of no-duplicate ordering */

  ORD_RANDOM_DUP                /* random order with duplicates */
};

/* fill array[] with random integers in [lower,upper] with duplicates allowed */
static void
random_integers(const size_t n, const int lower, const int upper,
                int array[], gsl_rng * r)
{
  size_t i;

  for (i = 0; i < n; ++i)
    array[i] = (int) ((upper - lower) * gsl_rng_uniform(r) + lower);
}

/* fills array[] with a random permutation of the integers between 0 and n - 1 */
static void
random_permuted_integers (const size_t n, int array[], gsl_rng * r)
{
  size_t i;

  for (i = 0; i < n; i++)
    array[i] = i;

  for (i = 0; i < n; i++)
    {
      size_t j = i + (unsigned) (gsl_rng_uniform(r) * (n - i));
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
}

static int
compare_ints(const void *pa, const void *pb, void *params)
{
  const int *a = pa;
  const int *b = pb;

  (void) params;

  return (*a < *b) ? -1 : (*a > *b);
}

/* Generates a list of integers that produce a balanced tree when
   inserted in order into a binary tree in the usual way.
   |min| and |max| inclusively bound the values to be inserted.
   Output is deposited starting at |*array|. */
static void
gen_balanced_tree (const int min, const int max, int **array)
{
  int i;

  if (min > max)
    return;

  i = (min + max + 1) / 2;
  *(*array)++ = i;
  gen_balanced_tree (min, i - 1, array);
  gen_balanced_tree (i + 1, max, array);
}

/* generates a permutation of the integers |0| to |n - 1| */
static void
gen_int_array (const size_t n, const enum array_order order, int array[], gsl_rng * r)
{
  size_t i;

  switch (order)
    {
      case ORD_RANDOM:
        random_permuted_integers (n, array, r);
        break;

      case ORD_ASCENDING:
        for (i = 0; i < n; i++)
          array[i] = i;
        break;

      case ORD_DESCENDING:
        for (i = 0; i < n; i++)
          array[i] = n - i - 1;
        break;

      case ORD_BALANCED:
        gen_balanced_tree (0, n - 1, &array);
        break;

      case ORD_ZIGZAG:
        for (i = 0; i < n; i++)
          {
            if (i % 2 == 0)
              array[i] = i / 2;
            else
              array[i] = n - i / 2 - 1;
          }
        break;

      case ORD_ASCENDING_SHIFTED:
        for (i = 0; i < n; i++)
          {
            array[i] = i + n / 2;
            if ((size_t) array[i] >= n)
              array[i] -= n;
          }
        break;

      case ORD_RANDOM_DUP:
        random_integers(n, -10, 10, array, r);
        break;

      default:
        assert (0);
    }
}

static void
check_traverser(const size_t n, const enum array_order order, gsl_bst_trav * trav, int data,
                const char *desc, const gsl_bst_workspace * w)
{
  int *prev, *cur, *next;

  prev = gsl_bst_trav_prev(trav);
  if (prev != NULL)
    {
      gsl_test(*prev > data, "bst %s[n=%zu,order=%d] %s traverser ahead of %d, but should be ahead of %d",
               gsl_bst_name(w), n, order, desc, *prev, data);
    }
  gsl_bst_trav_next(trav);

  cur = gsl_bst_trav_cur(trav);
  gsl_test(*cur != data, "bst %s[n=%zu,order=%d] %s traverser at %d, but should be at %d",
           gsl_bst_name(w), n, order, desc, *cur, data);

  next = gsl_bst_trav_next(trav);
  if (next != NULL)
    {
      gsl_test(*next < data, "bst %s[n=%zu,order=%d] %s traverser behind %d, but should be behind %d",
               gsl_bst_name(w), n, order, desc, *next, data);
    }
  gsl_bst_trav_prev(trav);
}

static void
test_bst_int(const size_t n, const gsl_bst_type * T, const enum array_order order, gsl_rng * r)
{
  int *data = malloc(n * sizeof(int));
  int *data_delete = malloc(n * sizeof(int));
  int *sorted_data = malloc(n * sizeof(int));
  gsl_bst_workspace * w = gsl_bst_alloc(T, NULL, compare_ints, NULL);
  gsl_bst_trav trav;
  int *p;
  int i;
  size_t nodes;

  /* generate data to be inserted in tree */
  gen_int_array(n, order, data, r);

  for (i = 0; i < (int) n; ++i)
    sorted_data[i] = data[i];

  gsl_sort_int(sorted_data, 1, n);

  if (order != ORD_RANDOM_DUP)
    {
      /* generate random order to delete data from tree */
      gen_int_array(n, ORD_RANDOM, data_delete, r);
    }
  else
    {
      for (i = 0; i < (int) n; ++i)
        data_delete[i] = sorted_data[i];
    }

  /* insert data */
  for (i = 0; i < (int) n; ++i)
    {
      p = gsl_bst_insert(&data[i], w);
      gsl_test(p != NULL, "bst_int %s[n=%zu,order=%d] insert i=%d", gsl_bst_name(w), n, order, i);
    }

  if (order != ORD_RANDOM_DUP)
    {
      nodes = gsl_bst_nodes(w);
      gsl_test(nodes != n, "bst_int %s[n=%zu,order=%d] after insertion count = %zu/%zu",
               gsl_bst_name(w), n, order, nodes, n);
    }

  /* test data was inserted and can be found */
  for (i = 0; i < (int) n; ++i)
    {
      p = gsl_bst_find(&data[i], w);
      gsl_test(*p != data[i], "bst_int %s[n=%zu,order=%d] find [%d,%d]",
               gsl_bst_name(w), n, order, *p, data[i]);

      p = gsl_bst_trav_find(&data[i], &trav, w);
      gsl_test(p == NULL, "bst_int %s[n=%zu,order=%d] trav_find unable to find item %d",
               gsl_bst_name(w), n, order, data[i]);

      check_traverser(n, order, &trav, data[i], "post-insertion", w);
    }

  /* traverse tree in-order */
  p = gsl_bst_trav_first(&trav, w);
  i = 0;
  while (p != NULL)
    {
      int *q = gsl_bst_trav_cur(&trav);

      gsl_test(*p != sorted_data[i], "bst_int %s[n=%zu,order=%d] traverse i=%d [%d,%d]",
               gsl_bst_name(w), n, order, i, *p, sorted_data[i]);

      gsl_test(*p != *q, "bst_int %s[n=%zu,order=%d] traverse cur i=%d [%d,%d]",
               gsl_bst_name(w), n, order, i, *p, *q);

      p = gsl_bst_trav_next(&trav);
      ++i;
    }

  gsl_test(i != (int) n, "bst_int %s[n=%zu,order=%d] traverse number=%d",
           gsl_bst_name(w), n, order, i);

  /* traverse tree in reverse order */
  p = gsl_bst_trav_last(&trav, w);
  i = n - 1;
  while (p != NULL)
    {
      int *q = gsl_bst_trav_cur(&trav);

      gsl_test(*p != sorted_data[i], "bst_int %s[n=%zu,order=%d] traverse reverse i=%d [%d,%d]",
               gsl_bst_name(w), n, order, i, *p, sorted_data[i]);

      gsl_test(*p != *q, "bst_int %s[n=%zu,order=%d] traverse reverse cur i=%d [%d,%d]",
               gsl_bst_name(w), n, order, i, *p, *q);

      p = gsl_bst_trav_prev(&trav);
      --i;
    }

  gsl_test(i != -1, "bst_int %s[n=%zu,order=%d] traverse reverse number=%d",
           gsl_bst_name(w), n, order, i);

  /* test traversal during tree modifications */
  for (i = 0; i < (int) n; ++i)
    {
      gsl_bst_trav x, y, z;

      gsl_bst_trav_find(&data[i], &x, w);
      check_traverser(n, order, &x, data[i], "pre-deletion", w);

      if (data[i] == data_delete[i])
        continue;

      p = gsl_bst_remove(&data_delete[i], w);
      gsl_test(*p != data_delete[i], "bst_int %s[n=%zu,order=%d] remove i=%d [%d,%d]",
               gsl_bst_name(w), n, order, i, *p, data_delete[i]);

      p = gsl_bst_trav_copy(&y, &x);
      gsl_test(*p != data[i], "bst_int %s[n=%zu,order=%d] copy i=%d [%d,%d]",
               gsl_bst_name(w), n, order, i, *p, data[i]);

      /* re-insert item */
      p = gsl_bst_trav_insert(&data_delete[i], &z, w);

      check_traverser(n, order, &x, data[i], "post-deletion", w);
      check_traverser(n, order, &y, data[i], "copied", w);
      check_traverser(n, order, &z, data_delete[i], "insertion", w);

#if 0
      /* delete again */
      gsl_bst_remove(&data[i], w);
#endif
    }

  /* emmpty tree */
  gsl_bst_empty(w);

  nodes = gsl_bst_nodes(w);
  gsl_test(nodes != 0, "bst_int %s[n=%zu,order=%d] empty count = %zu",
           gsl_bst_name(w), n, order, nodes);

  gsl_bst_free(w);
  free(data);
  free(data_delete);
  free(sorted_data);
}

static void
test_bst(const gsl_bst_type * T, gsl_rng * r)
{
  enum array_order order;

  for (order = 0; order < ORD_END_NODUP; ++order)
    {
      test_bst_int(50, T, order, r);
      test_bst_int(100, T, order, r);
      test_bst_int(500, T, order, r);
    }
}

int
main(void)
{
  gsl_rng * r = gsl_rng_alloc(gsl_rng_default);

  test_bst(gsl_bst_avl, r);
  test_bst(gsl_bst_rb, r);

  gsl_rng_free(r);

  exit (gsl_test_summary());
}
