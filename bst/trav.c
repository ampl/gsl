/* trav.c
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

int
gsl_bst_trav_init(gsl_bst_trav * trav, const gsl_bst_workspace * w)
{
  int status = (w->type->trav_init)((void *) &trav->trav_data, (const void *) &w->table);
  trav->type = w->type;
  return status;
}

/*
gsl_bst_trav_first()
  Initialize traverser to least-valued item in tree and return
a pointer to it. Return NULL if tree has no nodes.
*/

void *
gsl_bst_trav_first(gsl_bst_trav * trav, const gsl_bst_workspace * w)
{
  trav->type = w->type;
  return (w->type->trav_first)((void *) &trav->trav_data, (const void *) &w->table);
}

/*
gsl_bst_trav_last()
  Initializes |trav| for |tree| and selects and returns a pointer
to its greatest-valued item. Returns NULL if |tree| contains no nodes.
*/

void *
gsl_bst_trav_last (gsl_bst_trav * trav, const gsl_bst_workspace * w)
{
  trav->type = w->type;
  return (w->type->trav_last)((void * ) &trav->trav_data, (const void *) &w->table);
}

/*
gsl_bst_trav_find()
  Searches for |item| in tree. If found, initializes |trav| to
the item found and returns the item as well.  If there is no matching
item, initializes |trav| to the null item and returns |NULL|.
*/

void *
gsl_bst_trav_find (const void * item, gsl_bst_trav * trav, const gsl_bst_workspace * w)
{
  trav->type = w->type;
  return (w->type->trav_find)(item, (void * ) &trav->trav_data, (const void *) &w->table);
}

/*
gsl_bst_trav_insert()
  Attempts to insert |item| into tree.  If |item| is inserted
successfully, it is returned and |trav| is initialized to its location.
If a duplicate is found, it is returned and |trav| is initialized to
its location.  No replacement of the item occurs.
If a memory allocation failure occurs, |NULL| is returned and |trav|
is initialized to the null item.
*/

void *
gsl_bst_trav_insert (void * item, gsl_bst_trav * trav, gsl_bst_workspace * w)
{
  trav->type = w->type;
  return (w->type->trav_insert)(item, (void * ) &trav->trav_data, (void *) &w->table);
}

/*
gsl_bst_trav_copy()
  Copy traverser 'src' into 'dest'
*/

void *
gsl_bst_trav_copy(gsl_bst_trav * dest, const gsl_bst_trav * src)
{
  dest->type = src->type;
  return (src->type->trav_copy)((void * ) &dest->trav_data, (const void *) &src->trav_data);
}

/*
gsl_bst_trav_next()
  Update traverser to point to next sequential node, and
return a pointer to its data
*/

void *
gsl_bst_trav_next(gsl_bst_trav * trav)
{
  return (trav->type->trav_next)((void *) &trav->trav_data);
}

/*
gsl_bst_trav_prev()
  Update traverser to point to previous sequential node, and
return a pointer to its data
*/

void *
gsl_bst_trav_prev(gsl_bst_trav * trav)
{
  return (trav->type->trav_prev)((void *) &trav->trav_data);
}

/*
gsl_bst_trav_cur()
  Return a pointer to data of current traverser node
*/

void *
gsl_bst_trav_cur(const gsl_bst_trav * trav)
{
  return (trav->type->trav_cur)((const void *) &trav->trav_data);
}

/*
gsl_bst_trav_replace()
  Replace current item in trav with new_item and returns the item
replaced. The new item must not change the ordering of the tree.
*/

void *
gsl_bst_trav_replace (gsl_bst_trav * trav, void * new_item)
{
  return (trav->type->trav_replace)((void * ) &trav->trav_data, new_item);
}
