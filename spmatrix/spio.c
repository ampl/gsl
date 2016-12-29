/* spio.c
 * 
 * Copyright (C) 2016 Patrick Alken
 * Copyright (C) 2016 Alexis Tantet
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

/*
gsl_spmatrix_fprintf()
  Print sparse matrix to file in MatrixMarket format:

M  N  NNZ
I1 J1 A(I1,J1)
...

Note that indices start at 1 and not 0
*/

int
gsl_spmatrix_fprintf(FILE *stream, const gsl_spmatrix *m,
                     const char *format)
{
  int status;

  /* print header */
  status = fprintf(stream, "%%%%MatrixMarket matrix coordinate real general\n");
  if (status < 0)
    {
      GSL_ERROR("fprintf failed for header", GSL_EFAILED);
    }

  /* print rows,columns,nnz */
  status = fprintf(stream, "%u\t%u\t%u\n",
                   (unsigned int) m->size1,
                   (unsigned int) m->size2,
                   (unsigned int) m->nz);
  if (status < 0)
    {
      GSL_ERROR("fprintf failed for dimension header", GSL_EFAILED);
    }

  if (GSL_SPMATRIX_ISTRIPLET(m))
    {
      size_t n;

      for (n = 0; n < m->nz; ++n)
        {
          status = fprintf(stream, "%u\t%u\t", (unsigned int) m->i[n] + 1, (unsigned int) m->p[n] + 1);
          if (status < 0)
            {
              GSL_ERROR("fprintf failed", GSL_EFAILED);
            }

          status = fprintf(stream, format, m->data[n]);
          if (status < 0)
            {
              GSL_ERROR("fprintf failed", GSL_EFAILED);
            }

          status = putc('\n', stream);
          if (status == EOF)
            {
              GSL_ERROR("putc failed", GSL_EFAILED);
            }
        }
    }
  else if (GSL_SPMATRIX_ISCCS(m))
    {
      size_t j, p;

      for (j = 0; j < m->size2; ++j)
        {
          for (p = m->p[j]; p < m->p[j + 1]; ++p)
            {
              status = fprintf(stream, "%u\t%u\t", (unsigned int) m->i[p] + 1, (unsigned int) j + 1);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

              status = fprintf(stream, format, m->data[p]);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

              status = putc('\n', stream);
              if (status == EOF)
                {
                  GSL_ERROR("putc failed", GSL_EFAILED);
                }
            }
        }
    }
  else if (GSL_SPMATRIX_ISCRS(m))
    {
      size_t i, p;

      for (i = 0; i < m->size1; ++i)
        {
          for (p = m->p[i]; p < m->p[i + 1]; ++p)
            {
              status = fprintf(stream, "%u\t%u\t", (unsigned int) i + 1, (unsigned int) m->i[p] + 1);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

              status = fprintf(stream, format, m->data[p]);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

              status = putc('\n', stream);
              if (status == EOF)
                {
                  GSL_ERROR("putc failed", GSL_EFAILED);
                }
            }
        }
    }
  else
    {
      GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
    }

  return GSL_SUCCESS;
}

gsl_spmatrix *
gsl_spmatrix_fscanf(FILE *stream)
{
  gsl_spmatrix *m;
  unsigned int size1, size2, nz;
  char buf[1024];
  int found_header = 0;

  /* read file until we find rows,cols,nz header */
  while (fgets(buf, 1024, stream) != NULL)
    {
      int c;

      /* skip comments */
      if (*buf == '%')
        continue;

      c = sscanf(buf, "%u %u %u",
                 &size1, &size2, &nz);
      if (c == 3)
        {
          found_header = 1;
          break;
        }
    }

  if (!found_header)
    {
      GSL_ERROR_NULL ("fscanf failed reading header", GSL_EFAILED);
    }

  m = gsl_spmatrix_alloc_nzmax((size_t) size1, (size_t) size2, (size_t) nz, GSL_SPMATRIX_TRIPLET);
  if (!m)
    {
      GSL_ERROR_NULL ("error allocating m", GSL_ENOMEM);
    }

  {
    unsigned int i, j;
    double val;

    while (fgets(buf, 1024, stream) != NULL)
      {
        int c = sscanf(buf, "%u %u %lg", &i, &j, &val);
        if (c < 3 || (i == 0) || (j == 0))
          {
            GSL_ERROR_NULL ("error in input file format", GSL_EFAILED);
          }
        else if ((i > size1) || (j > size2))
          {
            GSL_ERROR_NULL ("element exceeds matrix dimensions", GSL_EBADLEN);
          }
        else
          {
            /* subtract 1 from (i,j) since indexing starts at 1 */
            gsl_spmatrix_set(m, i - 1, j - 1, val);
          }
      }
  }

  return m;
}

int
gsl_spmatrix_fwrite(FILE *stream, const gsl_spmatrix *m)
{
  size_t items;

  /* write header: size1, size2, nz */

  items = fwrite(&(m->size1), sizeof(size_t), 1, stream);
  if (items != 1)
    {
      GSL_ERROR("fwrite failed on size1", GSL_EFAILED);
    }

  items = fwrite(&(m->size2), sizeof(size_t), 1, stream);
  if (items != 1)
    {
      GSL_ERROR("fwrite failed on size2", GSL_EFAILED);
    }

  items = fwrite(&(m->nz), sizeof(size_t), 1, stream);
  if (items != 1)
    {
      GSL_ERROR("fwrite failed on nz", GSL_EFAILED);
    }

  /* write m->i and m->data which are size nz in all storage formats */

  items = fwrite(m->i, sizeof(size_t), m->nz, stream);
  if (items != m->nz)
    {
      GSL_ERROR("fwrite failed on row indices", GSL_EFAILED);
    }

  items = fwrite(m->data, sizeof(double), m->nz, stream);
  if (items != m->nz)
    {
      GSL_ERROR("fwrite failed on data", GSL_EFAILED);
    }

  if (GSL_SPMATRIX_ISTRIPLET(m))
    {
      items = fwrite(m->p, sizeof(size_t), m->nz, stream);
      if (items != m->nz)
        {
          GSL_ERROR("fwrite failed on column indices", GSL_EFAILED);
        }
    }
  else if (GSL_SPMATRIX_ISCCS(m))
    {
      items = fwrite(m->p, sizeof(size_t), m->size2 + 1, stream);
      if (items != m->size2 + 1)
        {
          GSL_ERROR("fwrite failed on column indices", GSL_EFAILED);
        }
    }
  else if (GSL_SPMATRIX_ISCRS(m))
    {
      items = fwrite(m->p, sizeof(size_t), m->size1 + 1, stream);
      if (items != m->size1 + 1)
        {
          GSL_ERROR("fwrite failed on column indices", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

int
gsl_spmatrix_fread(FILE *stream, gsl_spmatrix *m)
{
  size_t size1, size2, nz;
  size_t items;

  /* read header: size1, size2, nz */

  items = fread(&size1, sizeof(size_t), 1, stream);
  if (items != 1)
    {
      GSL_ERROR("fread failed on size1", GSL_EFAILED);
    }

  items = fread(&size2, sizeof(size_t), 1, stream);
  if (items != 1)
    {
      GSL_ERROR("fread failed on size2", GSL_EFAILED);
    }

  items = fread(&nz, sizeof(size_t), 1, stream);
  if (items != 1)
    {
      GSL_ERROR("fread failed on nz", GSL_EFAILED);
    }

  if (m->size1 != size1)
    {
      GSL_ERROR("matrix has wrong size1", GSL_EBADLEN);
    }
  else if (m->size2 != size2)
    {
      GSL_ERROR("matrix has wrong size2", GSL_EBADLEN);
    }
  else if (nz > m->nzmax)
    {
      GSL_ERROR("matrix nzmax is too small", GSL_EBADLEN);
    }
  else
    {
      /* read m->i and m->data arrays, which are size nz for all formats */

      items = fread(m->i, sizeof(size_t), nz, stream);
      if (items != nz)
        {
          GSL_ERROR("fread failed on row indices", GSL_EFAILED);
        }

      items = fread(m->data, sizeof(double), nz, stream);
      if (items != nz)
        {
          GSL_ERROR("fread failed on data", GSL_EFAILED);
        }

      m->nz = nz;

      if (GSL_SPMATRIX_ISTRIPLET(m))
        {
          items = fread(m->p, sizeof(size_t), nz, stream);
          if (items != nz)
            {
              GSL_ERROR("fread failed on column indices", GSL_EFAILED);
            }

          /* build binary search tree for m */
          gsl_spmatrix_tree_rebuild(m);
        }
      else if (GSL_SPMATRIX_ISCCS(m))
        {
          items = fread(m->p, sizeof(size_t), size2 + 1, stream);
          if (items != size2 + 1)
            {
              GSL_ERROR("fread failed on row pointers", GSL_EFAILED);
            }
        }
      else if (GSL_SPMATRIX_ISCRS(m))
        {
          items = fread(m->p, sizeof(size_t), size1 + 1, stream);
          if (items != size1 + 1)
            {
              GSL_ERROR("fread failed on column pointers", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
