/* spmatrix/file_source.c
 * 
 * Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Patrick Alken
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

/*
gsl_spmatrix_fprintf()
  Print sparse matrix to file in MatrixMarket format:

M  N  NNZ
I1 J1 A(I1,J1)
...

Note that indices start at 1 and not 0
*/

int
FUNCTION (gsl_spmatrix, fprintf) (FILE * stream, const TYPE (gsl_spmatrix) * m,
                                   const char * format)
{
  int status;

  /* print header */

#if defined(BASE_GSL_COMPLEX_LONG) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)
  status = fprintf(stream, "%%%%MatrixMarket matrix coordinate complex general\n");
#else
  status = fprintf(stream, "%%%%MatrixMarket matrix coordinate real general\n");
#endif

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

  if (GSL_SPMATRIX_ISCOO(m))
    {
      size_t n;

      for (n = 0; n < m->nz; ++n)
        {
          status = fprintf(stream, "%d\t%d\t", m->i[n] + 1, m->p[n] + 1);
          if (status < 0)
            {
              GSL_ERROR("fprintf failed", GSL_EFAILED);
            }

          status = fprintf(stream, format, m->data[MULTIPLICITY * n]);
          if (status < 0)
            {
              GSL_ERROR("fprintf failed", GSL_EFAILED);
            }

#if defined(BASE_GSL_COMPLEX_LONG) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

          status = putc('\t', stream);
          if (status == EOF)
            {
              GSL_ERROR("putc failed", GSL_EFAILED);
            }

          status = fprintf(stream, format, m->data[MULTIPLICITY * n + 1]);
          if (status < 0)
            {
              GSL_ERROR("fprintf failed", GSL_EFAILED);
            }

#endif

          status = putc('\n', stream);
          if (status == EOF)
            {
              GSL_ERROR("putc failed", GSL_EFAILED);
            }
        }
    }
  else if (GSL_SPMATRIX_ISCSC(m))
    {
      size_t j;
      int p;

      for (j = 0; j < m->size2; ++j)
        {
          for (p = m->p[j]; p < m->p[j + 1]; ++p)
            {
              status = fprintf(stream, "%d\t%u\t", m->i[p] + 1, (unsigned int) (j + 1));
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

              status = fprintf(stream, format, m->data[MULTIPLICITY * p]);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

#if defined(BASE_GSL_COMPLEX_LONG) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

              status = putc('\t', stream);
              if (status == EOF)
                {
                  GSL_ERROR("putc failed", GSL_EFAILED);
                }

              status = fprintf(stream, format, m->data[MULTIPLICITY * p + 1]);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

#endif

              status = putc('\n', stream);
              if (status == EOF)
                {
                  GSL_ERROR("putc failed", GSL_EFAILED);
                }
            }
        }
    }
  else if (GSL_SPMATRIX_ISCSR(m))
    {
      size_t i;
      int p;

      for (i = 0; i < m->size1; ++i)
        {
          for (p = m->p[i]; p < m->p[i + 1]; ++p)
            {
              status = fprintf(stream, "%u\t%d\t", (unsigned int) (i + 1), m->i[p] + 1);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

              status = fprintf(stream, format, m->data[MULTIPLICITY * p]);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

#if defined(BASE_GSL_COMPLEX_LONG) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

              status = putc('\t', stream);
              if (status == EOF)
                {
                  GSL_ERROR("putc failed", GSL_EFAILED);
                }

              status = fprintf(stream, format, m->data[MULTIPLICITY * p + 1]);
              if (status < 0)
                {
                  GSL_ERROR("fprintf failed", GSL_EFAILED);
                }

#endif

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

TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, fscanf) (FILE * stream)
{
  TYPE (gsl_spmatrix) * m;
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

      c = sscanf(buf, "%u %u %u", &size1, &size2, &nz);
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

  m = FUNCTION (gsl_spmatrix, alloc_nzmax) ((size_t) size1, (size_t) size2, (size_t) nz, GSL_SPMATRIX_COO);
  if (!m)
    {
      GSL_ERROR_NULL ("error allocating m", GSL_ENOMEM);
    }

#if defined(BASE_GSL_COMPLEX_LONG) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

  {
    unsigned int i, j;
    ATOMIC_IO xr, xi;
    BASE x;

    while (fgets(buf, 1024, stream) != NULL)
      {
        int c = sscanf(buf, "%u %u " IN_FORMAT " " IN_FORMAT, &i, &j, &xr, &xi);
        if (c < 4 || i == 0 || j == 0)
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
            GSL_REAL(x) = xr;
            GSL_IMAG(x) = xi;
            FUNCTION (gsl_spmatrix, set) (m, i - 1, j - 1, x);
          }
      }
  }

#else

  {
    unsigned int i, j;
    ATOMIC_IO tmp;

    while (fgets(buf, 1024, stream) != NULL)
      {
        int c = sscanf(buf, "%u %u " IN_FORMAT, &i, &j, &tmp);
        if (c < 3 || i == 0 || j == 0)
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
            FUNCTION (gsl_spmatrix, set) (m, i - 1, j - 1, tmp);
          }
      }
  }

#endif

  return m;
}

int
FUNCTION (gsl_spmatrix, fwrite) (FILE * stream, const TYPE (gsl_spmatrix) * m)
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

  items = fwrite(m->i, sizeof(int), m->nz, stream);
  if (items != m->nz)
    {
      GSL_ERROR("fwrite failed on row indices", GSL_EFAILED);
    }

  items = fwrite(m->data, MULTIPLICITY * sizeof(ATOMIC), m->nz, stream);
  if (items != m->nz)
    {
      GSL_ERROR("fwrite failed on data", GSL_EFAILED);
    }

  if (GSL_SPMATRIX_ISCOO(m))
    {
      items = fwrite(m->p, sizeof(int), m->nz, stream);
      if (items != m->nz)
        {
          GSL_ERROR("fwrite failed on column indices", GSL_EFAILED);
        }
    }
  else if (GSL_SPMATRIX_ISCSC(m))
    {
      items = fwrite(m->p, sizeof(int), m->size2 + 1, stream);
      if (items != m->size2 + 1)
        {
          GSL_ERROR("fwrite failed on column indices", GSL_EFAILED);
        }
    }
  else if (GSL_SPMATRIX_ISCSR(m))
    {
      items = fwrite(m->p, sizeof(int), m->size1 + 1, stream);
      if (items != m->size1 + 1)
        {
          GSL_ERROR("fwrite failed on column indices", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_spmatrix, fread) (FILE * stream, TYPE (gsl_spmatrix) * m)
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

      items = fread(m->i, sizeof(int), nz, stream);
      if (items != nz)
        {
          GSL_ERROR("fread failed on row indices", GSL_EFAILED);
        }

      items = fread(m->data, MULTIPLICITY * sizeof(ATOMIC), nz, stream);
      if (items != nz)
        {
          GSL_ERROR("fread failed on data", GSL_EFAILED);
        }

      m->nz = nz;

      if (GSL_SPMATRIX_ISCOO(m))
        {
          items = fread(m->p, sizeof(int), nz, stream);
          if (items != nz)
            {
              GSL_ERROR("fread failed on column indices", GSL_EFAILED);
            }

          /* build binary search tree for m */
          FUNCTION (gsl_spmatrix, tree_rebuild) (m);
        }
      else if (GSL_SPMATRIX_ISCSC(m))
        {
          items = fread(m->p, sizeof(int), size2 + 1, stream);
          if (items != size2 + 1)
            {
              GSL_ERROR("fread failed on row pointers", GSL_EFAILED);
            }
        }
      else if (GSL_SPMATRIX_ISCSR(m))
        {
          items = fread(m->p, sizeof(int), size1 + 1, stream);
          if (items != size1 + 1)
            {
              GSL_ERROR("fread failed on column pointers", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
