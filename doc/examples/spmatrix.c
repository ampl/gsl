#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_spmatrix.h>

int
main()
{
  gsl_spmatrix *A = gsl_spmatrix_alloc(5, 4); /* triplet format */
  gsl_spmatrix *B, *C;
  size_t i, j;

  /* build the sparse matrix */
  gsl_spmatrix_set(A, 0, 2, 3.1);
  gsl_spmatrix_set(A, 0, 3, 4.6);
  gsl_spmatrix_set(A, 1, 0, 1.0);
  gsl_spmatrix_set(A, 1, 2, 7.2);
  gsl_spmatrix_set(A, 3, 0, 2.1);
  gsl_spmatrix_set(A, 3, 1, 2.9);
  gsl_spmatrix_set(A, 3, 3, 8.5);
  gsl_spmatrix_set(A, 4, 0, 4.1);

  printf("printing all matrix elements:\n");
  for (i = 0; i < 5; ++i)
    for (j = 0; j < 4; ++j)
      printf("A(%zu,%zu) = %g\n", i, j,
             gsl_spmatrix_get(A, i, j));

  /* print out elements in triplet format */
  printf("matrix in triplet format (i,j,Aij):\n");
  gsl_spmatrix_fprintf(stdout, A, "%.1f");

  /* convert to compressed column format */
  B = gsl_spmatrix_ccs(A);

  printf("matrix in compressed column format:\n");
  printf("i = [ ");
  for (i = 0; i < B->nz; ++i)
    printf("%d, ", B->i[i]);
  printf("]\n");

  printf("p = [ ");
  for (i = 0; i < B->size2 + 1; ++i)
    printf("%d, ", B->p[i]);
  printf("]\n");

  printf("d = [ ");
  for (i = 0; i < B->nz; ++i)
    printf("%g, ", B->data[i]);
  printf("]\n");

  /* convert to compressed row format */
  C = gsl_spmatrix_crs(A);

  printf("matrix in compressed row format:\n");
  printf("i = [ ");
  for (i = 0; i < C->nz; ++i)
    printf("%d, ", C->i[i]);
  printf("]\n");

  printf("p = [ ");
  for (i = 0; i < C->size1 + 1; ++i)
    printf("%d, ", C->p[i]);
  printf("]\n");

  printf("d = [ ");
  for (i = 0; i < C->nz; ++i)
    printf("%g, ", C->data[i]);
  printf("]\n");

  gsl_spmatrix_free(A);
  gsl_spmatrix_free(B);
  gsl_spmatrix_free(C);

  return 0;
}
