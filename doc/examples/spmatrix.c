#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_spmatrix.h>

int
main()
{
  gsl_spmatrix *A = gsl_spmatrix_alloc(5, 4); /* triplet format */
  gsl_spmatrix *C;
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
  for (i = 0; i < A->nz; ++i)
    printf("(%zu, %zu, %.1f)\n", A->i[i], A->p[i], A->data[i]);

  /* convert to compressed column format */
  C = gsl_spmatrix_compcol(A);

  printf("matrix in compressed column format:\n");
  printf("i = [ ");
  for (i = 0; i < C->nz; ++i)
    printf("%zu, ", C->i[i]);
  printf("]\n");

  printf("p = [ ");
  for (i = 0; i < C->size2 + 1; ++i)
    printf("%zu, ", C->p[i]);
  printf("]\n");

  printf("d = [ ");
  for (i = 0; i < C->nz; ++i)
    printf("%g, ", C->data[i]);
  printf("]\n");

  gsl_spmatrix_free(A);
  gsl_spmatrix_free(C);

  return 0;
}
