#include <stdio.h>
#include <gsl/gsl_block.h>

int
main (void)
{
  gsl_block * b = gsl_block_alloc (100);
  
  printf ("length of block = "F_ZU"\n", b->size);
  printf ("block data address = %p\n", b->data);

  gsl_block_free (b);
  return 0;
}
