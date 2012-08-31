/* memmove.c -- copy memory.
   Copy LENGTH bytes from SOURCE to DEST.  Does not null-terminate.
   In the public domain.
   By David MacKenzie <djm@gnu.ai.mit.edu>.  */

#if HAVE_CONFIG_H
#include <config.h>
#endif

void *
memmove (destaddr, sourceaddr, length)
     void *destaddr;
     const void *sourceaddr;
     unsigned length;
{
  char *dest = destaddr;
  const char *source = sourceaddr;
  if (source < dest)
    /* Moving from low mem to hi mem; start at end.  */
    for (source += length, dest += length; length; --length)
      *--dest = *--source;
  else if (source != dest)
    /* Moving from hi mem to low mem; start at beginning.  */
    for (; length; --length)
      *dest++ = *source++;

  return destaddr;
}
