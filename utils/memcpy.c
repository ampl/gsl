/* Copy LEN bytes starting at SRCADDR to DESTADDR.  Result undefined
   if the source overlaps with the destination.
   Return DESTADDR. */

#if HAVE_CONFIG_H
#include <config.h>
#endif

void *
memcpy (destaddr, srcaddr, len)
     void *destaddr;
     const void *srcaddr;
     unsigned int len;
{
  char *dest = (char *) destaddr;
  const char *src = srcaddr;
  while (len-- > 0)
    *dest++ = *src++;
  return destaddr;
}
