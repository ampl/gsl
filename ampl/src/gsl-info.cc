/*
 This program generates AMPL declarations for the functions provided
 by the amplgsl library and prints the library version information.

 Copyright (C) 2020 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <stdio.h>
#include <string.h>
#include "funcadd.h"

#undef fopen
#undef fclose
#undef fprintf
#undef printf

static FILE *out;

#define UNUSED(x) (void)(x)


/* See AddFunc in funcadd.h */
static void declare_func(const char *name, rfunc f,
    int type, int nargs, void *funcinfo, AmplExports *ae) {
  const char *attr = "";
  UNUSED(f);
  UNUSED(type);
  UNUSED(nargs);
  UNUSED(funcinfo);
  UNUSED(ae);
  if ((type & FUNCADD_RANDOM_VALUED) != 0)
    attr = " random";
  else if ((type & FUNCADD_STRING_VALUED) != 0)
    attr = " symbolic";
  if (strcmp(name, "gsl_version") == 0) {
    typedef const char *(*Func)(arglist *al);
    Func get_version = (Func)f;
    printf("AMPLGSL version %d, GSL version %s\n", LIBDATE, get_version(NULL));
    fprintf(out, "# AMPLGSL version %d, GSL version %s\nload amplgsl.dll;\n", LIBDATE, get_version(NULL));
  }
  if (strcmp(name, "gsl_sort") ==0)
    return;
  fprintf(out, "function %s%s;\n", name, attr);
}

static void dummy_at_reset(AmplExports *ae, Exitfunc *f, void *data) {
  UNUSED(ae);
  UNUSED(f);
  UNUSED(data);
}

int main() {
  AmplExports ae = AmplExports();
  ae.Addfunc = declare_func;
  ae.AtReset = dummy_at_reset;
  out = fopen("gsl.ampl.tmp", "w");
  if (!out) {
    printf("Can't open gsl.ampl.tmp");
    return 1;
  }
  fprintf(out,
      "# Automatically generated AMPL declarations for the GSL functions.\n");
  funcadd_ASL(&ae);
  fclose(out);
  remove("gsl.ampl");
  if (rename("gsl.ampl.tmp", "gsl.ampl")) {
    printf("Can't rename gsl.ampl.tmp to gsl.ampl");
    return 1;
  }
  return 0;
}
