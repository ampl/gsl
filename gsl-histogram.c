/* histogram/gsl-histogram.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <gsl/gsl_histogram.h>

void
print_help(void)
{
  fprintf (stderr, "Usage: gsl-histogram [-u] xmin xmax [n]\n");
  fprintf (stderr, "Computes a histogram of the data on stdin using n bins from xmin to xmax.\n"
                   "If n is unspecified then bins of integer width are used.\n"
                   "If -u is given, histogram is normalized so that sum of all bins is unity.\n");
}

int
main (int argc, char **argv)
{
  double a = 0.0, b = 1.0;
  size_t n = 10;
  int unit = 0;
  int c;

  while ((c = getopt(argc, argv, "u")) != (-1))
    {
      switch (c)
        {
          case 'u':
            unit = 1;
            break;

          default:
            print_help();
            exit(0);
            break;
        }
    }

  if (argc - optind < 2)
    {
      print_help();
      exit(0);
    }

  a = atof (argv[optind++]);
  b = atof (argv[optind++]);

  if (argc - optind > 0)
    n = atoi (argv[optind++]);
  else
    n = (int)(b - a) ;

  {
    double x;
    gsl_histogram *h = gsl_histogram_alloc (n);

    gsl_histogram_set_ranges_uniform (h, a, b);

    while (fscanf(stdin, "%lg", &x) == 1)
      {
        gsl_histogram_increment(h, x);
      }

#ifdef DISPLAY_STATS
    {
      double mean = gsl_histogram_mean (h);
      double sigma = gsl_histogram_sigma (h);
      fprintf (stdout, "# mean = %g\n", mean);
      fprintf (stdout, "# sigma = %g\n", sigma);
    }
#endif

    /* normalize histogram if needed */
    if (unit)
      {
        double sum = gsl_histogram_sum(h);

        if (sum > 0.0)
          gsl_histogram_scale(h, 1.0 / sum);
      }

    gsl_histogram_fprintf (stdout, h, "%g", "%g");

    gsl_histogram_free (h);
  }

  return 0;
}
