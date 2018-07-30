/* statistics/test_int_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jim Davies, Brian Gough
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

void FUNCTION (test, func) (const size_t stridea, const size_t strideb);

void
FUNCTION (test, func) (const size_t stridea, const size_t strideb)
{
  /* sample sets of integers */
  size_t i;
  const size_t ina = 20, inb = 20;

  const BASE raw1[] = {1, 2, 3, 4, 5, 6} ;
  
  const BASE irawa[] =
  {17, 18, 16, 18, 12,
   20, 18, 20, 20, 22,
   20, 10, 8, 12, 16,
   16, 18, 20, 18, 21};

  const BASE irawb[] =
  {19, 20, 22, 24, 10,
   25, 20, 22, 21, 23,
   20, 10, 12, 14, 12,
   20, 22, 24, 23, 17};

  BASE * sorted ;

  BASE * test1 = (BASE *) malloc (stridea * 6 * sizeof(BASE));
  BASE * igroupa = (BASE *) malloc (stridea * ina * sizeof(BASE));
  BASE * igroupb = (BASE *) malloc (strideb * inb * sizeof(BASE));

  double rel = 1e-10 ;

  for (i = 0 ; i < ina ; i++)
    igroupa[i * stridea] = irawa[i] ;

  for (i = 0 ; i < inb ; i++)
    igroupb[i * strideb] = irawb[i] ;

  for (i = 0 ; i < 6 ; i++)
    test1[i * stridea] = raw1[i] ;



  {
    double mean = FUNCTION(gsl_stats,mean) (igroupa, stridea, ina);
    double expected = 17.0;
    gsl_test_rel (mean,expected, rel, NAME(gsl_stats) "_mean (integer)");
  }

  {
    double mean = FUNCTION(gsl_stats,mean) (test1, stridea, 6);
    double expected = 3.5;
    gsl_test_rel (mean,expected, rel, NAME(gsl_stats) "_mean (fractional)");
  }

  {
    double mean = FUNCTION(gsl_stats,mean) (igroupa, stridea, ina);
    double var = FUNCTION(gsl_stats,variance_with_fixed_mean) (igroupa, stridea, ina, mean);
    double expected = 13.7;
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_variance_with_fixed_mean");
  }

  {
    double mean = FUNCTION(gsl_stats,mean) (igroupa, stridea, ina);
    double sd = FUNCTION(gsl_stats,sd_with_fixed_mean) (igroupa, stridea, ina, mean);
    double expected = 3.70135110466435;
    gsl_test_rel (sd, expected, rel, NAME(gsl_stats) "_sd_with_fixed_mean");
  }

  {
    double var = FUNCTION(gsl_stats,variance) (igroupa, stridea, ina);
    double expected = 14.4210526315789;
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_variance");
  }

  {
    double sd_est = FUNCTION(gsl_stats,sd) (igroupa, stridea, ina);
    double expected = 3.79750610685209;
    gsl_test_rel (sd_est, expected, rel, NAME(gsl_stats) "_sd");
  }

  {
    double absdev = FUNCTION(gsl_stats,absdev) (igroupa, stridea, ina);
    double expected = 2.9;
    gsl_test_rel (absdev, expected, rel, NAME(gsl_stats) "_absdev");
  }

  {
    double skew = FUNCTION(gsl_stats,skew) (igroupa, stridea, ina);
    double expected = -0.909355923168064;
    gsl_test_rel (skew, expected, rel, NAME(gsl_stats) "_skew");
  }

  {
    double kurt = FUNCTION(gsl_stats,kurtosis) (igroupa, stridea, ina);
    double expected = -0.233692524908094 ;
    gsl_test_rel (kurt, expected, rel, NAME(gsl_stats) "_kurtosis");
  }

  {
    double c = FUNCTION(gsl_stats,covariance) (igroupa, stridea, igroupb, strideb, inb);
    double expected = 14.5263157894737;
    gsl_test_rel (c, expected, rel, NAME(gsl_stats) "_covariance");
  }

  {
    double r = FUNCTION(gsl_stats,correlation) (igroupa, stridea, igroupb, strideb, inb);
    double expected = 0.793090350710101;
    gsl_test_rel (r, expected, rel, NAME(gsl_stats) "_correlation");
  }

  {
    double *work = malloc(2 * ina * sizeof(double));
    double r = FUNCTION(gsl_stats,spearman) (igroupa, stridea, igroupb, strideb, inb, work);
    double expected = 0.644239670474577;
    gsl_test_rel (r, expected, rel, NAME(gsl_stats) "_spearman");
    free(work);
  }

  {
    double pv = FUNCTION(gsl_stats,pvariance) (igroupa, stridea, ina, igroupb, strideb, inb);
    double expected = 18.8421052631579;
    gsl_test_rel (pv, expected, rel, NAME(gsl_stats) "_pvariance");
  }

  {
    double t = FUNCTION(gsl_stats,ttest) (igroupa, stridea, ina, igroupb, strideb, inb);
    double expected = -1.45701922702927;
    gsl_test_rel (t, expected, rel, NAME(gsl_stats) "_ttest");
  }

  {
    int max = FUNCTION(gsl_stats,max) (igroupa, stridea, ina);
    int expected = 22;
    gsl_test (max != expected,
              NAME(gsl_stats) "_max (%d observed vs %d expected)", max, expected);
  }

  {
    int min = FUNCTION(gsl_stats,min) (igroupa, stridea, ina);
    int expected = 8;
    gsl_test (min != expected,
              NAME(gsl_stats) "_min (%d observed vs %d expected)", min, expected);
  }

  {
    BASE min, max;
    BASE expected_max = 22;
    BASE expected_min = 8;
    
    FUNCTION(gsl_stats,minmax) (&min, &max, igroupa, stridea, ina);
 
    gsl_test  (max != expected_max,
               NAME(gsl_stats) "_minmax max (" OUT_FORMAT " observed vs " OUT_FORMAT " expected)", 
               max, expected_max);
    gsl_test  (min != expected_min,
               NAME(gsl_stats) "_minmax min (" OUT_FORMAT " observed vs " OUT_FORMAT " expected)", 
               min, expected_min);
  }

  {
    int max_index = FUNCTION(gsl_stats,max_index) (igroupa, stridea, ina);
    int expected = 9 ;
    gsl_test (max_index != expected,
              NAME(gsl_stats) "_max_index (%d observed vs %d expected)",
              max_index, expected);
  }

  {
    int min_index = FUNCTION(gsl_stats,min_index) (igroupa, stridea, ina);
    int expected = 12 ;
    gsl_test (min_index != expected,
              NAME(gsl_stats) "_min_index (%d observed vs %d expected)",
              min_index, expected);
  }

  {
    size_t min_index, max_index;
    size_t expected_max_index = 9;
    size_t expected_min_index = 12;

    FUNCTION(gsl_stats,minmax_index) (&min_index, &max_index, igroupa, stridea, ina);

    gsl_test  (max_index != expected_max_index,
               NAME(gsl_stats) "_minmax_index max (%u observed vs %u expected)", 
               max_index, expected_max_index);
    gsl_test  (min_index != expected_min_index,
               NAME(gsl_stats) "_minmax_index min (%u observed vs %u expected)", 
               min_index, expected_min_index);
  }


  sorted = (BASE *) malloc(stridea * ina * sizeof(BASE)) ;

  for (i = 0 ; i < ina ; i++)
    sorted[stridea * i] = igroupa[stridea * i] ;


  TYPE(gsl_sort)(sorted, stridea, ina) ;

  {
    double median = FUNCTION(gsl_stats,median_from_sorted_data)(sorted, stridea, ina) ;
    double expected = 18;
    gsl_test_rel (median,expected, rel,
                  NAME(gsl_stats) "_median_from_sorted_data (even)");
  }

  {
    double median = FUNCTION(gsl_stats,median_from_sorted_data)(sorted, stridea, ina - 1) ;
    double expected = 18;
    gsl_test_rel (median,expected, rel,
                  NAME(gsl_stats) "_median_from_sorted_data (odd)");
  }

  {
    BASE * work = (BASE *) malloc (stridea * ina * sizeof(BASE));
    double expected = 18;
    double median;

    for (i = 0; i < ina; i++)
      work[i * stridea] = (BASE) irawa[i];

    median = FUNCTION(gsl_stats,median)(work, stridea, ina) ;

    gsl_test_rel (median, expected, rel,
                  NAME(gsl_stats) "_median (even)");

    free(work);
  }

  {
    BASE * work = (BASE *) malloc (stridea * (ina - 1) * sizeof(BASE));
    double expected = 18;
    double median;

    for (i = 0; i < ina - 1; i++)
      work[i * stridea] = sorted[i * stridea];

    median = FUNCTION(gsl_stats,median)(work, stridea, ina - 1) ;

    gsl_test_rel (median, expected, rel,
                  NAME(gsl_stats) "_median (odd)");

    free(work);
  }

  {
    const double trim = 0.31;
    double trmean = FUNCTION(gsl_stats,trmean_from_sorted_data)(trim, sorted, stridea, ina) ;
    double expected = 17.875;
    gsl_test_rel  (trmean, expected, rel,
                   NAME(gsl_stats) "_trmean_from_sorted_data (even)");
  }

  {
    const double trim = 0.3;
    double trmean = FUNCTION(gsl_stats,trmean_from_sorted_data)(trim, sorted, stridea, ina - 1) ;
    double expected = 17.66666666666666;
    gsl_test_rel  (trmean, expected, rel,
                   NAME(gsl_stats) "_trmean_from_sorted_data (odd)");
  }

  {
    double gastwirth = FUNCTION(gsl_stats,gastwirth_from_sorted_data)(sorted, stridea, ina) ;
    double expected = 17.9;
    gsl_test_rel  (gastwirth, expected, rel,
                   NAME(gsl_stats) "_gastwirth_from_sorted_data (even)");
  }

  {
    double gastwirth = FUNCTION(gsl_stats,gastwirth_from_sorted_data)(sorted, stridea, ina - 1) ;
    double expected = 17.4;
    gsl_test_rel  (gastwirth, expected, rel,
                   NAME(gsl_stats) "_gastwirth_from_sorted_data (odd)");
  }

  {
    double zeroth = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, ina, 0.0) ;
    double expected = 8;
    gsl_test_rel (zeroth,expected, rel,
                  NAME(gsl_stats) "_quantile_from_sorted_data (0)");
  }

  {
    double top = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, ina, 1.0) ;
    double expected = 22;
    gsl_test_rel (top,expected, rel,
                  NAME(gsl_stats) "_quantile_from_sorted_data (100)");
  }

  {
    double median = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, ina, 0.5) ;
    double expected = 18;
    gsl_test_rel (median,expected, rel,
                  NAME(gsl_stats) "_quantile_from_sorted_data (50, even)");
  }

  {
    double median = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, ina - 1, 0.5);
    double expected = 18;
    gsl_test_rel (median,expected, rel,
                  NAME(gsl_stats) "_quantile_from_sorted_data (50, odd)");
  }

  {
    size_t k;
    BASE * work = (BASE *) malloc (stridea * ina * sizeof(BASE));

    for (k = 0; k < ina; ++k)
      {
        int expected = sorted[k * stridea];
        int kselect;

        /* copy irawa[] for each k, since gsl_stats_select() changes input array */
        for (i = 0; i < ina; i++)
          work[i * stridea] = (BASE) irawa[i];

        kselect = FUNCTION(gsl_stats,select)(work, stridea, ina, k);

        gsl_test (kselect != expected,
                  NAME(gsl_stats) "_select (%d observed vs %d expected)",
                  kselect, expected);
      }

    free(work);
  }

  {
    double * work = (double *) malloc (ina * sizeof(double));
    double expected = 2.0;
    double mad0 = FUNCTION(gsl_stats,mad0)(igroupa, stridea, ina, work);
    gsl_test_rel (mad0, expected, rel, NAME(gsl_stats) "_mad0 (even)");
    free(work);
  }

  {
    double * work = (double *) malloc ((ina - 1) * sizeof(double));
    double expected = 2.0;
    double mad0 = FUNCTION(gsl_stats,mad0)(igroupa, stridea, ina - 1, work);
    gsl_test_rel (mad0, expected, rel, NAME(gsl_stats) "_mad0 (odd)");
    free(work);
  }

  {
    BASE * work = malloc(ina * sizeof(BASE));
    double sn = FUNCTION(gsl_stats,Sn_from_sorted_data)(sorted, stridea, ina, work) ;
    double expected = 2.3852;
    gsl_test_rel (sn, expected, rel,
                  NAME(gsl_stats) "_Sn_from_sorted_data (even)");
    free(work);
  }

  {
    BASE * work = malloc((ina - 1) * sizeof(BASE));
    double sn = FUNCTION(gsl_stats,Sn_from_sorted_data)(sorted, stridea, ina - 1, work) ;
    double expected = 2.503801104972376;
    gsl_test_rel (sn, expected, rel,
                  NAME(gsl_stats) "_Sn_from_sorted_data (odd)");
    free(work);
  }

  {
    BASE * work = malloc(3 * ina * sizeof(BASE));
    int * work_int = malloc(5 * ina * sizeof(int));
    double qn = FUNCTION(gsl_stats,Qn_from_sorted_data)(sorted, stridea, ina, work, work_int) ;
    double expected = 3.732513488036874;
    gsl_test_rel (qn, expected, rel,
                  NAME(gsl_stats) "_Qn_from_sorted_data (even)");
    free(work);
    free(work_int);
  }

  {
    BASE * work = malloc(3 * (ina - 1) * sizeof(BASE));
    int * work_int = malloc(5 * (ina - 1) * sizeof(int));
    double qn = FUNCTION(gsl_stats,Qn_from_sorted_data)(sorted, stridea, ina - 1, work, work_int) ;
    double expected = 4.118443402621429;
    gsl_test_rel (qn, expected, rel,
                  NAME(gsl_stats) "_Qn_from_sorted_data (odd)");
    free(work);
    free(work_int);
  }

  free (sorted);
  free (igroupa);
  free (igroupb);
  free (test1);
}
