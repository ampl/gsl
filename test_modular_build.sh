#!/bin/bash

dirs="cblas const  err        ieee-utils     monte               multimin      permutation  rstat     spmatrix    utils
blas     cdf          deriv  fft        integration    movstat             multiroots    poly         siman     statistics  vector
block    cheb         dht    filter     interpolation  multifit            multiset      qrng         sort      sum         wavelet
bspline  diff   fit        linalg         multifit_nlinear    ntuple        randist      spblas    
bst      combination     gsl        matrix         multilarge          ode-initval   rng          specfunc  sys
complex      eigen  histogram  min            multilarge_nlinear  ode-initval2  roots        splinalg "


cd build/wsl

for t in ${dirs}; do
  cmake ../.. -DBUILDLIBS=$t -DNO_AMPL_BINDINGS=1 || exit 1 
  make -j16 || exit 1
done

cd ../..