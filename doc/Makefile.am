## Process this file with automake to produce Makefile.in

SUBDIRS = . examples

info_TEXINFOS = gsl-ref.texi
noinst_TEXINFOS = gsl-design.texi
gsl_ref_TEXINFOS = err.texi cblas.texi blas.texi min.texi fft.texi rng.texi randist.texi roots.texi rstat.texi statistics.texi specfunc.texi specfunc-airy.texi specfunc-bessel.texi specfunc-clausen.texi specfunc-coulomb.texi specfunc-coupling.texi specfunc-dawson.texi specfunc-debye.texi specfunc-dilog.texi specfunc-elementary.texi specfunc-ellint.texi specfunc-elljac.texi specfunc-erf.texi specfunc-exp.texi specfunc-expint.texi specfunc-fermi-dirac.texi specfunc-gamma.texi specfunc-gegenbauer.texi specfunc-hyperg.texi specfunc-lambert.texi specfunc-laguerre.texi specfunc-legendre.texi specfunc-log.texi specfunc-mathieu.texi specfunc-pow-int.texi specfunc-psi.texi specfunc-synchrotron.texi specfunc-transport.texi specfunc-trig.texi specfunc-zeta.texi siman.texi vectors.texi debug.texi histogram.texi ode-initval.texi integration.texi ieee754.texi montecarlo.texi sum.texi intro.texi usage.texi dwt.texi dht.texi interp.texi poly.texi linalg.texi eigen.texi multiroots.texi sort.texi permutation.texi combination.texi multiset.texi complex.texi math.texi fitting.texi multifit.texi const.texi ntuple.texi diff.texi qrng.texi cheb.texi multimin.texi gpl.texi fdl.texi autoconf.texi bspline.texi spblas.texi spmatrix.texi splinalg.texi

man_MANS = gsl.3 gsl-config.1 gsl-randist.1 gsl-histogram.1

figures = multimin.eps siman-test.eps siman-energy.eps 12-cities.eps initial-route.eps final-route.eps fft-complex-radix2-f.eps fft-complex-radix2-t.eps fft-complex-radix2.eps fft-real-mixedradix.eps roots-bisection.eps roots-false-position.eps roots-newtons-method.eps roots-secant-method.eps histogram.eps histogram2d.eps min-interval.eps fit-wlinear.eps fit-wlinear2.eps fit-exp.eps ntuple.eps qrng.eps cheb.eps vdp.eps interp2.eps interp2d.eps rand-beta.tex rand-binomial.tex rand-cauchy.tex rand-chisq.tex rand-erlang.tex rand-exponential.tex rand-fdist.tex rand-flat.tex rand-gamma.tex rand-gaussian.tex rand-geometric.tex rand-laplace.tex rand-logarithmic.tex rand-logistic.tex rand-lognormal.tex rand-pareto.tex rand-poisson.tex rand-hypergeometric.tex rand-nbinomial.tex rand-pascal.tex rand-bivariate-gaussian.tex rand-rayleigh.tex rand-rayleigh-tail.tex rand-tdist.tex rand-weibull.tex random-walk.tex randplots.gnp rand-exppow.tex rand-levy.tex rand-levyskew.tex rand-gumbel.tex rand-bernoulli.tex rand-gaussian-tail.tex rand-gumbel1.tex rand-gumbel2.tex landau.dat rand-landau.tex dwt-orig.eps dwt-samp.eps interpp2.eps bspline.eps robust.eps sparse_poisson.eps interp_compare.eps lcurve.eps multilarge.eps

noinst_DATA = $(figures)

EXTRA_DIST = $(man_MANS) $(noinst_DATA) gsl-design.texi fftalgorithms.tex fftalgorithms.bib algorithm.sty algorithmic.sty calc.sty statnotes.tex

# pdf disabled, use postscript and ps2pdf
.PHONY: pdf
pdf:
