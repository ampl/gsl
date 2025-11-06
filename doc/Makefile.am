## Process this file with automake to produce Makefile.in

SUBDIRS = . examples

#info_TEXINFOS = gsl-ref.texi
noinst_TEXINFOS = gsl-design.texi

man_MANS = gsl.3 gsl-config.1 gsl-randist.1 gsl-histogram.1

figures =                              \
  images/bspline_deriv.png             \
  images/bspline_gram1.png             \
  images/bspline_gram2.png             \
  images/bspline_interp.png            \
  images/bspline_knots.png             \
  images/bspline_lsbreak.png           \
  images/bspline_lsend.png             \
  images/bspline_lsorder.png           \
  images/bspline_per.png               \
  images/bspline_proj.png              \
  images/cheb.png                      \
  images/dwt.png                       \
  images/fft-complex-radix2.png        \
  images/fft-real-mixedradix.png       \
  images/filt_edge.png                 \
  images/fit-exp.png                   \
  images/fit-wlinear2.png              \
  images/fit-wlinear.png               \
  images/gaussfilt.png                 \
  images/gaussfilt2.png                \
  images/histogram.png                 \
  images/histogram2d.png               \
  images/impulse.png                   \
  images/interp2d.png                  \
  images/interp_compare.png            \
  images/interp.png                    \
  images/interpp.png                   \
  images/min-interval.png              \
  images/movstat1.png                  \
  images/movstat2.png                  \
  images/movstat3.png                  \
  images/multilarge.png                \
  images/multimin.png                  \
  images/nlfit2b.png                   \
  images/nlfit2.png                    \
  images/nlfit3.png                    \
  images/ntuple.png                    \
  images/ode-vdp.png                   \
  images/qrng.png                      \
  images/rand-bernoulli.png            \
  images/rand-beta.png                 \
  images/rand-binomial.png             \
  images/rand-bivariate-gaussian.png   \
  images/rand-cauchy.png               \
  images/rand-chisq.png                \
  images/rand-erlang.png               \
  images/rand-exponential.png          \
  images/rand-exppow.png               \
  images/rand-fdist.png                \
  images/rand-flat.png                 \
  images/rand-gamma.png                \
  images/rand-gaussian.png             \
  images/rand-gaussian-tail.png        \
  images/rand-geometric.png            \
  images/rand-gumbel1.png              \
  images/rand-gumbel2.png              \
  images/rand-hypergeometric.png       \
  images/rand-landau.png               \
  images/rand-laplace.png              \
  images/rand-levy.png                 \
  images/rand-levyskew.png             \
  images/rand-logarithmic.png          \
  images/rand-logistic.png             \
  images/rand-lognormal.png            \
  images/rand-nbinomial.png            \
  images/random-walk.png               \
  images/rand-pareto.png               \
  images/rand-pascal.png               \
  images/rand-poisson.png              \
  images/rand-rayleigh.png             \
  images/rand-rayleigh-tail.png        \
  images/rand-tdist.png                \
  images/rand-weibull.png              \
  images/regularized2.png              \
  images/regularized.png               \
  images/robust.png                    \
  images/siman-12-cities.png           \
  images/siman-energy.png              \
  images/siman-final-route.png         \
  images/siman-initial-route.png       \
  images/siman-test.png                \
  images/sparse_poisson.png

rst_files =                            \
  autoconf.rst                         \
  blas.rst                             \
  bspline.rst                          \
  cblas.rst                            \
  cheb.rst                             \
  combination.rst                      \
  complex.rst                          \
  const.rst                            \
  contrib.rst                          \
  debug.rst                            \
  dht.rst                              \
  diff.rst                             \
  dwt.rst                              \
  eigen.rst                            \
  err.rst                              \
  fdl.rst                              \
  fft.rst                              \
  filter.rst                           \
  gpl.rst                              \
  histogram.rst                        \
  ieee754.rst                          \
  include.rst                          \
  index.rst                            \
  integration.rst                      \
  interp.rst                           \
  intro.rst                            \
  linalg.rst                           \
  lls.rst                              \
  math.rst                             \
  min.rst                              \
  montecarlo.rst                       \
  movstat.rst                          \
  multimin.rst                         \
  multiroots.rst                       \
  multiset.rst                         \
  nls.rst                              \
  ntuple.rst                           \
  ode-initval.rst                      \
  permutation.rst                      \
  poly.rst                             \
  qrng.rst                             \
  randist.rst                          \
  rng.rst                              \
  roots.rst                            \
  rstat.rst                            \
  siman.rst                            \
  sort.rst                             \
  spblas.rst                           \
  specfunc-airy.rst                    \
  specfunc-bessel.rst                  \
  specfunc-clausen.rst                 \
  specfunc-coulomb.rst                 \
  specfunc-coupling.rst                \
  specfunc-dawson.rst                  \
  specfunc-debye.rst                   \
  specfunc-dilog.rst                   \
  specfunc-elementary.rst              \
  specfunc-ellint.rst                  \
  specfunc-elljac.rst                  \
  specfunc-erf.rst                     \
  specfunc-expint.rst                  \
  specfunc-exp.rst                     \
  specfunc-fermi-dirac.rst             \
  specfunc-gamma.rst                   \
  specfunc-gegenbauer.rst              \
  specfunc-hermite.rst                 \
  specfunc-hyperg.rst                  \
  specfunc-laguerre.rst                \
  specfunc-lambert.rst                 \
  specfunc-legendre.rst                \
  specfunc-log.rst                     \
  specfunc-mathieu.rst                 \
  specfunc-pow-int.rst                 \
  specfunc-psi.rst                     \
  specfunc.rst                         \
  specfunc-synchrotron.rst             \
  specfunc-transport.rst               \
  specfunc-trig.rst                    \
  specfunc-zeta.rst                    \
  splinalg.rst                         \
  spmatrix.rst                         \
  statistics.rst                       \
  sum.rst                              \
  usage.rst                            \
  vectors.rst

static_files =                         \
  _static/fdl.txt                      \
	_static/gpl.txt

noinst_DATA = $(figures) $(rst_files) $(static_files)

EXTRA_DIST = $(man_MANS) $(noinst_DATA) gsl-design.texi gsl-ref.info conf.py

SPHINX_OPTS        =
SPHINX_BUILD       = sphinx-build
SPHINX_SOURCEDIR   = .
SPHINX_BUILDDIR    = _build

# targets supported by sphinx
SPHINX_TARGETS = html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf latexpdfja text man texinfo gettext changes xml pseudoxml linkcheck doctest coverage

$(SPHINX_TARGETS) :
	@$(SPHINX_BUILD) -M $@ "$(SPHINX_SOURCEDIR)" "$(SPHINX_BUILDDIR)" $(SPHINX_OPTS) $(O)

# special case info target
info:
	@$(SPHINX_BUILD) -M $@ "$(SPHINX_SOURCEDIR)" "$(SPHINX_BUILDDIR)" $(SPHINX_OPTS) $(O)
	cp "$(SPHINX_BUILDDIR)/texinfo/gsl-ref.info" .

install-data-local: $(srcdir)/gsl-ref.info
	$(MKDIR_P) "$(DESTDIR)$(infodir)"
	$(INSTALL_DATA) $(srcdir)/gsl-ref.info $(DESTDIR)$(infodir)

uninstall-local:
	rm -f "$(DESTDIR)$(infodir)/gsl-ref.info"

clean-local:
	rm -rf "$(SPHINX_BUILDDIR)"
