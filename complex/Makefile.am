noinst_LTLIBRARIES = libgslcomplex.la 

pkginclude_HEADERS = gsl_complex.h gsl_complex_math.h

INCLUDES = -I$(top_srcdir)

libgslcomplex_la_SOURCES = math.c inline.c

TESTS = $(check_PROGRAMS)
check_PROGRAMS = test

test_SOURCES = test.c results.h results1.h results2.h results_real.h results_zreal.h
test_LDADD = libgslcomplex.la ../err/libgslerr.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la  ../utils/libutils.la

results:
	emacs -batch -l test.el -f test1 | sed  's/00\+e/0e/g' > results1.h
	emacs -batch -l test.el -f test2 | sed  's/00\+e/0e/g' > results.h
	emacs -batch -l test.el -f test3 | sed  's/00\+e/0e/g' > results_real.h
	emacs -batch -l test.el -f test4 | sed  's/00\+e/0e/g' > results2.h
	emacs -batch -l test.el -f test5 | sed  's/00\+e/0e/g' > results_zreal.h
