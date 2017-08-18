.. index::
   single: autoconf, using with GSL

.. _chap_autoconf-macros:

***************
Autoconf Macros
***************

.. include:: include.rst

For applications using :code:`autoconf` the standard macro
:code:`AC_CHECK_LIB` can be used to link with GSL automatically
from a :code:`configure` script.  The library itself depends on the
presence of a |cblas| and math library as well, so these must also be
located before linking with the main :code:`libgsl` file.  The following
commands should be placed in the :file:`configure.ac` file to perform
these tests::

  AC_CHECK_LIB([m],[cos])
  AC_CHECK_LIB([gslcblas],[cblas_dgemm])
  AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

It is important to check for :code:`libm` and :code:`libgslcblas` before
:code:`libgsl`, otherwise the tests will fail.  Assuming the libraries
are found the output during the configure stage looks like this::

  checking for cos in -lm... yes
  checking for cblas_dgemm in -lgslcblas... yes
  checking for gsl_blas_dgemm in -lgsl... yes

If the library is found then the tests will define the macros
:code:`HAVE_LIBGSL`, :code:`HAVE_LIBGSLCBLAS`, :code:`HAVE_LIBM` and add
the options :code:`-lgsl -lgslcblas -lm` to the variable :code:`LIBS`.

The tests above will find any version of the library.  They are suitable
for general use, where the versions of the functions are not important.
An alternative macro is available in the file :file:`gsl.m4` to test for
a specific version of the library.  To use this macro simply add the
following line to your :file:`configure.in` file instead of the tests
above::

  AX_PATH_GSL(GSL_VERSION,
             [action-if-found],
             [action-if-not-found])

The argument :macro:`GSL_VERSION` should be the two or three digit
:code:`major.minor` or :code:`major.minor.micro` version number of the release
you require. A suitable choice for :code:`action-if-not-found` is::

  AC_MSG_ERROR(could not find required version of GSL)

Then you can add the variables :macro:`GSL_LIBS` and :macro:`GSL_CFLAGS` to
your Makefile.am files to obtain the correct compiler flags.
:macro:`GSL_LIBS` is equal to the output of the :code:`gsl-config --libs`
command and :macro:`GSL_CFLAGS` is equal to :code:`gsl-config --cflags`
command. For example::

  libfoo_la_LDFLAGS = -lfoo $(GSL_LIBS) -lgslcblas

Note that the macro :macro:`AX_PATH_GSL` needs to use the C compiler so it
should appear in the :file:`configure.in` file before the macro
:macro:`AC_LANG_CPLUSPLUS` for programs that use C++.

To test for :code:`inline` the following test should be placed in your
:file:`configure.in` file::

  AC_C_INLINE

  if test "$ac_cv_c_inline" != no ; then
    AC_DEFINE(HAVE_INLINE,1)
    AC_SUBST(HAVE_INLINE)
  fi

and the macro will then be defined in the compilation flags or by
including the file :file:`config.h` before any library headers.  

The following autoconf test will check for :code:`extern inline`::

  dnl Check for "extern inline", using a modified version
  dnl of the test for AC_C_INLINE from acspecific.mt
  dnl
  AC_CACHE_CHECK([for extern inline], ac_cv_c_extern_inline,
  [ac_cv_c_extern_inline=no
  AC_TRY_COMPILE([extern $ac_cv_c_inline double foo(double x);
  extern $ac_cv_c_inline double foo(double x) { return x+1.0; };
  double foo (double x) { return x + 1.0; };], 
  [  foo(1.0)  ],
  [ac_cv_c_extern_inline="yes"])
  ])

  if test "$ac_cv_c_extern_inline" != no ; then
    AC_DEFINE(HAVE_INLINE,1)
    AC_SUBST(HAVE_INLINE)
  fi

The substitution of portability functions can be made automatically if
you use :code:`autoconf`. For example, to test whether the BSD function
:func:`hypot` is available you can include the following line in the
configure file :file:`configure.in` for your application::

  AC_CHECK_FUNCS(hypot)

and place the following macro definitions in the file
:file:`config.h.in`::

  /* Substitute gsl_hypot for missing system hypot */

  #ifndef HAVE_HYPOT
  #define hypot gsl_hypot
  #endif

The application source files can then use the include command
:code:`#include <config.h>` to substitute :func:`gsl_hypot` for each
occurrence of :func:`hypot` when :func:`hypot` is not available.
