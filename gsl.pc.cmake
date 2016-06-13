prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=@CMAKE_INSTALL_PREFIX@
libdir=${prefix}/lib
includedir=${prefix}/include
GSL_CBLAS_LIB=-lgslcblas

Name: GSL
Description: GNU Scientific Library
Version: @PACKAGE_VERSION@
Libs: -L${libdir} -lgsl ${GSL_CBLAS_LIB} -lm -lm
Cflags: -I${includedir}
