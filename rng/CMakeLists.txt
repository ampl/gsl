add_library(gslrng
  borosh13.c cmrg.c coveyou.c default.c file.c fishman18.c fishman20.c
  fishman2x.c gfsr4.c knuthran2.c knuthran.c knuthran2002.c lecuyer21.c
  minstd.c mrg.c mt.c r250.c ran0.c ran1.c ran2.c ran3.c rand48.c rand.c
  random.c randu.c ranf.c ranlux.c ranlxd.c ranlxs.c ranmar.c rng.c slatec.c
  taus.c taus113.c transputer.c tt.c types.c uni32.c uni.c vax.c waterman14.c
  zuf.c inline.c)

add_executable(gslrng-test test.c)
target_link_libraries(gslrng-test
  gslrng gslieeeutils gslerr gsltest gslsys utils ${CMAKE_REQUIRED_LIBRARIES})
