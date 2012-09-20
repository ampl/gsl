file(GLOB GSL_HEADER_PATHS "${GSL_SOURCE_DIR}/gsl*.h"
  "${GSL_SOURCE_DIR}/*/gsl*.h")

foreach (path ${GSL_HEADER_PATHS})
  get_filename_component(name ${path} NAME)
  set(GSL_HEADERS ${GSL_HEADERS} ${name})
endforeach ()

configure_file(copy-headers.cmake.in copy-headers.cmake @ONLY)
add_custom_command(OUTPUT ${GSL_HEADERS}
  COMMAND ${CMAKE_COMMAND} -P copy-headers.cmake)
add_custom_target(copy-headers ALL DEPENDS ${GSL_HEADERS})
