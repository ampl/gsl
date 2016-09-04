# CMake initialization code that should be run before the project command.

include(CMakeParseArguments)

# Joins arguments and sets the result to <var>.
# Usage:
#   join(<var> [<arg>...])
function (join var)
  unset(result)
  foreach (arg ${ARGN})
    if (DEFINED result)
      set(result "${result} ${arg}")
    else ()
      set(result "${arg}")
    endif ()
  endforeach ()
  set(${var} "${result}" PARENT_SCOPE)
endfunction ()

# Sets cache variable <var> to the value <value>. The arguments
# following <type> are joined into a single docstring which allows
# breaking long documentation into smaller strings.
# Usage:
#   set_cache(<var> <value> <type> docstring... [FORCE])
function (set_cache var value type)
  cmake_parse_arguments(set_cache FORCE "" "" ${ARGN})
  unset(force)
  if (set_cache_FORCE)
    set(force FORCE)
  endif ()
  join(docstring ${set_cache_UNPARSED_ARGUMENTS})
  set(${var} ${value} CACHE ${type} "${docstring}" ${force})
endfunction ()

if (NOT CMAKE_BUILD_TYPE)
  # Set the default CMAKE_BUILD_TYPE to Release.
  # This should be done before the project command since the latter sets
  # CMAKE_BUILD_TYPE itself.
  set_cache(CMAKE_BUILD_TYPE Release STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or"
    "CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif ()

function (override var file)
  if (EXISTS "${file}")
    set(${var} ${file} PARENT_SCOPE)
  endif ()
endfunction ()

# Set compiler options (MSVC)
# This should be done before the project command.
override(CMAKE_USER_MAKE_RULES_OVERRIDE
  ${CMAKE_CURRENT_LIST_DIR}/c_flag_overrides.cmake)
override(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
  ${CMAKE_CURRENT_LIST_DIR}/cxx_flag_overrides.cmake)

if (MP_WINSDK)
  # Find Windows SDK.
  set(winsdk_key
    "HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Microsoft SDKs\\Windows")
  find_program(MP_SETENV NAMES SetEnv.cmd
    PATHS "[${winsdk_key};CurrentInstallFolder]/bin")
  if (MP_SETENV)
    if (NOT CMAKE_GENERATOR MATCHES Win64)
      set(setenv_arg "/x86")
    endif ()

    # Call SetEnv.cmd and set environment variables accordingly.
    message(STATUS "Found SetEnv: ${MP_SETENV}")
    set(run_setenv "CMakeFiles\\run-setenv.bat")
    file(WRITE ${run_setenv} "call %*\nset\n")
    execute_process(COMMAND ${run_setenv} "${MP_SETENV}" "${setenv_arg}"
      OUTPUT_VARIABLE out ERROR_VARIABLE err)
    string(REPLACE ";" "\;" out "${out}")
    string(REGEX MATCHALL "[^\n]+\n" out "${out}")
    foreach (env ${out})
      if (env MATCHES "([^=]+)=(.*)\n")
        set(ENV{${CMAKE_MATCH_1}} "${CMAKE_MATCH_2}")
      endif ()
    endforeach ()

    # If Microsoft SDK is installed create script run-msbuild.bat that
    # calls SetEnv.cmd to set up build environment and runs msbuild.
    # It is useful when building Visual Studio projects with the SDK
    # toolchain rather than Visual Studio.
    # Set FrameworkPathOverride to get rid of MSB3644 warnings.
    file(WRITE "${CMAKE_BINARY_DIR}/run-msbuild.bat" "
      call \"${MP_SETENV}\" ${setenv_arg}
      msbuild -p:FrameworkPathOverride=^\"C:\\Program Files^
\\Reference Assemblies\\Microsoft\\Framework\\.NETFramework\\v4.0^\" %*")
  endif ()
endif ()
