# A CMake script to run CMake in a Microsoft SDK build environment.

set(winsdk_key
  "HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Microsoft SDKs\\Windows")
find_program(WINSDK_SETENV NAMES SetEnv.cmd
  PATHS "[${winsdk_key};CurrentInstallFolder]/bin")
if (WINSDK_SETENV AND ARGS)
  if (NOT ARGS MATCHES Win64)
    set(setenv_arg "/x86")
  endif ()
  # If Microsoft SDK is installed create script run-msbuild.bat that
  # calls SetEnv.cmd to to set up build environment and runs msbuild.
  # It is useful when building Visual Studio projects with the SDK
  # toolchain rather than Visual Studio.
  # Set FrameworkPathOverride to get rid of MSB3644 warnings.
  file(WRITE "${CMAKE_BINARY_DIR}/run-msbuild.bat" "
    call \"${WINSDK_SETENV}\" ${setenv_arg}
    msbuild -p:FrameworkPathOverride=^\"C:\\Program Files^
\\Reference Assemblies\\Microsoft\\Framework\\.NETFramework\\v4.0^\" %*")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E echo "\"${WINSDK_SETENV}\" ${setenv_arg}")
endif ()
