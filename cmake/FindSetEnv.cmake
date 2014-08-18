# A CMake script to find SetEnv.cmd.

set(winsdk_key
  "HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Microsoft SDKs\\Windows")
find_program(WINSDK_SETENV NAMES SetEnv.cmd
  PATHS "[${winsdk_key};CurrentInstallFolder]/bin")
if (WINSDK_SETENV AND ARGS)
  if (NOT ARGS MATCHES Win64)
    set(setenv_arg "/x86")
  endif ()
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E echo "\"${WINSDK_SETENV}\" ${setenv_arg}")
endif ()
