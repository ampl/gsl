@echo on
rem Use this script if you want to use the Windows SDK toolchain.
rem It configures build environment and runs CMake forwarding all command-line
rem arguments to it, for example:
rem
rem   run-cmake -G "Visual Studio 10 Win64" .

set args=%*
set args=%args:"=\"%
for /F "delims=" %%i IN ('cmake "-DARGS=%args%" -P %~dp0/RunCMake.cmake') DO (
  call %%i
)
cmake %*
