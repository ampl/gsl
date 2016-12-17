REM @echo off

Echo GSL Windows Build NuGet
REM set INCREMENTDISABLE=true

REM # XEON x64 Build Vars #
set _SCRIPT_DRIVE=%~d0
set _SCRIPT_FOLDER=%~dp0
set SRC=%CD%
set BUILDTREE=%SRC%\build-win\
SET tbs_arch=x64
SET vcvar_arg=x86_amd64
SET ms_build_suffix=Bin\amd64
SET cmake_platform="Visual Studio 14 2015 Win64"
SET VS14="C:\Program Files (x86)\Microsoft Visual Studio 14.0\"
SET MSB14="C:\Program Files (x86)\MSBuild\14.0\"
SET MSBPath=%MSB14%%ms_build_suffix%
set PATH=%MSBPath%;%SRC%;%PATH%

REM # VC Vars #
call %VS14%\VC\vcvarsall %vcvar_arg%
@echo on

REM # Clean Build Tree #
if defined INCREMENTDISABLE ( 
	echo "Incremental Build disabled"
    rd /s /q %BUILDTREE%
) else (
	echo "Incremental Build enabled"
)
mkdir %BUILDTREE%
cd %BUILDTREE%

:shared_GSL
REM # GSL SHARED #
if exist %BUILDTREE%Release\gsl.lib (
    if exist %BUILDTREE%Release\gslcblas.lib (
		ECHO GSL Libs Found
		GOTO:copy_files
	)
)

ECHO %cmake_platform% SHARED
cmake -G %cmake_platform% ^
-DBUILD_SHARED_LIBS:BOOL=ON ^
-DCMAKE_CXX_FLAGS_RELEASE="/MD" ^
-DCMAKE_CXX_FLAGS_DEBUG="/MDd" ^
-DCMAKE_C_FLAGS_RELEASE="/MD" ^
-DCMAKE_C_FLAGS_DEBUG="/MDd" ^
-DCMAKE_BUILD_TYPE="Release" %SRC%
msbuild gsl.sln /p:Configuration=Release /m

:copy_files
set BINDIR=%SRC%\build-nuget\
echo %BINDIR%
xcopy %BUILDTREE%Release\gsl* %BINDIR%
xcopy %BUILDTREE%bin\Release\gsl.dll %BINDIR%
del %BINDIR%gsl
xcopy /I %BUILDTREE%gsl %BINDIR%gsl

:nuget_req
REM # make nuget packages from binaries #
nuget pack %SRC%\gsl-msvc14-x64.nuspec

GOTO:eof

REM --- exit ----
GOTO:eof
