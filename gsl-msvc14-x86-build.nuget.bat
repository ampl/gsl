REM @echo off

Echo GSL Windows Build NuGet
set INCREMENTDISABLE=true

REM # XEON x86 Build Vars #
set _SCRIPT_DRIVE=%~d0
set _SCRIPT_FOLDER=%~dp0
set SRC=%CD%
set BUILDTREE=%SRC%\build-win\
SET tbs_arch=x86
SET vcvar_arg=x86
SET ms_build_suffix=Bin
SET cmake_platform="Visual Studio 14 2015"
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
rd /s /q %BINDIR%
mkdir %BINDIR%
echo %BINDIR%
xcopy %BUILDTREE%Release\gsl* %BINDIR%
xcopy %BUILDTREE%bin\Release\gsl.dll %BINDIR%
del %BINDIR%gsl
xcopy /I %BUILDTREE%gsl %BINDIR%gsl
copy %SRC%\gsl-msvc14.targets %BINDIR%\gsl-msvc14-x86.targets

mkdir %BUILDTREE%Static
cd %BUILDTREE%Static

:static_GSL
REM # GSL STATIC #
if exist %BUILDTREE%Static\Release\gsl.lib (
    if exist %BUILDTREE%Static\Release\gslcblas.lib (
		ECHO GSL Libs Found
		GOTO:copy_static_files
	)
)

ECHO %cmake_platform% STATIC
cmake -G %cmake_platform% ^
-DBUILD_SHARED_LIBS:BOOL=OFF ^
-DCMAKE_CXX_FLAGS_RELEASE="/MD" ^
-DCMAKE_CXX_FLAGS_DEBUG="/MDd" ^
-DCMAKE_C_FLAGS_RELEASE="/MD" ^
-DCMAKE_C_FLAGS_DEBUG="/MDd" ^
-DCMAKE_BUILD_TYPE="Release" %SRC%
msbuild gsl.sln /p:Configuration=Release /m

:copy_static_files
set BINDIR=%SRC%\build-nuget\static
rd /s /q %BINDIR%
mkdir %BINDIR%
echo %BINDIR%
xcopy %BUILDTREE%Static\Release\gsl* %BINDIR%

:nuget_req
REM # make nuget packages from binaries #
cd %BUILDTREE%
nuget pack %SRC%\gsl-msvc14-x86.nuspec
cd %SRC%
REM --- exit ----
GOTO:eof
