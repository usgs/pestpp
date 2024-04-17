@echo off

set first_path=%cd%

cd "%~dp0\.."

rem fast but without fortran
rmdir /Q /S bin
rmdir /Q /S build
mkdir build
rem call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\bin\compilervars.bat" intel64
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 --force
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icx ..
ninja
cpack -G ZIP
copy /y *.zip ..\

cd ..
rem call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
rmdir /Q /S build
mkdir build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..
ninja
cpack -G ZIP
copy /y *.zip ..\

cd %first_path%
pause
