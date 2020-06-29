@echo off

set first_path=%cd%
cd "%~dp0\.."

cd ..
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
rmdir /Q /S build
mkdir build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..
ninja
cpack -G ZIP
copy /y *.zip ..\

cd %first_path%
pause
