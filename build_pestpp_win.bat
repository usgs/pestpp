rmdir /Q /S build
mkdir build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icl ..
ninja
cpack -G ZIP
copy /y *.zip ..\

cd ..
rmdir /Q /S build
mkdir build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..
ninja
cpack -G ZIP
copy /y *.zip ..\
cd ..