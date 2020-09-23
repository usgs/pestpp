mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpc ..
make -j install
cpack -G TGZ
