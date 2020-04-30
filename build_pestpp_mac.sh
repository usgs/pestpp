rm -r build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpc ..
make -j 5 
cpack
cp *.tar.gz ../