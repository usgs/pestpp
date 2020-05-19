#!/bin/sh

set -e

full_path=$(realpath $0)
script_path=$(dirname $full_path)
cd "$script_path"/..

rm -r build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpc ..
make -j
cpack -G TGZ
cp *.tar.gz ../
