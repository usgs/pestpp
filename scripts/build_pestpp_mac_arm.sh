#!/bin/sh
realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}
set -e

full_path=$(realpath $0)
script_path=$(dirname $full_path)
cd "$script_path"/..

rm -rfd build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -DFORCE_STATIC=OFF ..
make -j
cpack -G TGZ
cp *.tar.gz ../
