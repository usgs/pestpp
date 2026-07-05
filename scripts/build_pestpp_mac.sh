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
# NB: macOS/Apple clang cannot produce fully static executables (no static
# libc / crt0.o), so FORCE_STATIC must stay OFF here - leaving it on fails the
# executable link with "ld: library not found for -lcrt0.o".
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ ..
make -j
cpack -G TGZ
cp *.tar.gz ../
