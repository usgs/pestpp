#!/bin/bash
# This script runs within docker, using the command (which may require sudo):
#   docker run -t -i --rm -v $(dirname `pwd`):/io phusion/holy-build-box-64:latest bash /io/scripts/build_pestpp_linux_hbb.sh

set -e

# Activate Holy Build Box environment
source /hbb_exe/activate

#unset LIBRARY_PATH
unset LDFLAGS
unset CXXFLAGS

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/hbb_exe -DINSTALL_LOCAL=OFF /io
make -j
cpack -G TGZ

# Copy result to host, changing ownership to same as directory
FILE=$(ls *.tar.gz)
cp $FILE /io/
chown $(stat /io -c %u:%g) /io/$FILE

echo "See result $FILE"
