#!/bin/bash
# This script runs within docker, using the command (which may require sudo):
#   docker run -t -i --rm -v $(dirname `pwd`):/io phusion/holy-build-box-64:latest bash /io/scripts/build_pestpp_linux_hbb.sh

set -e

# Activate Holy Build Box environment
hbb_prefix=/hbb_exe
source /$hbb_prefix/activate

# Dependencies
yum install -y lapack-devel

#unset LIBRARY_PATH
unset LDFLAGS
unset CXXFLAGS

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_Fortran=ON -DCMAKE_INSTALL_PREFIX=$hbb_prefix -DINSTALL_LOCAL=OFF /io
make -j
cpack -G TGZ

# Copy result to host, changing ownership to same as directory
FILE=$(ls *.tar.gz)
cp $FILE /io/
chown $(stat /io -c %u:%g) /io/$FILE

echo "See result $FILE"
