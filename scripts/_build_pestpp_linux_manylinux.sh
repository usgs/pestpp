#!/bin/bash
# Runs INSIDE the manylinux container. Use the wrapper rather than calling this directly:
#   scripts/build_pestpp_linux_manylinux.sh
#
# Replaces the holy-build-box pair for the same job: a linux build that runs somewhere other
# than the machine it was compiled on. manylinux gets used rather than hbb because the wheel
# has to be built this way regardless - pypi rejects a plain linux_x86_64 tag - so one image
# covers the tarball and the wheel, and auditwheel gives a mechanical check that hbb's libcheck
# only advises on.
#
# What makes the result portable is two separate things:
#   - the IMAGE sets the glibc floor. manylinux_2_28 is glibc 2.28 (rhel 8, ubuntu 18.10,
#     debian 10 and newer). glibc stays dynamic on purpose - see the note in CMakeLists.txt
#     about NSS and getaddrinfo, which PANTHER needs to resolve the master by hostname.
#   - CMakeLists.txt links libstdc++ and libgcc STATICALLY on linux, so the binaries do not ask
#     the host for a GLIBCXX_3.4.xx symbol version it may not have. That is the failure mode
#     that bites more often than glibc itself.
set -euo pipefail

# the image ships several pythons, none of them on PATH, and no system cmake
export PATH="/opt/python/cp312-cp312/bin:$PATH"
pip install --quiet cmake ninja

cd /io
rm -rf build_manylinux
cmake -S . -B build_manylinux -G Ninja -DCMAKE_BUILD_TYPE=Release -DINSTALL_LOCAL=OFF
cmake --build build_manylinux -j"$(nproc)"

cd build_manylinux
cpack -G TGZ

echo
echo "=== what the binaries actually require ==="
for exe in bin/pestpp-ies bin/pestpp-glm; do
    [ -f "$exe" ] || continue
    echo "--- $exe"
    # the highest glibc symbol version needed IS the floor - anything older can run it
    readelf -V "$exe" 2>/dev/null | grep -oE 'GLIBC_2\.[0-9]+' | sort -u -t. -k2 -n | tail -1
    # libstdc++ must NOT appear here; if it does, the static link did not take
    ldd "$exe" 2>/dev/null | grep -E 'libstdc\+\+|libgcc_s' \
        && echo "  WARNING: c++ runtime is still dynamic" \
        || echo "  c++ runtime: static (good)"
done

FILE=$(ls pestpp-*.tar.gz)
cp "$FILE" /io/
# the container runs as root, so the copy would otherwise land in the host tree owned by root
chown "$(stat -c %u:%g /io)" "/io/$FILE"
echo
echo "wrote $FILE"
