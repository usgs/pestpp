#!/bin/bash
# Build a portable linux tarball locally, in the same image CI uses.
#
#   scripts/build_pestpp_linux_manylinux.sh
#
# Run it from anywhere in the repo; the tarball lands in the repo root.
#
# On apple silicon the image is emulated - it is x86-64 only - so this is slow and the binaries
# it produces will not run on the host. Fine for checking linkage and the glibc floor, painful
# as a development loop. CI is the place for the real builds.
set -euo pipefail

IMAGE="${PESTPP_MANYLINUX_IMAGE:-quay.io/pypa/manylinux_2_28_x86_64}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

docker_cmd=docker
command -v docker >/dev/null 2>&1 || docker_cmd=podman
command -v "$docker_cmd" >/dev/null 2>&1 || {
    echo "need docker or podman on PATH" >&2
    exit 1
}

# --platform because the image is x86-64 and the host may not be
exec "$docker_cmd" run --rm -t \
    --platform linux/amd64 \
    -v "$REPO":/io \
    "$IMAGE" \
    bash /io/scripts/_build_pestpp_linux_manylinux.sh
