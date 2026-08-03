# This CMake file is part of PEST++

add_executable(pestpp-selftest selftest.cpp)

target_compile_options(pestpp-selftest PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-selftest
  rm_serial
  rm_yamr
  rm_agent
  rm_external
)

# Installed, so that bin/ holds a CURRENT selftest.
#
# It was previously left out on the grounds that it is a CI-only binary run from the build
# tree - which CI does do. What that missed is that a stale bin/pestpp-selftest, left behind by
# an older build, never got overwritten, while local_install's `cp pestpp-*` glob copied it
# into bin/<ostag>/ on every build. The result was a binary with a fresh timestamp and
# five-week-old contents, which reported PASS (31/31) where the build tree reported 1278 - a
# green result that meant nothing. Installing it keeps bin/ in step with the source.
#
# The race the old comment warned about is real and is handled at the other end: pestpp-selftest
# is in local_install's DEPENDS, so it is built before cmake_install.cmake runs. Adding an
# install() rule WITHOUT that entry would let ninja install ahead of the build.
install(TARGETS pestpp-selftest RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
