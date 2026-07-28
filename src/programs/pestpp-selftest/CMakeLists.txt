# This CMake file is part of PEST++

add_executable(pestpp-selftest selftest.cpp)

target_compile_options(pestpp-selftest PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-selftest
  rm_serial
  rm_yamr
  rm_agent
  rm_external
)

# intentionally NOT installed: this is a CI-only test binary, run from the build tree.
# Installing it would also add it to local_install's cmake_install.cmake without being in
# that target's DEPENDS list, racing ahead of its own build under ninja.
