# This CMake file is part of PEST++

add_executable(pestpp-selftest selftest.cpp)

target_compile_options(pestpp-selftest PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-selftest
  rm_serial
  rm_yamr
  rm_external
)

install(TARGETS pestpp-selftest RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
