# This CMake file is part of PEST++

add_executable(pestpp-opt pestpp-opt.cpp)

target_compile_options(pestpp-opt PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-opt
  opt
  rm_external
  rm_serial
  rm_yamr
)

install(TARGETS pestpp-opt RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
