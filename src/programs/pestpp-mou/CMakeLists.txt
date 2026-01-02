# This CMake file is part of PEST++

add_executable(pestpp-mou pestpp-mou.cpp)

target_compile_options(pestpp-mou PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-mou
  rm_serial
  rm_yamr
  rm_external
)

install(TARGETS pestpp-mou RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
