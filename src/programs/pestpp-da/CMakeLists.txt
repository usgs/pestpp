# This CMake file is part of PEST++

add_executable(pestpp-da pestpp-da.cpp)

target_compile_options(pestpp-da PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-da
  rm_serial
  rm_yamr
  rm_external
)

install(TARGETS pestpp-da RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
