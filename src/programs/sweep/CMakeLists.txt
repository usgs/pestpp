# This CMake file is part of PEST++

add_executable(pestpp-swp sweep.cpp)

target_compile_options(pestpp-swp PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-swp
  rm_serial
  rm_yamr
  rm_external
)

install(TARGETS pestpp-swp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
