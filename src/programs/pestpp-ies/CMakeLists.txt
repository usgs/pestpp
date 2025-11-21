# This CMake file is part of PEST++

add_executable(pestpp-ies pestpp-ies.cpp)

target_compile_options(pestpp-ies PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-ies
  rm_serial
  rm_yamr
  rm_external
)

install(TARGETS pestpp-ies RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
