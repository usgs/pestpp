# This CMake file is part of PEST++

add_executable(pestpp-glm pest++.cpp)

target_compile_options(pestpp-glm PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-glm
  rm_external
  rm_serial
  rm_yamr
)

install(TARGETS pestpp-glm RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
