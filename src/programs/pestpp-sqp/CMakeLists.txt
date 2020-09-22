# This CMake file is part of PEST++

add_executable(pestpp-sqp pestpp-sqp.cpp)

target_compile_options(pestpp-sqp PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-sqp
  rm_serial
  rm_yamr
)

install(TARGETS pestpp-sqp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
