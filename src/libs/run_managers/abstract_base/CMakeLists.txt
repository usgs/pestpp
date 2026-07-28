# This CMake file is part of PEST++

add_library(rm_abstract STATIC
  debug.cpp
  linpackc.cpp
  model_interface.cpp
  RunManagerAbstract.cpp
  RunStorage.cpp
  Serialization.cpp
)

target_include_directories(rm_abstract INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(rm_abstract PRIVATE ${PESTPP_CXX_WARN_FLAGS})

# rm_abstract sits BELOW pestpp_com in the stack: it depends only on the shared low-level
# pieces. Linking it back to pestpp_com was the other half of the dependency cycle that
# made any shared-library build impossible to even configure.
target_link_libraries(rm_abstract common)

# PIC is set globally now (CMAKE_POSITION_INDEPENDENT_CODE in the top-level CMakeLists)
