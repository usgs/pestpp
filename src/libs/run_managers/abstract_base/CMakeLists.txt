# This CMake file is part of PEST++

# This library depends on pestpp_com; avoid cyclic dependency
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

target_link_libraries(rm_abstract pestpp_com)

if(BUILD_SHARED_LIBS)
  set_property(TARGET rm_abstract PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
