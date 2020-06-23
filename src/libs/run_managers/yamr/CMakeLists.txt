# This CMake file is part of PEST++

add_library(rm_yamr
  PantherAgent.cpp
  RunManagerPanther.cpp
)

target_include_directories(rm_yamr INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(rm_yamr PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(rm_yamr rm_abstract)

if(BUILD_SHARED_LIBS)
  install(TARGETS rm_yamr LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
