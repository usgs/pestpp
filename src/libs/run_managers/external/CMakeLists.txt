# This CMake file is part of PEST++

add_library(rm_external RunManagerExternal.cpp)

target_include_directories(rm_external INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(rm_external PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(rm_external rm_abstract)

if(BUILD_SHARED_LIBS)
  install(TARGETS rm_external LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
