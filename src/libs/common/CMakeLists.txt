# This CMake file is part of PEST++

add_library(common
  fortran_wrappers.cpp
  network_package.cpp
  network_wrapper.cpp
  pest_error.cpp
  system_variables.cpp
  Transformable.cpp
  utilities.cpp
)

target_include_directories(common INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(common PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(common
  Eigen
  Threads::Threads
)

target_compile_features(common PUBLIC cxx_std_11)
set_target_properties(common PROPERTIES CXX_EXTENSIONS OFF)

if(WIN32)
  # Winsock library required
  target_link_libraries(common ws2_32)
endif()

if(BUILD_SHARED_LIBS)
  install(TARGETS common LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
