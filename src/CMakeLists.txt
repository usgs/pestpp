# This CMake file is part of PEST++

# Guard against attempting running CMake starting from this directory
# as was the case with the previous GNU Makefile setup
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  message(FATAL_ERROR "Run cmake from the base directory, not 'src'")
endif()

add_subdirectory(libs)
add_subdirectory(programs)
add_subdirectory(utilities)
