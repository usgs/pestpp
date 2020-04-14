# This CMake file is part of PEST++

# Create a header-only target for Eigen that other targets can use
add_library(Eigen INTERFACE)

target_include_directories(Eigen INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/Eigen")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # silence these warnings when building CXX objects
  target_compile_options(Eigen
    INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wno-deprecated-declarations>
  )
endif()

if(WIN32)
  add_compile_definitions(_WIN32_WINNT=0x0600)
endif()

add_subdirectory(common)
add_subdirectory(opt)
add_subdirectory(pestpp_common)
add_subdirectory(run_managers)

if(Fortran_ENABLED)
  # These libraries are not linked anywhere
  # add_subdirectory(mio)
  # add_subdirectory(propack)
endif()
