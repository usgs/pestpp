# This CMake file is part of PEST++

add_subdirectory(abstract_base)
add_subdirectory(external)
add_subdirectory(serial)
add_subdirectory(yamr)

if(Fortran_ENABLED)
  add_subdirectory(wrappers)
endif()
