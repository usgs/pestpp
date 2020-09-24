# This CMake file is part of PEST++

add_subdirectory(gsa)
add_subdirectory(pest++)
add_subdirectory(pestpp-ies)
add_subdirectory(pestpp-sqp)
add_subdirectory(pestpp-opt)
add_subdirectory(sweep)
add_subdirectory(pestpp-mou)
if(Fortran_ENABLED)
  add_subdirectory(pestpp-pso)
endif()
