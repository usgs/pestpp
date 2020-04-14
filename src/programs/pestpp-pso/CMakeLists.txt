# This CMake file is part of PEST++

add_executable(pestpp-pso
  bestsubs.f90
  estimation.f90
  objsubs.f90
  pareto.f90
  paretosubs.f90
  parsubs.f90
  parunc.f90
  pertsubs.f90
  psomain.f90
  readsubs.f90
  rmsubs.f90
  utilsubs.f90
  writesubs.f90
)

target_link_libraries(pestpp-pso
  rm_wrappers
  ${LAPACK_LIBRARIES}
)

if(NOT "${BUILD_SHARED_LIBS}" AND
    "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  set_property(TARGET pestpp-pso PROPERTY LINKER_LANGUAGE Fortran)
endif()

install(TARGETS pestpp-pso RUNTIME DESTINATION ${BINDIR})
