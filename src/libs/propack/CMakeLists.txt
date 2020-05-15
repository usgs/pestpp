# This CMake file is part of PEST++

add_library(propack
  dgemm_ovwr.F
  dgetu0.F
  dlanbpro.F
  dlansvd.F
  dmgs.pentium.F
  dreorth.F
  dritzvec.F
  dsafescal.F
  propack_misc.for
)

target_link_libraries(propack ${LAPACK_LIBRARIES})

if(WIN32)
  target_compile_definitions(opt PRIVATE __windows__)
endif()

if(BUILD_SHARED_LIBS)
  install(TARGETS propack LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
  # pthread required for static linking
  target_link_libraries(propack Threads::Threads)
endif()
