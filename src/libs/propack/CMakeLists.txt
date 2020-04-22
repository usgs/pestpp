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

if(BUILD_SHARED_LIBS)
  install(TARGETS propack LIBRARY DESTINATION ${LIBDIR})
else()
  # pthread required for static linking
  target_link_libraries(propack Threads::Threads)
endif()
