# This CMake file is part of PEST++

add_library(mio mio.F90)

target_link_libraries(mio common)

if(BUILD_SHARED_LIBS)
  install(TARGETS mio LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
