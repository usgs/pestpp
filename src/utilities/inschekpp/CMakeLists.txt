# This CMake file is part of PEST++

add_executable(inschekpp
  inschek.F
  nblnk.F
  pgetcl.F
  space.F
)

install(TARGETS inschekpp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
