# This CMake file is part of PEST++

add_executable(pestpp-sen
  GsaAbstractBase.cpp
  main.cpp
  MorrisMethod.cpp
  PooledVariance.cpp
  sobol.cpp
  Stats.cpp
  # TornadoPlot.cpp
)

target_compile_options(pestpp-sen PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp-sen
  rm_serial
  rm_yamr
  rm_external
)

install(TARGETS pestpp-sen RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
