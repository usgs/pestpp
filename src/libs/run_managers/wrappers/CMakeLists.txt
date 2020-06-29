# This CMake file is part of PEST++

# This library is only linked by pestpp-pso, so force it static
add_library(rm_wrappers STATIC
  RunManagerCWrapper.cpp
  RunManagerFortranWrapper.cpp
)

if(WIN32)
  target_compile_definitions(rm_wrappers PRIVATE _SCL_SECURE_NO_WARNINGS EXTERNAL_UPPER)
endif()

target_compile_options(rm_wrappers PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(rm_wrappers
  rm_serial
  rm_yamr
)
