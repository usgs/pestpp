# This CMake file is part of PEST++

# The C ABI. Built shared by default regardless of BUILD_SHARED_LIBS: a static C ABI would
# be pointless, since the whole purpose is to be dlopen'd by ctypes.
add_library(pestpp_capi SHARED
  pestpp_capi.cpp
)

target_include_directories(pestpp_capi
  INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_compile_options(pestpp_capi PRIVATE ${PESTPP_CXX_WARN_FLAGS})

# The entry points are marked with explicit visibility attributes in the header, so this
# target keeps the project-wide "hidden" default and exports only the C ABI.
target_link_libraries(pestpp_capi
  pestpp_com
  rm_abstract
  rm_serial
  rm_yamr
)

install(TARGETS pestpp_capi LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
