# This CMake file is part of PEST++

add_library(rm_serial RunManagerSerial.cpp)

target_include_directories(rm_serial INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(rm_serial PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(rm_serial rm_abstract)

if(BUILD_SHARED_LIBS)
  install(TARGETS rm_serial LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
