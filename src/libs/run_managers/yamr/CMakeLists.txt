# This CMake file is part of PEST++

# The master (RunManagerPanther) sits BELOW pestpp_com - the tools drive it. The agent
# (PantherAgent) sits ABOVE it - it processes a control file itself, so it needs Pest.
# Keeping them in one target forced a pestpp_com <-> rm_yamr cycle and made a shared build
# impossible to configure, so they are separate targets at their true levels.
add_library(rm_yamr
  RunManagerPanther.cpp
)

add_library(rm_agent
  PantherAgent.cpp
)

target_include_directories(rm_yamr INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(rm_yamr PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(rm_yamr rm_abstract)

target_include_directories(rm_agent INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(rm_agent PRIVATE ${PESTPP_CXX_WARN_FLAGS})
target_link_libraries(rm_agent rm_abstract pestpp_com)

if(BUILD_SHARED_LIBS)
  install(TARGETS rm_yamr LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
