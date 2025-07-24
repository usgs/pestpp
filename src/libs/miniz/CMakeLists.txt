# This CMake file is part of PEST++
add_library(miniz
        miniz.c
)

target_include_directories(miniz INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(miniz PRIVATE ${PESTPP_CXX_WARN_FLAGS})

if(FORCE_STATIC)
    target_link_libraries(miniz -static)
endif()

target_compile_features(miniz PUBLIC cxx_std_11)
set_target_properties(miniz PROPERTIES CXX_EXTENSIONS OFF)

if(WIN32)
    # Winsock library required
    target_link_libraries(miniz ws2_32)
endif()

if(BUILD_SHARED_LIBS)
    install(TARGETS common LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

