# This CMake file is part of PEST++
add_library(zlib
        adler32.c
        compress.c
        crc32.c
        deflate.c
        gzclose.c
        gzlib.c
        gzread.c
        gzwrite.c
        inflate.c
        infback.c
        inftrees.c
        inffast.c
        trees.c
        uncompr.c
        zutil.c
)

target_include_directories(zlib INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_options(zlib PRIVATE ${PESTPP_CXX_WARN_FLAGS})

if(FORCE_STATIC)
    target_link_libraries(zlib -static)
endif()

target_compile_features(zlib PUBLIC cxx_std_11)
set_target_properties(zlib PROPERTIES CXX_EXTENSIONS OFF)

if(WIN32)
    # Winsock library required
    target_link_libraries(zlib ws2_32)
endif()

if(BUILD_SHARED_LIBS)
    install(TARGETS common LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

