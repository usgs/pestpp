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

# Ships as pestpp-api.<so|dll|dylib>, matching the pestpp-* naming of the executables it sits
# beside. PREFIX is cleared because posix would otherwise produce libpestpp-api.so, which
# would be the only artifact in bin/ wearing a "lib" prefix. The CMake target keeps its
# pestpp_capi name - that is internal, and renaming it would churn every link line for no
# gain. A side effect worth knowing: the bin/<ostag> copy globs pestpp-*, which this now
# matches on its own.
set_target_properties(pestpp_capi PROPERTIES
  OUTPUT_NAME "pestpp-api"
  PREFIX "")

target_compile_options(pestpp_capi PRIVATE ${PESTPP_CXX_WARN_FLAGS})

# PESTPP_API_BUILD flips PESTPP_API from dllimport to dllexport, and it is PRIVATE on purpose:
# this is the only target that may export these symbols. A downstream consumer including the
# installed header must NOT define it - on MSVC that would re-export our symbols from their
# binary, and for the data symbols (PESTPP_NAME_LEN, PESTPP_MESSAGE_LEN) it is a link error.
target_compile_definitions(pestpp_capi PRIVATE PESTPP_API_BUILD)

# Export the C ABI and nothing else, stated on the target rather than inherited.
#
# The project-wide default is "hidden", but only when BUILD_SHARED_LIBS is OFF - the top-level
# CMakeLists relaxes it to "default" for a shared build. This target is SHARED either way, so
# inheriting that would mean the export surface of the C ABI depends on how the rest of the
# project was configured: with BUILD_SHARED_LIBS=ON it would export every symbol in
# pestpp_com, Eigen and the vendored Coin-OR solvers. That is a real hazard, not a tidiness
# one - the whole point of this library is to be loaded alongside other native libraries
# (libmf6, say) in one process, and a wide export table is how symbol interposition happens.
set_target_properties(pestpp_capi PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  WINDOWS_EXPORT_ALL_SYMBOLS OFF)

# On elf, also refuse to re-export anything that came in from a static archive.
#
# The visibility preset above governs how THIS target's objects are compiled; it says nothing
# about the archives linked into it. Those are compiled with the project-wide preset, so today
# they are hidden too - but that makes the export surface depend on how the rest of the
# project happened to be configured, which is exactly the coupling this target is trying to
# get away from. --exclude-libs states it locally instead: whatever the archives contain, none
# of it is re-exported from here. Safe because this library is a leaf loaded by ctypes with
# all its dependencies linked in statically, so nothing outside it needs those symbols.
#
# GNU ld and lld only; apple's ld64 has no equivalent and does not need one, since mach-o does
# not re-export archive members by default. Probed rather than assumed: this is hardening on
# top of a guarantee the visibility preset already provides, so a linker that does not take the
# option should skip it, not fail to build.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_LINK_OPTIONS "-Wl,--exclude-libs,ALL")
  check_cxx_source_compiles("int main(){return 0;}" PESTPP_LD_HAS_EXCLUDE_LIBS)
  unset(CMAKE_REQUIRED_LINK_OPTIONS)
  if(PESTPP_LD_HAS_EXCLUDE_LIBS)
    target_link_options(pestpp_capi PRIVATE "LINKER:--exclude-libs,ALL")
  endif()

  # ...and a version script, which is the part that actually settles it.
  #
  # --exclude-libs governs static archives; the visibility preset governs how this target's
  # objects are compiled. Neither covers the case that turned up in CI: libstdc++ declares
  # `namespace std _GLIBCXX_VISIBILITY(default)`, so template instantiations this file uses
  # are emitted weak-and-visible out of our OWN object regardless of -fvisibility=hidden.
  # A version script is applied at link time to everything, so it does not care.
  set(_pestpp_api_version_script "${CMAKE_CURRENT_SOURCE_DIR}/pestpp-api.version")
  set(CMAKE_REQUIRED_LINK_OPTIONS "-Wl,--version-script=${_pestpp_api_version_script}")
  check_cxx_source_compiles("int main(){return 0;}" PESTPP_LD_HAS_VERSION_SCRIPT)
  unset(CMAKE_REQUIRED_LINK_OPTIONS)
  if(PESTPP_LD_HAS_VERSION_SCRIPT)
    target_link_options(pestpp_capi PRIVATE
      "LINKER:--version-script=${_pestpp_api_version_script}")
    # relink when the script changes, which cmake would not otherwise notice
    set_property(TARGET pestpp_capi APPEND PROPERTY
      LINK_DEPENDS "${_pestpp_api_version_script}")
  endif()
endif()

target_link_libraries(pestpp_capi
  pestpp_com
  rm_abstract
  rm_serial
  rm_yamr
  rm_external
)

# Ships in bin/ next to the executables, on every platform.
#
# The default would put it in lib/ on posix and bin/ on windows, because a .dll is a RUNTIME
# artifact while a .so/.dylib is a LIBRARY one. That split is right for a library you link
# against and wrong for this one: it is loaded at runtime by ctypes, pest++ distributes a flat
# bin/ directory, and a python caller should not have to know which of two directories their
# platform chose. The import library is the exception - it is link-time only, so it stays in
# lib/ where it does not clutter the distribution.
install(TARGETS pestpp_capi
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}    # windows .dll
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINDIR}    # posix .so / .dylib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})   # windows import .lib

# The header goes with it: a shared library with a C ABI is not much use to anyone writing
# against it from C or C++ without the declarations.
install(FILES pestpp-api.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# the python layer ships too, into python/ next to bin/ and include/. without this the release
# tarball had the library and the header but nothing to call them with, so anybody wanting the
# friendly layer had to go clone the repo and get the commit right.
#
# the files stay flat and importable rather than being made into a package - there is no pip
# install yet on purpose (see python/README.md), so the way you use these is to put the
# directory on sys.path. find_library() looks in ../bin relative to itself, which is exactly
# where the library lands in this layout, so it works out of the tarball with nothing set.
install(FILES
  ${CMAKE_SOURCE_DIR}/python/pestpp.py
  ${CMAKE_SOURCE_DIR}/python/pestpp_lib.py
  ${CMAKE_SOURCE_DIR}/python/pestpp_progress.py
  ${CMAKE_SOURCE_DIR}/python/README.md
  DESTINATION python)
