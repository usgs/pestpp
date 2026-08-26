# This is the top-level CMake file for PEST++
cmake_minimum_required(VERSION 3.15)
set(CMAKE_VERBOSE_MAKEFILE,ON)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_definitions(_DEBUG)
endif()

# Get the version from the common source header, which is the single definition of it.
#
# The version may carry a release-candidate suffix ("5.2.29rc1"), and that has to be split:
# project(VERSION) accepts only numeric major.minor.patch, and CPack's version fields are
# numeric too. The FULL string is what names the package, so an rc archive is distinguishable
# from the release at a glance.
#
# Held in PESTPP_VERSION_FULL rather than PESTPP_VERSION on purpose: project() SETS
# ${PROJECT_NAME}_VERSION, so a variable of that name is silently overwritten a few lines
# below with the numeric form - and the suffix would vanish with it.
file(READ "${CMAKE_SOURCE_DIR}/src/libs/common/config_os.h" _file)
string(REGEX MATCH "PESTPP_VERSION \"([0-9]+\\.[0-9]+\\.[0-9]+[A-Za-z0-9._-]*)\"" _ ${_file})
set(PESTPP_VERSION_FULL ${CMAKE_MATCH_1})
if("${PESTPP_VERSION_FULL}" STREQUAL "")
  message(SEND_ERROR
    "Could not find a PESTPP_VERSION of the form \"X.Y.Z\" or \"X.Y.Zrc<n>\" in "
    "src/libs/common/config_os.h")
endif()
string(REGEX REPLACE "^([0-9]+\\.[0-9]+\\.[0-9]+).*$" "\\1"
  PESTPP_VERSION_NUMERIC "${PESTPP_VERSION_FULL}")

message(STATUS
  "Configuring CMake ${CMAKE_VERSION} to build PESTPP ${PESTPP_VERSION_FULL}")

# This logic needs to be considered before project()
#
# /bigobj used to live here and did nothing. MSVC is set by compiler identification, which does
# not happen until project(), so this if() is false on a fresh configure - unlike WIN32 below,
# which comes from the host and IS set this early. It went unnoticed because SQP.cpp sat just
# under the object-file section limit without it; adding /Zi pushed it over and the build died
# with C1128. It is now in the MSVC block AFTER project(). Do not move it back up here.
set(_change_MSVC_flags FALSE)
if(WIN32)
  add_compile_definitions(NOMINMAX)
  if(CMAKE_VERSION VERSION_LESS 3.15.0)
    set(_change_MSVC_flags TRUE)
  else()
    # MSVC runtime. /MT (static CRT) is the default because the shipped executables are
    # meant to run without a redistributable. A CPython extension module, however, must be
    # built against the DLL runtime (/MD) to share a heap and CRT state with the
    # interpreter, so make this selectable rather than hard-coded.
    cmake_policy(SET CMP0091 NEW)
    option(PESTPP_MSVC_STATIC_RUNTIME "Link the static MSVC runtime (/MT). OFF gives /MD, required for Python extension modules." ON)
    if(PESTPP_MSVC_STATIC_RUNTIME)
      set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    else()
      set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL$<$<CONFIG:Debug>:Debug>")
    endif()
    message(STATUS "PESTPP_MSVC_STATIC_RUNTIME=${PESTPP_MSVC_STATIC_RUNTIME} -> ${CMAKE_MSVC_RUNTIME_LIBRARY}")
    # /MT gives pestpp-api.<dll> a PRIVATE C runtime: its own heap and its own file descriptor
    # table, not shared with a host process that dlopen's it. That is fine, and deliberately
    # left as the default, because the C ABI never lets ownership cross that boundary - no
    # entry point returns memory the caller frees (the one double** is a borrowed view of the
    # library's own buffer, returned via pestpp_release_view), and no FILE* or descriptor
    # appears in the ABI. Keep it that way; see the ownership note in pestpp-api.h.
    #
    # This was NOT always true. pestpp_redirect_output() used to _dup2 file descriptor 1, which
    # under /MT closed a console handle the host was still using - python's next write failed
    # with WinError 6 and the freed handle value was recycled onto an open output file. That is
    # fixed in the library rather than by changing runtimes, because switching would have to be
    # whole-project anyway: pestpp_com and the run managers are STATIC archives
    # (BUILD_SHARED_LIBS defaults OFF), so /MD on the shared target alone would link /MT
    # archives into an /MD module - two runtimes inside one dll, which is the same bug moved.
    if(NOT PESTPP_MSVC_STATIC_RUNTIME)
      message(STATUS "  (shared runtime: required if you link the static libs into a CPython extension)")
    endif()
  endif()
endif()

project(PESTPP
  VERSION ${PESTPP_VERSION_NUMERIC}
  HOMEPAGE_URL "https://github.com/usgs/pestpp/"
  LANGUAGES CXX)

# Prevent running cmake in source dir
if("${PESTPP_SOURCE_DIR}" STREQUAL "${PESTPP_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed")
endif()

# Debug symbols for the SHIPPED windows build.
#
# Without a .pdb, a crash report from a released binary is undiagnosable: a windows crash dump
# records addresses and the module identity, nothing else, so every pestpp frame in the stack
# reads pestpp-ies.exe+0x3f2a1 and stops there. The pdb is the table that turns that back into a
# file and a line. It is never read while the program runs - only by whoever opens the dump.
#
# The linker stamps a matching guid into the exe and the pdb, so a pdb only ever works for the
# exact build that produced it. That is why these have to be kept at release time; they cannot
# be regenerated for a binary someone downloaded months ago.
#
# /Zi is compiler-side debug info and /DEBUG makes the linker write the pdb. Neither changes
# codegen, so this is the same executable we ship today. /DEBUG DOES turn off the linker's
# default /OPT:REF and /OPT:ICF though, so put both back explicitly - otherwise the exe grows
# and stops folding identical code, which would be a real change and not just an extra file.
if(MSVC)
  # every config, not just Release: SQP.cpp exceeds the default section limit once any debug
  # info is emitted, and the diagnostic RelWithDebInfo build emits it for every file
  add_compile_options(/bigobj)
  add_compile_options($<$<CONFIG:Release>:/Zi>)
  add_link_options($<$<CONFIG:Release>:/DEBUG>
                   $<$<CONFIG:Release>:/OPT:REF>
                   $<$<CONFIG:Release>:/OPT:ICF>)
endif()

# Let the binary say which build it is. startup_report() prints this into every record file, so
# a user sending one back identifies their own build without being asked. $<CONFIG> rather than
# CMAKE_BUILD_TYPE so it is still right under a multi-config generator, where the build type is
# not known until build time.
add_compile_definitions(PESTPP_BUILD_CONFIG="$<CONFIG>")

# Fortran is not required (anymore), but can be enabled using either:
#  1. ENABLE_Fortran=ON and use whatever default compiler is available
#  2. CMAKE_Fortran_COMPILER=ifort (or other compiler)
if(DEFINED ENABLE_Fortran)
  set(_default_ENABLE_Fortran ${ENABLE_Fortran})
else()
  if(DEFINED CMAKE_Fortran_COMPILER)
    set(_default_ENABLE_Fortran ON)
  else()
    set(_default_ENABLE_Fortran OFF)
  endif()
endif()
option(ENABLE_Fortran "Enable Fortran in PEST++" ${_default_ENABLE_Fortran})
if(ENABLE_Fortran)
  enable_language(Fortran)
  find_package(LAPACK REQUIRED)
  message(STATUS "BLA_VENDOR=${BLA_VENDOR}")
  message(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")
  message(STATUS "LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
  set(local_install_pso pestpp-pso)
  set(Fortran_ENABLED TRUE)
else()
  set(Fortran_ENABLED FALSE)
endif()

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
message(STATUS "CMAKE_THREAD_LIBS_INIT=${CMAKE_THREAD_LIBS_INIT}")

# C++17. Not C++20: the bundled Eigen 3.3.9 uses std::result_of, removed in C++20, and the
# remaining bare throw() specifications become errors there. Both are fixable, but they are
# a separate job (Eigen 3.4+ plus a noexcept sweep).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_TRY_COMPILE_TARGET_TYPE)

add_compile_definitions(_HAS_STD_BYTE=0)

# Position-independent code for every target, not just shared builds. Objects that are not
# PIC cannot be linked into a shared library or a Python extension module on Linux, so
# gating this on BUILD_SHARED_LIBS meant the static libs could never be reused that way.
# The cost on a normal static build is negligible.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(BUILD_SHARED_LIBS "Build PEST++ shared (default is static, OFF)" OFF)
message(STATUS "BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}")

# Symbol visibility. Hidden is right for the executables - smaller, faster to load - but a
# shared library that hides everything exports nothing and is useless to link against.
# Exporting all symbols is the pragmatic choice: annotating the API with export macros is a
# separate job, and the point of the shared build today is that it works at all.
# See https://gcc.gnu.org/wiki/Visibility
if(BUILD_SHARED_LIBS)
  set(CMAKE_CXX_VISIBILITY_PRESET default)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN OFF)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
endif()

# Global warning flags
set(PESTPP_CXX_WARN_FLAGS "${PESTPP_CXX_WARN_FLAGS}"
  CACHE STRING "C++ flags used to compile PEST++ targets")

# Global definitions
if(WIN32)
  add_compile_definitions(_WIN32_WINNT=0x0600)
endif()

option(FORCE_STATIC "Force '-static' flag to link applications (default OFF)" OFF)
message(STATUS "FORCE_STATIC=${FORCE_STATIC}")

# Default install directories
include(GNUInstallDirs)

# Build main targets in src
add_subdirectory(src)

# Define an custom OS tag for CPACK_SYSTEM_NAME and (for now) local bin sub-directory
message(STATUS
  "CMAKE compiler ID ${CMAKE_CXX_COMPILER_ID}")
if(WIN32)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
    set(_ostag "iwin")
  else()
    set(_ostag "win")
  endif()
elseif(APPLE)
  set(_ostag "mac")
elseif(UNIX)
  set(_ostag "linux")
else()
  set(_ostag "unknown")
endif()

# Detect if using a 32-bit compiler, show warning and modify _ostag
if(CMAKE_SIZEOF_VOID_P STREQUAL 4)
  message(WARNING
    "A 32-bit compiler has been detected. It is recommended that "
    "PEST++ should be compiled using a 64-bit compiler.")
  set(_ostag "${_ostag}32")
endif()

# Use a custom command to locally install binaries after compiling
# https://www.linux.com/training-tutorials/cmake-recipe-2-install-local-folder-build-dir-testing/
option(INSTALL_LOCAL "Install local executables in ./bin" ON)
if(INSTALL_LOCAL)
  set(_OPTION_INSTALL_LOCAL_ALL ALL)
endif()
# bin/${_ostag}/ is NOT optional and must not be "cleaned up" away. pyemu puts
# os.path.abspath("../bin/<plat>") on PATH at import time (pyemu/utils/os_utils.py), so every
# benchmark script that starts a tool by bare name - `start_workers(t_d, "pestpp-mou", ...)` -
# resolves it through this directory. That includes the cloned pestpp-*_benchmarks repos, which
# are not ours to edit. Deleting this layout breaks them; it was tried, and it did.
#
# What IS fixed here is how the copy is made. It used to be a shell command run straight after
# install():
#     cd bin && rm -rf ${_ostag} && mkdir ${_ostag} && cp pestpp-* ${_ostag}/
# which re-read ten multi-megabyte files tens of milliseconds after install() wrote them. On one
# macOS runner every file above ~5 MB died with `cp: fcopyfile failed: Input/output error` while
# the two below ~4 MB survived - a clean size ordering, not a random fault - leaving bin/mac/
# full of truncated binaries that then failed to exec. The `rm -rf` also blanked the directory
# for the length of the copy, so anything running concurrently saw it empty.
#
# Copying each target from the BUILD TREE avoids both: the build-tree binary was written at link
# time rather than milliseconds ago, copy_if_different is one defined operation per file with no
# glob, and nothing is ever removed, so the directory is never momentarily empty.
set(_local_install_bin "${PESTPP_SOURCE_DIR}/bin/${_ostag}")
set(_local_install_copies
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${_local_install_bin}")
foreach(_tgt
    pestpp-glm pestpp-ies pestpp-opt pestpp-sen pestpp-swp
    pestpp-mou pestpp-sqp pestpp-da pestpp-selftest pestpp_capi)
  list(APPEND _local_install_copies
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
      "$<TARGET_FILE:${_tgt}>" "${_local_install_bin}/")
  # ship the pdb beside the exe so the archive carries its own symbols. OPTIONAL because
  # TARGET_PDB_FILE is only written when the linker was asked for one, and MSVC only - the
  # generator expression is not valid for targets that cannot produce a pdb at all
  if(MSVC)
    install(FILES "$<TARGET_PDB_FILE:${_tgt}>"
            DESTINATION ${CMAKE_INSTALL_BINDIR} OPTIONAL)
  endif()
endforeach()

add_custom_target(local_install
  ${_OPTION_INSTALL_LOCAL_ALL}
  "${CMAKE_COMMAND}"
    -D CMAKE_INSTALL_PREFIX=${PESTPP_SOURCE_DIR}
    -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
  ${_local_install_copies}
  DEPENDS
    pestpp-glm
    pestpp-ies
    pestpp-opt
    pestpp-sen
    pestpp-swp
    pestpp-mou
    pestpp-sqp
    pestpp-da
    pestpp_capi
    # not shipped, but installed so bin/ never serves a stale selftest - and listed HERE
    # because its install() rule would otherwise let ninja run cmake_install.cmake before the
    # binary exists. Anything given an install() rule belongs in this list.
    pestpp-selftest
  COMMENT "Installing local executables and the shared library to ${PESTPP_SOURCE_DIR}/bin"
)

# Packaging
string(TOLOWER ${PROJECT_NAME} CPACK_PACKAGE_NAME)
set(CPACK_PACKAGE_VERSION_MAJOR ${PESTPP_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PESTPP_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PESTPP_VERSION_PATCH})
# the FULL version, suffix and all, so a release candidate names itself in the archive
set(CPACK_PACKAGE_VERSION ${PESTPP_VERSION_FULL})
set(CPACK_SYSTEM_NAME ${_ostag})
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY OFF)
include(CPack)
