# This is the top-level CMake file for PEST++
cmake_minimum_required(VERSION 3.9)

# Get PESTPP_VERSION from common source header
file(READ "${CMAKE_SOURCE_DIR}/src/libs/common/config_os.h" _file)
string(REGEX MATCH "PESTPP_VERSION \"([0-9\.]+)\"" _ ${_file})
set(PESTPP_VERSION ${CMAKE_MATCH_1})
if("${PESTPP_VERSION}" STREQUAL "")
  message(SEND_ERROR
    "Could not find PESTPP_VERSION in src/libs/common/config_os.h")
endif()

project(PESTPP
  VERSION ${PESTPP_VERSION}
  LANGUAGES CXX)

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
  set(Fortran_ENABLED TRUE)
else()
  set(Fortran_ENABLED FALSE)
endif()

find_package(Threads REQUIRED)

# Ensure C++11 is used
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Global warning flags
set(PESTPP_CXX_WARN_FLAGS "${PESTPP_CXX_WARN_FLAGS}"
  CACHE STRING "C++ flags used to compile PEST++ targets")

option(BUILD_SHARED_LIBS "Build PEST++ shared (default is static)" OFF)
message(STATUS "BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}")

# Default install directory is custom layout in source dir ./bin
if(WIN32)
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(_osdir "iwin")
  else()
    set(_osdir "win")
  endif()
elseif(APPLE)
  set(_osdir "mac")
elseif(UNIX)
  set(_osdir "linux")
else()
  message(SEND_ERROR "Could not detect OS")
endif()
set(_default_INSTALL_PREFIX "${PESTPP_SOURCE_DIR}/bin/${_osdir}")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${_default_INSTALL_PREFIX}"
    CACHE PATH "Default custom install prefix" FORCE)
endif()

if("${CMAKE_INSTALL_PREFIX}" STREQUAL "${_default_INSTALL_PREFIX}")
  if(BUILD_SHARED_LIBS AND UNIX)
    message(WARNING
      "Shared libraries may not be found correctly when installed to default "
      "source dir, unless LD_LIBRARY_PATH is also set to that path. "
      "It is advised to use BUILD_SHARED_LIBS=OFF or set "
      "CMAKE_INSTALL_PREFIX to a common system prefix.")
  endif()
  message(STATUS
    "Install location is default, to local source tree: "
    "${CMAKE_INSTALL_PREFIX}")
  set(BINDIR "${CMAKE_INSTALL_PREFIX}")
  set(LIBDIR "${CMAKE_INSTALL_PREFIX}")
  # Remove old files before installing new ones
  install(CODE
    "MESSAGE(\"Cleaning up previous files in ${CMAKE_INSTALL_PREFIX}\")")
  install(CODE "file(REMOVE_RECURSE ${CMAKE_INSTALL_PREFIX})")
else()
  message(STATUS
    "Install prefix set by CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
  set(BINDIR "bin")
  set(LIBDIR "lib")
endif()

# Build main targets in src
add_subdirectory(src)
