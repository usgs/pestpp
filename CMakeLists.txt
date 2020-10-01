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

message(STATUS
  "Configuring CMake ${CMAKE_VERSION} to build PESTPP ${PESTPP_VERSION}")

# This logic needs to be considered before project()
set(_change_MSVC_flags FALSE)
if(WIN32)
  if(CMAKE_VERSION VERSION_LESS 3.15.0)
    set(_change_MSVC_flags TRUE)
  else()
    # Set MSVC runtime to MultiThreaded (/MT)
    cmake_policy(SET CMP0091 NEW)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  endif()
endif()

project(PESTPP
  VERSION ${PESTPP_VERSION}
  LANGUAGES CXX)

# Prevent running cmake in source dir
if("${PESTPP_SOURCE_DIR}" STREQUAL "${PESTPP_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed")
endif()

if(_change_MSVC_flags)
  # Modify compile flags to change MSVC runtime from /MD to /MT
  set(_re_match "([\\/\\-]M)D")
  set(_re_replace "\\1T")
  string(REGEX REPLACE ${_re_match} ${_re_replace}
    CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  string(REGEX REPLACE ${_re_match} ${_re_replace}
    CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
  string(REGEX REPLACE ${_re_match} ${_re_replace}
    CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
  string(REGEX REPLACE ${_re_match} ${_re_replace}
    CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL}")
  string(REGEX REPLACE ${_re_match} ${_re_replace}
    CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

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

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
message(STATUS "CMAKE_THREAD_LIBS_INIT=${CMAKE_THREAD_LIBS_INIT}")

# Ensure C++11 is used
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Use global "-fvisibility=hidden" see https://gcc.gnu.org/wiki/Visibility
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Global warning flags
set(PESTPP_CXX_WARN_FLAGS "${PESTPP_CXX_WARN_FLAGS}"
  CACHE STRING "C++ flags used to compile PEST++ targets")

# Global definitions
if(WIN32)
  add_compile_definitions(_WIN32_WINNT=0x0600)
endif()

option(BUILD_SHARED_LIBS "Build PEST++ shared (default is static, OFF)" OFF)
message(STATUS "BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}")

option(FORCE_STATIC "Force '-static' flag to link applications (default OFF)" OFF)
message(STATUS "FORCE_STATIC=${FORCE_STATIC}")

# Default install directories
include(GNUInstallDirs)

# Build main targets in src
add_subdirectory(src)

# Define an custom OS tag for CPACK_SYSTEM_NAME and (for now) local bin sub-directory
if(WIN32)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
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
# Use an extra command to copy binaries into a sub-directory within bin -- please remove soon
if(WIN32)
  set(_local_install_COMMAND_to_be_removed
    cd "${PESTPP_SOURCE_DIR}\\bin" &&
    md ${_ostag} &
    copy /b pestpp-* ${_ostag}\\
  )
else()
  set(_local_install_COMMAND_to_be_removed
    cd "${PESTPP_SOURCE_DIR}/bin" &&
    rm -rf ${_ostag} &&
    mkdir ${_ostag} &&
    cp pestpp-* ${_ostag}/
  )
endif()
add_custom_target(local_install
  ${_OPTION_INSTALL_LOCAL_ALL}
  "${CMAKE_COMMAND}"
    -D CMAKE_INSTALL_PREFIX=${PESTPP_SOURCE_DIR}
    -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
  COMMAND ${_local_install_COMMAND_to_be_removed}
  DEPENDS
    pestpp-glm
    pestpp-ies
    pestpp-opt
    pestpp-sen
    pestpp-swp
    pestpp-mou
    pestpp-sqp
  COMMENT "Installing local executables to ${PESTPP_SOURCE_DIR}/bin"
)

# Packaging
string(TOLOWER ${PROJECT_NAME} CPACK_PACKAGE_NAME)
set(CPACK_PACKAGE_VERSION_MAJOR ${PESTPP_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PESTPP_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PESTPP_VERSION_PATCH})
set(CPACK_SYSTEM_NAME ${_ostag})
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY OFF)
include(CPack)
