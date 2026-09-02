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

# Keep the users manual's title page in step with config_os.h.
#
# The docs build already injects this when it renders (docs/preprocess_manual.py), so the
# PUBLISHED manual is never wrong. This is for the copy in the repo - the one github renders
# when someone browses documentation/ - which otherwise sits at whatever version it was last
# hand-edited to. It had drifted to 5.2.28 while the code moved on, and nothing said so.
#
# The version is a line of its own in the manual, deliberately, so this touches those few
# characters and nothing else. It used to be the first HEADING, which meant the table of
# contents mirrored it and two places in two syntaxes had to agree.
#
# Only writes when the string actually differs, so a repeat configure leaves the file alone -
# same mtime, same size, nothing in git status. QUOTE the file contents everywhere below: the
# manual is full of style="width:6.26806in;height:1.68194in" and an unquoted dereference would
# treat every semicolon as a cmake list separator and shred the img tags.
# CMAKE_SOURCE_DIR, not PESTPP_SOURCE_DIR: the latter is set by project(), which has not run
# yet up here, so it would expand to nothing and the EXISTS check would quietly never fire
set(_manual "${CMAKE_SOURCE_DIR}/documentation/pestpp_users_manual.md")
if(EXISTS "${_manual}")
  file(READ "${_manual}" _manual_text)
  string(REGEX REPLACE "\nVersion [^\n]*\n" "\nVersion ${PESTPP_VERSION_FULL}\n"
         _manual_new "${_manual_text}")
  if(NOT _manual_new STREQUAL _manual_text)
    file(WRITE "${_manual}" "${_manual_new}")
    message(STATUS "  users manual version -> ${PESTPP_VERSION_FULL}")
  endif()
  unset(_manual_text)
  unset(_manual_new)
endif()
unset(_manual)

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
# Linux: fold the C++ runtime into the binaries, the way holy-build-box does.
#
# NOT glibc. glibc stays dynamic, deliberately. Statically linking it breaks NSS, which glibc
# implements by dlopen'ing libnss_*.so at runtime - and getaddrinfo goes through NSS, which is
# how a PANTHER agent resolves the master by hostname (PantherAgent.cpp). A static glibc build
# would either fail those lookups or end up needing the very glibc it was trying to avoid.
#
# libstdc++ and libgcc are the ones worth folding in. A binary built with a newer gcc asks the
# host for a GLIBCXX_3.4.xx symbol version it may not have, and that is the usual way a linux
# build dies on someone else's machine - more often than glibc itself. Folding them in removes
# the question. The remaining floor is glibc alone, set by whatever image the build ran in.
#
# Safe for the shared library too: the C ABI passes no C++ types, so a host process with its own
# libstdc++ never meets ours. That is the same argument that makes /MT safe on windows - see the
# ownership note in pestpp-api.h.
if(UNIX AND NOT APPLE)
  add_link_options(-static-libstdc++ -static-libgcc)
endif()

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

# A small, deliberately narrow set of warnings that are on by default.
#
# Not -Wall: a build with -Wall -Wextra reports about 1100 warnings in our own code, nearly all
# of them hygiene (signed/unsigned comparisons, unused variables, ignored qualifiers). At that
# volume nobody reads any of them, which is how the three below went unnoticed for years:
#
#   OperSys::gets_s      called itself forever on mac and linux
#   get_euclidean_fitness returned an uninitialised double
#   two self-assignments in Ensemble::from_binary that did nothing
#
# Each of these is a defect with no plausible innocent reading, and each is caught by exactly
# one of these three flags. The categories are now clean, so anything they report is new.
#
# Errors, not warnings. A warning in a build that is otherwise green does not get read, so as
# warnings these would record the rule without enforcing it - and the whole reason the three
# defects above survived so long is that nobody was reading build output. All three categories
# are at zero, so this cannot fail on anything that exists today: it can only fail on something
# newly written, which is the point.
#
# Checked in the -Werror= form because that is what has to be accepted. Neither Eigen nor the
# vendored Coin-OR trips these - Eigen never fires any of the three, and Coin-OR's three
# self-assignments are turned off at that target in src/libs/opt/CMakeLists.txt.
#
# To back out, configure with -DPESTPP_CXX_WARN_FLAGS= to clear it, or drop the offending entry
# from the list below.
#
# One asymmetry to know about: on gcc, -Wuninitialized also turns on -Wmaybe-uninitialized,
# which is a guess rather than a proof - it fires where gcc cannot see that a loop always runs.
# It earned its keep once: it found w_sendall and w_recvall returning an unassigned int on a
# zero-length transfer, on the first CI build after this went in.
#
# But it cannot stay an error. Building the whole tree with gcc 16 turns up five more, and not
# one of them is ours to fix: two are inside libstdc++/libc++ headers (unique_ptr.h, string) and
# three are in vendored Coin-OR. A warning that fires inside the standard library cannot gate a
# release build - a compiler upgrade would break the build over code nobody here wrote.
#
# So it is turned back off below, AFTER the -Werror flags so it wins. What remains is the
# definite check, which proves rather than guesses, and which is the one that caught the real
# defects. Our own code is clean under gcc 16 for both.
#
# Detected rather than assumed. -Wself-assign is clang-only, -Winfinite-recursion needs gcc 12
# or newer, and gcc treats an unknown -W option as a hard error rather than ignoring it - so
# hardcoding these would break the build on the very compilers we ship from. MSVC accepts none
# of them and simply gets nothing added.
include(CheckCXXCompilerFlag)
foreach(_pestpp_warn_flag -Werror=uninitialized -Werror=infinite-recursion -Werror=self-assign)
  string(MAKE_C_IDENTIFIER "PESTPP_HAS${_pestpp_warn_flag}" _pestpp_warn_ok)
  check_cxx_compiler_flag("${_pestpp_warn_flag}" ${_pestpp_warn_ok})
  if(${_pestpp_warn_ok})
    list(APPEND PESTPP_CXX_WARN_FLAGS "${_pestpp_warn_flag}")
  endif()
endforeach()

# after the loop on purpose - later on the command line, so it overrides the -Werror above for
# the guessing half of -Wuninitialized. see the note above for why this cannot be an error.
# clang has no such warning and skips this.
check_cxx_compiler_flag("-Wno-maybe-uninitialized" PESTPP_HAS_Wno_maybe_uninitialized)
if(PESTPP_HAS_Wno_maybe_uninitialized)
  list(APPEND PESTPP_CXX_WARN_FLAGS "-Wno-maybe-uninitialized")
endif()

message(STATUS "PESTPP_CXX_WARN_FLAGS=${PESTPP_CXX_WARN_FLAGS}")

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
# Ship .dSYM bundles with the binaries (macOS). OFF because they are large - see the note in
# the target loop below. The diagnostic build turns it on; ordinary builds should not.
option(PESTPP_MAC_DSYM "Run dsymutil and install .dSYM bundles alongside the binaries" OFF)

set(_local_install_bin "${PESTPP_SOURCE_DIR}/bin/${_ostag}")
set(_local_install_copies
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${_local_install_bin}")
foreach(_tgt
    pestpp-glm pestpp-ies pestpp-opt pestpp-sen pestpp-swp
    pestpp-mou pestpp-sqp pestpp-da pestpp-selftest pestpp_capi)
  list(APPEND _local_install_copies
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
      "$<TARGET_FILE:${_tgt}>" "${_local_install_bin}/")
  # NOT in Release: the shipped archive stays as lean as it is today, and the release symbols
  # travel separately - their own zip, uploaded as a ci artifact and attached to the release.
  # The diagnostic build DOES carry them, since a debugging binary without symbols is pointless.
  #
  # OPTIONAL because TARGET_PDB_FILE is only written when the linker was asked for one, and
  # MSVC only - the generator expression is not valid for a target that cannot produce a pdb.
  if(MSVC)
    install(FILES "$<TARGET_PDB_FILE:${_tgt}>"
            CONFIGURATIONS RelWithDebInfo Debug
            DESTINATION ${CMAKE_INSTALL_BINDIR} OPTIONAL)
  endif()

  # macOS: the same problem, solved differently. -g leaves the DWARF in the .o files and puts
  # only a map of paths into the binary, so once the build tree is gone the symbols are
  # unreachable - a mac "diagnostic" build with no dsymutil gives you assertions and an
  # unreadable stack. dsymutil collects them into a .dSYM bundle keyed to the binary by uuid.
  #
  # Behind an OPTION rather than a config check because it is expensive: the bundles run to
  # roughly 14x the size of the binary, and doing that on every ordinary build would cost
  # hundreds of megabytes and the time to write them. Only the diagnostic build asks for it.
  #
  # Done at INSTALL time, not POST_BUILD: add_custom_command(TARGET ...) only works in the
  # directory that created the target, and every one of these is created in a subdirectory.
  # Install time is the better moment anyway - the bundle is written straight into the install
  # tree, so cpack picks it up with no separate install(DIRECTORY) rule to keep in step.
  if(APPLE AND PESTPP_MAC_DSYM)
    install(CODE "execute_process(COMMAND dsymutil \"$<TARGET_FILE:${_tgt}>\" -o \"\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/$<TARGET_FILE_NAME:${_tgt}>.dSYM\")")
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
