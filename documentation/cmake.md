# Building PEST++ with CMake

## Introduction

CMake is a cross-platform tool used to build software on most operating systems and compilers. CMake can be installed using one of several methods, such as from [official sources](https://cmake.org/download/), from Anaconda/Miniconda (`conda install cmake`), or from PyPI (`pip install cmake`). Newest versions are recommended, with a minimum 3.9 enforced for PEST++.

## General guidance

The normal approach to compile a project with CMake is outside of the source tree, in a directory named `build` (or similar). Therefore, before running `cmake`, ensure this is done in an empty directory (above the base directory):

    mkdir build
    cd build
    cmake ..

Variables and options are passed to CMake using `-D<var>=<value>` options to the command-line interface.

* `BUILD_SHARED_LIBS=OFF` - Static libraries are built by default, but shared libraries can be enabled by setting value to `ON`.
* `CMAKE_BUILD_TYPE` - Specify the build type for single-configuration generators, including Makefile and Ninja. Possible values are `Debug`, `Release`, `RelWithDebInfo`, or `MinSizeRel`.
* `CMAKE_INSTALL_PREFIX` - If CMake is used to install, specify an install prefix, such as `$HOME/.local`.
* `CMAKE_CXX_COMPILER` - Normally this is detected, but it can be overridden to use a different C++ compiler.
* `CMAKE_Fortran_COMPILER` - If specified (such as `ifort`), this enables Fortran support to build additional targets.
* `ENABLE_Fortran=OFF` - If set `ON`, enable Fortran support, using default Fortran compiler to build additional targets.
* `INSTALL_LOCAL=ON` - By default, executables are installed locally in ./bin after they are built, which is handy for testing.

## Linux and macOS

If a generator (`-G`) is not specified, the default for Linux is Makefile, which simple to use.

### Basic usage

Create an empty directory:

    mkdir build && cd build

Configure a release build using the default C++ compiler:

    cmake -DCMAKE_BUILD_TYPE=Release ..

Build all targets with (e.g.) 8 cores:

    make -j8

Or build a single executable, and any dependencies for that target:

    make pestpp-ies

Package the binaries into a tar.gz archive file with:

    cpack -G TGZ

Optionally, install to `CMAKE_INSTALL_PREFIX` (default usually `/usr/local`) with:

    make install

### Configure for different compilers

The configure step can be re-run, but CMake usually requires a "clean" build directory when switching compilers:

    cd .. && rm -rf build && mkdir build && cd build

To configure with Intel C++ and Fortran:

    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort ..

Or Clang++ version 8:

    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-8 ..

## Windows

### Using Ninja

[Ninja](https://ninja-build.org/) is a fast generator that works well on Windows with CMake using either Visual Studio or Intel compilers. Ninja can similarly can be installed from [official sources](https://github.com/ninja-build/ninja/releases), from Anaconda/Miniconda (`conda install ninja`), or from PyPI (`pip install ninja`). If Intel Fortran support is anticipated, it is recommended to install the [Kitware fork of Ninja](https://github.com/Kitware/ninja/releases).

Ninja needs to be run from a command prompt that has been configured for a particular compiler, such as Visual Studio or Intel compilers.

For instance, from a "x64 Native Tools Command Prompt for VS2017", configure with a specified generator:

    cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..

Or to configure for Intel C++/Fortran using an Intel Compiler prompt:

    cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icl -DCMAKE_Fortran_COMPILER=ifort

Build with Ninja normally takes advantage of as many cores as are available:

    ninja

Or specify a number of cores to use, e.g. 5:

    ninja -j5

And package the executables to a ZIP file:

    cpack -G ZIP

### Visual Studio

Visual Studio generators create a project that can be opened with the Visual Studio GUI. As they have multiple configurations, `CMAKE_BUILD_TYPE` does not need to be specified, but `--config`/`-C` needs to be specified for some steps.

Configure using Visual Studio 2019:

    cmake -G"Visual Studio 16 2019" -A x64 ..

Or configure using Visual Studio 2017:

    cmake -G"Visual Studio 15 2017 Win64" ..

Build all targets with "Release" configuration:

    cmake --build . --config Release

If you have CMake version 3.12 or later, parallel builds can be specified using either `-j` to use maximum number of concurrent processes or (e.g.) `-j5` for 5 jobs:

    cmake --build . --config Release -j5

Or choose only one target to compile:

    cmake --build . --config Release --target pestpp-ies

And package the executables to a ZIP file:

    cpack -G ZIP -C Release

### NMake Makefiles

If Ninja is not available, a basic generator that works for both Visual Studio and Intel compilers is to specify "NMake Makefiles". This method is the slowest, as it cannot be parallelized. For example, configure for Intel C++/Fortran:

    cmake -G"NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icl -DCMAKE_Fortran_COMPILER=ifort ..

Then follow similar steps to build:

    cmake --build .

And package the executables to a ZIP file:

    cpack -G ZIP -C Release
