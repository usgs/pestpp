# Building PEST++ with CMake

## Introduction

CMake is a cross-platform tool used to build software on most operating systems and compilers. CMake can be installed using one of several methods, such as from [official sources](https://cmake.org/download/), from Anaconda/Miniconda (`conda install cmake`), or from PyPI (`pip install cmake`). Newest versions are recommended, with a minimum 3.9 enforced for PEST++.

## General guidance

The normal approach to compile a project with CMake is outside of the source tree, in a directory named `build` (or similar). Therefore, before running `cmake`, ensure this is done in an empty directory:

    mkdir build
    cd build
    cmake ..

Variables and options are passed to CMake using `-D<var>=<value>` options to the command-line interface.

* `BUILD_SHARED_LIBS=OFF` - Static libraries are built by default, but shared libraries can be enabled by setting value to `ON`.
* `CMAKE_BUILD_TYPE` - Specify the build type for single-configuration generators, including Makefile and Ninja. Possible values are `Debug`, `Release`, `RelWithDebInfo`, or `MinSizeRel`.
* `CMAKE_INSTALL_PREFIX` - Default is to install executable to the project source tree in `./bin`, but a different prefix (such as `$HOME/.local`) can be specified.
* `CMAKE_CXX_COMPILER` - Normally this is detected, but it can be overridden to use a different C++ compiler.
* `CMAKE_Fortran_COMPILER` - If specified (such as `ifort`), this enables Fortran support to build additional targets.
* `ENABLE_Fortran=OFF` - If set `ON`, enable Fortran support, using default Fortran compiler to build additional targets.

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

Install to desired location, which will normally be in `./bin/linux` or `./bin/mac`:

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

[Ninja](https://ninja-build.org/) is a fast generator that works well with CMake and Visual Studio compilers, and similarly can be installed from [official sources](https://github.com/ninja-build/ninja/releases), from Anaconda/Miniconda (`conda install ninja`), or from PyPI (`pip install cmake`).

Configuring with a specified generator:

    cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..

Build normally takes advantage of as many cores as are available:

    cmake --build .

And install:

    cmake --install .

Build and install steps can be combined with one command:

    cmake --build . --target install

### Visual Studio

Configure using Visual Studio 2019:

    cmake -G"Visual Studio 16 2019" -A Win64 ..

Or configure using Visual Studio 2017:

    cmake -G"Visual Studio 15 2017 Win64" ..

Build all targets with "Release" configuration:

    cmake --build . --config Release

Or choose only one target to compile:

    cmake --build . --config Release --target pestpp-ies

Install executables:

    cmake --build . --config Release --target install

### Intel C++/Fortran

To configure using Intel C++/Fortran, use "NMake Makefiles" generator:

    cmake -G"NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icl -DCMAKE_Fortran_COMPILER=ifort ..

Then follow similar steps to build and install:

    cmake --build . --config Release --target install
