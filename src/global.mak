# Global options for GNU Make to build PEST++ libraries and programs
#
# Choose one of each:
# SYSTEM ?= linux (default)
#        ?= mac
#        ?= win
# COMPILER ?= gcc (default, if ifort not available)
#          ?= intel
# STATIC ?= -static (default)
#        ?= no (shared dynamic linking)
# These can be kept in local.mak
-include $(top_builddir)/local.mak

ifndef SYSTEM  # then autodetect this
    ifeq ($(OS),Windows_NT)
        SYSTEM ?= win
    else
        UNAME_S := $(shell uname -s)
        ifeq ($(UNAME_S),Linux)
            SYSTEM ?= linux
        endif
        ifeq ($(UNAME_S),Darwin)
            SYSTEM ?= mac
        endif
    endif
endif


ifeq ($(SYSTEM),mac)  # macOS
    bindir ?= $(top_builddir)/../bin/mac/
else ifeq ($(SYSTEM),linux)  # GNU Linux
    bindir ?= $(top_builddir)/../bin/linux/
else ifeq ($(SYSTEM),win)  # Microsoft Windows
    bindir ?= $(top_builddir)/../bin/windows/
else
    $(error SYSTEM not understood: $(SYSTEM). Use one of mac, linux or win.)
endif

ifndef COMPILER  # also autodetect this
    # If ifort on PATH, then assume intel
    ifeq ($(SYSTEM),win)
        IFORT := $(shell where ifort 2> nul)
    else
        IFORT := $(shell which ifort 2> /dev/null)
    endif
    ifneq ($(IFORT),)
        COMPILER := intel
    else
        COMPILER := gcc
    endif
endif

# Determine static (default '-static') or shared dynamic linking
STATIC ?= -static
ifndef STATIC
    STATIC = no
endif

ifeq ($(SYSTEM),win)  # Microsoft Windows
    EXE_EXT = .exe
    OBJ_EXT = .obj
    LIB_EXT = .lib
    CP = copy
    RM = del /Q
    MKDIR = md
else  # POSIX (mac, linux)
    OBJ_EXT = .o
    LIB_PRE = lib
    LIB_EXT = .a
    CP = cp
    MKDIR = mkdir -p
endif

# Default optimization
OPT_FLAGS ?= -O2

ifeq ($(COMPILER),intel)  # Intel compilers
    ifeq ($(SYSTEM),win)
        # Warning: this build method is not well tested
        CXX = icl
        OPT_FLAGS ?= /nologo /Qmkl:sequential
        CXXFLAGS ?= $(OPT_FLAGS) /Qstd=c++11 /EHsc
        FFLAGS ?= $(OPT_FLAGS) /fpp
        FFREE = /free
    else  # mac,linux
        CXX = icpc
        OPT_FLAGS ?= -mkl=sequential
        CXXFLAGS ?= $(OPT_FLAGS) -std=c++11
        FFLAGS ?= $(OPT_FLAGS) -fpp
        FFREE = -free
    endif
    FC = ifort

    ifeq ($(SYSTEM),win)
        EXT_INCLUDES = -I"$(MKLROOT)"\include
        EXT_LIBS = \
            mkl_blas95_lp64.lib \
            mkl_lapack95_lp64.lib
        ifeq ($(STATIC),no)  # dynamic linking
            EXT_LIBS += \
                mkl_intel_lp64_dll.lib \
                mkl_sequential_dll.lib \
                mkl_core_dll.lib
        else  # static linking
            EXT_LIBS += \
                mkl_intel_lp64.lib \
                mkl_sequential.lib \
                mkl_core.lib
        endif
    else ifeq ($(SYSTEM),linux)
        MKLROOT ?= /opt/intel/compilers_and_libraries/linux/mkl
        EXT_INCLUDES = -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
        EXT_LIBS = \
            ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a \
            ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a
        ifeq ($(STATIC),no)  # dynamic linking
            EXT_LIBS += \
                -L${MKLROOT}/lib/intel64 \
                -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
        else  # static linking
            EXT_LIBS += \
                -Wl,--start-group \
                    ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
                    ${MKLROOT}/lib/intel64/libmkl_sequential.a \
                    ${MKLROOT}/lib/intel64/libmkl_core.a \
                -Wl,--end-group \
                -lifport
        endif
        EXT_LIBS += -lifcore -lpthread -lm -ldl
    else ifeq ($(SYSTEM),mac)
        MKLROOT ?= /opt/intel/compilers_and_libraries_2019.4.233/mac/mkl
        EXTRADIR = /opt/intel/compilers_and_libraries_2019.4.233/mac/compiler
        EXT_INCLUDES = -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
        EXT_LIBS = \
             ${MKLROOT}/lib/libmkl_lapack95_ilp64.a \
            ${MKLROOT}/lib/libmkl_blas95_ilp64.a
        ifeq ($(STATIC),no)  # dynamic linking
            EXT_LIBS += \
                -L${MKLROOT}/lib \
                -Wl,-rpath,${MKLROOT}/lib \
                -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
        else  # static linking
            EXT_LIBS += \
                ${MKLROOT}/lib/libmkl_intel_ilp64.a \
                ${MKLROOT}/lib/libmkl_sequential.a \
                ${MKLROOT}/lib/libmkl_core.a \
            ${EXTRADIR}/lib/libifcore.a
        endif
        EXT_LIBS += -lpthread -lm -ldl
    endif  # $(SYSTEM)
else  # $(COMPILER))
    # Assume GNU Compiler Collection
    ifeq ($(COMPILER),gcc)
        CXX = g++
        FC = gfortran
    else
        CXX ?= g++
        FC ?= gfortran
    endif
    CXXFLAGS ?= $(OPT_FLAGS) -std=c++11
    FFLAGS ?= $(OPT_FLAGS) -cpp
    FFREE = -ffree-form
    EXT_LIBS = -lpthread -llapack -lblas -lgfortran -lquadmath
# else
#     $(error COMPILER not understood: $(COMPILER). Use one of intel or gcc.)
endif  # $(COMPILER)

# Assume linker is the C++ compiler
LD = $(CXX)
LDFLAGS += -pthread
ifneq ($(STATIC),no)
    LDFLAGS += $(STATIC)
endif

# r=insert with replacement; c=create archive; s=add index
ARFLAGS := rcs

LIBS_DIR := $(top_builddir)/libs

PESTPP_INCLUDES := \
    -I $(LIBS_DIR)/Eigen \
    -I $(LIBS_DIR)/common \
    -I $(LIBS_DIR)/run_managers/abstract_base \
    -I $(LIBS_DIR)/run_managers/yamr \
    -I $(LIBS_DIR)/run_managers/serial \
    -I $(LIBS_DIR)/run_managers/external \
    -I $(LIBS_DIR)/run_managers/wrappers \
    -I $(LIBS_DIR)/pestpp_common \
    -I $(LIBS_DIR)/opt $(EXT_INCLUDES)

# Be careful with the order of library dependencies
PESTPP_LIBS := \
    -L$(LIBS_DIR)/pestpp_common -lpestpp_com \
    -L$(LIBS_DIR)/run_managers/wrappers -lrm_wrappers \
    -L$(LIBS_DIR)/run_managers/yamr -lrm_yamr \
    -L$(LIBS_DIR)/run_managers/serial -lrm_serial \
    -L$(LIBS_DIR)/run_managers/external -lrm_external \
    -L$(LIBS_DIR)/run_managers/abstract_base -lrm_abstract \
    -L$(LIBS_DIR)/mio -lmio \
    -L$(LIBS_DIR)/common -lcommon \
    -L$(LIBS_DIR)/propack -lpropack \
    -L$(LIBS_DIR)/opt -lopt \
     $(EXT_LIBS)


# Generic pattern rules

%$(OBJ_EXT): %.cpp
	$(CXX) -c $(PESTPP_INCLUDES) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%$(OBJ_EXT): %.f
	$(FC) -c $(FFLAGS) -o $@ $<

%$(OBJ_EXT): %.for
	$(FC) -c $(FFLAGS) -o $@ $<

%$(OBJ_EXT): %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%$(OBJ_EXT): %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $<

%$(OBJ_EXT): %.FOR
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $<

.SUFFIXES: .c .h .cpp .hpp .f .for .F .FOR .mod $(OBJ_EXT) $(LIB_EXT)

# Default rules for handling subdirectories

%-target:
	$(MAKE) -C $*

%-install:
	$(MAKE) -C $* install

%-clean:
	$(MAKE) -C $* clean
