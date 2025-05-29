# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

LIB := $(LIB_PRE)pestpp_com$(LIB_EXT)
OBJECTS := \
    FileManager \
    ParamTransformSeq \
    SVDPackage \
    Jacobian_1to1 \
    Pest \
    SVD_PROPACK \
    Jacobian \
    pest_data_structs \
    SVDSolver \
    eigen_tools \
    PriorInformation \
    TerminationController \
    ModelRunPP \
    QSqrtMatrix \
    Transformation \
    ObjectiveFunc \
    Regularization \
    RestartController \
    OutputFileWriter \
    PerformanceLog \
    SVDASolver \
    DifferentialEvolution \
    Ensemble \
    EnsembleSmoother \
    Localizer\
    logger\
    linear_analysis\
    covariance \
    constraints \
    EnsembleMethodUtils \
    sequential_lp
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
