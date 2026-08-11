# This CMake file is part of PEST++

# This library depends on rm_abstract; avoid cyclic dependency
add_library(pestpp_com
  constraints.cpp
  covariance.cpp
  eigen_tools.cpp
  Ensemble.cpp
  EnsembleMethodUtils.cpp
  EnsembleSmoother.cpp
  FileManager.cpp
  GLM.cpp
  Jacobian_1to1.cpp
  Jacobian.cpp
  linear_analysis.cpp
  Localizer.cpp
  logger.cpp
  ModelRunPP.cpp
  ObjectiveFunc.cpp
  OutputFileWriter.cpp
  ParamTransformSeq.cpp
  PerformanceLog.cpp
  Pest.cpp
  pest_data_structs.cpp
  pest_options_registry.cpp
  PriorInformation.cpp
  QSqrtMatrix.cpp
  Regularization.cpp
  RestartController.cpp
  SVDPackage.cpp
  SVDSolver.cpp
  TerminationController.cpp
  Transformation.cpp
  SQP.cpp
  MOEA.cpp
  sequential_lp.cpp
  DataAssimilator.cpp
  DaCycleDriver.cpp
)

# SVD_PROPACK.cpp is deliberately NOT built. It needs the Fortran PROPACK routines
# (dlamch_, dlansvd_sparce_), and nothing constructs SVD_PROPACK anyway - ++SVD_PACK(PROPACK)
# is deprecated and silently falls back to REDSVD. A static archive hid this because the
# object was never pulled in; a shared library must resolve every symbol.

target_include_directories(pestpp_com INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}" /src/libs/opt)

target_compile_options(pestpp_com PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp_com
  common
  rm_abstract
  opt
  rm_yamr
  rm_serial
)

# PIC is set globally now (CMAKE_POSITION_INDEPENDENT_CODE in the top-level CMakeLists)
