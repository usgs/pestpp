# This CMake file is part of PEST++

# This library depends on rm_abstract; avoid cyclic dependency
add_library(pestpp_com STATIC
  constraints.cpp
  covariance.cpp
  DifferentialEvolution.cpp
  eigen_tools.cpp
  Ensemble.cpp
  EnsembleMethodUtils.cpp
  EnsembleSmoother.cpp
  FileManager.cpp
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
  PriorInformation.cpp
  QSqrtMatrix.cpp
  Regularization.cpp
  RestartController.cpp
  SVDASolver.cpp
  SVDPackage.cpp
  SVD_PROPACK.cpp
  SVDSolver.cpp
  TerminationController.cpp
  Transformation.cpp
  SQP.cpp
  MOEA.cpp
  sequential_lp.cpp
)

target_include_directories(pestpp_com INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}" /src/libs/opt)

target_compile_options(pestpp_com PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(pestpp_com
  common
  rm_abstract
  opt
)

if(BUILD_SHARED_LIBS)
  set_property(TARGET pestpp_com PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
