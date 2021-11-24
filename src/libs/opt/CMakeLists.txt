# This CMake file is part of PEST++

# This library is only linked by pestpp-opt, so force it static
add_library(opt STATIC
  CbcOrClpParam.cpp
  Clp_ampl.cpp
  ClpCholeskyBase.cpp
  ClpCholeskyDense.cpp
  ClpCholeskyTaucs.cpp
  ClpCholeskyWssmp.cpp
  ClpCholeskyWssmpKKT.cpp
  Clp_C_Interface.cpp
  ClpConstraint.cpp
  ClpConstraintLinear.cpp
  ClpConstraintQuadratic.cpp
  ClpDualRowDantzig.cpp
  ClpDualRowPivot.cpp
  ClpDualRowSteepest.cpp
  ClpDummyMatrix.cpp
  ClpDynamicExampleMatrix.cpp
  ClpDynamicMatrix.cpp
  ClpEventHandler.cpp
  ClpFactorization.cpp
  ClpGubDynamicMatrix.cpp
  ClpGubMatrix.cpp
  ClpHelperFunctions.cpp
  ClpInterior.cpp
  ClpLinearObjective.cpp
  ClpLsqr.cpp
  # ClpMain.cpp -- remove for building library
  ClpMatrixBase.cpp
  ClpMessage.cpp
  ClpModel.cpp
  ClpNetworkBasis.cpp
  ClpNetworkMatrix.cpp
  ClpNode.cpp
  ClpNonLinearCost.cpp
  ClpObjective.cpp
  ClpPackedMatrix.cpp
  ClpPdcoBase.cpp
  ClpPdco.cpp
  ClpPlusMinusOneMatrix.cpp
  ClpPredictorCorrector.cpp
  ClpPresolve.cpp
  ClpPrimalColumnDantzig.cpp
  ClpPrimalColumnPivot.cpp
  ClpPrimalColumnSteepest.cpp
  ClpQuadraticObjective.cpp
  ClpSimplex.cpp
  ClpSimplexDual.cpp
  ClpSimplexNonlinear.cpp
  ClpSimplexOther.cpp
  ClpSimplexPrimal.cpp
  ClpSolve.cpp
  ClpSolver.cpp
  CoinAbcBaseFactorization1.cpp
  CoinAbcBaseFactorization2.cpp
  CoinAbcBaseFactorization3.cpp
  CoinAbcBaseFactorization4.cpp
  CoinAbcBaseFactorization5.cpp
  CoinAlloc.cpp
  CoinBuild.cpp
  CoinDenseFactorization.cpp
  CoinDenseVector.cpp
  CoinError.cpp
  CoinFactorization1.cpp
  CoinFactorization2.cpp
  CoinFactorization3.cpp
  CoinFactorization4.cpp
  CoinFileIO.cpp
  CoinFinite.cpp
  CoinIndexedVector.cpp
  CoinLpIO.cpp
  CoinMessage.cpp
  CoinMessageHandler.cpp
  CoinModel.cpp
  CoinModelUseful2.cpp
  CoinModelUseful.cpp
  CoinMpsIO.cpp
  CoinOslFactorization2.cpp
  CoinOslFactorization3.cpp
  CoinOslFactorization.cpp
  CoinPackedMatrix.cpp
  CoinPackedVectorBase.cpp
  CoinPackedVector.cpp
  CoinParam.cpp
  CoinParamUtils.cpp
  CoinPostsolveMatrix.cpp
  CoinPrePostsolveMatrix.cpp
  CoinPresolveDoubleton.cpp
  CoinPresolveDual.cpp
  CoinPresolveDupcol.cpp
  CoinPresolveEmpty.cpp
  CoinPresolveFixed.cpp
  CoinPresolveForcing.cpp
  CoinPresolveHelperFunctions.cpp
  CoinPresolveImpliedFree.cpp
  CoinPresolveIsolated.cpp
  CoinPresolveMatrix.cpp
  CoinPresolveMonitor.cpp
  CoinPresolvePsdebug.cpp
  CoinPresolveSingleton.cpp
  CoinPresolveSubst.cpp
  CoinPresolveTighten.cpp
  CoinPresolveTripleton.cpp
  CoinPresolveUseless.cpp
  CoinPresolveZeros.cpp
  CoinRational.cpp
  CoinSearchTree.cpp
  CoinShallowPackedVector.cpp
  CoinSimpFactorization.cpp
  CoinSnapshot.cpp
  CoinStructuredModel.cpp
  CoinWarmStartBasis.cpp
  CoinWarmStartDual.cpp
  CoinWarmStartPrimalDual.cpp
  CoinWarmStartVector.cpp
  Idiot.cpp
  IdiSolve.cpp
)

target_include_directories(opt INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_definitions(opt PRIVATE
  HAVE_CMATH
  HAVE_FLOAT_H
  USE_EKKWSSMP
)

if(WIN32)
  target_compile_definitions(opt PRIVATE
    CLP_BUILD
    _CRT_SECURE_NO_WARNINGS
  )
endif()

target_compile_options(opt PRIVATE ${PESTPP_CXX_WARN_FLAGS})

target_link_libraries(opt pestpp_com)
