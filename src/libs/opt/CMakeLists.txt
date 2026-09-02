# This CMake file is part of PEST++

# This library is only linked by pestpp-opt, so force it static
# Deliberately STATIC even when BUILD_SHARED_LIBS is ON. This is vendored Coin-OR, and the
# compiled subset references symbols from Clp's command-line driver (saveSolution, mainTest,
# CbcOrClpParam::set*Parameter) that are not built here. A static archive never pulls those
# objects in, so it links fine; a shared library must resolve every symbol and cannot.
# Global PIC means it still links cleanly into the shared pestpp_com.
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

# PESTPP_CXX_WARN_FLAGS is deliberately NOT applied here.
#
# This is vendored Coin-OR. Our warning policy exists to catch defects in code we write and can
# fix, and none of that applies to a third-party library we take as-is: every warning it emits
# is one nobody will ever act on, and with -Werror it stops being noise and starts breaking the
# build outright.
#
# It did exactly that. Coin-OR tripped -Werror=self-assign (three deliberate self-assignments,
# the old trick for silencing an unused-variable warning), and then CoinIndexedVector.cpp
# tripped -Werror=uninitialized on the RelWithDebInfo build - which the Release build had not
# hit, because gcc's analysis changes with the optimisation level. Silencing them one category
# at a time was chasing symptoms; the cause is linting code we do not maintain.
#
# Our own targets keep the full set, which is where it does some good.

# opt is vendored Coin-OR and uses nothing from pestpp_com; this edge only created a
# dependency cycle, which blocks any shared-library build
#target_link_libraries(opt pestpp_com)
