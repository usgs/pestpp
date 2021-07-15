<p align="center">
  <img src="documentation/pestpplogo.png" alt="pestpplogo image">
</p>
# PEST++

## a Software Suite for Parameter Estimation, Uncertainty Analysis, Management Optimization and Sensitivity Analysis

PEST++ is a software suite aimed at supporting complex numerical models in the decision-support context.  Much focus has been devoted to supporting environmental models (groundwater, surface water, etc) but these tools are readily applicable to any computer model.

[![Travis Status](https://travis-ci.org/usgs/pestpp.svg?branch=master)](https://travis-ci.org/usgs/pestpp)
[![Appveyor status](https://ci.appveyor.com/api/projects/status/rqadojcv8bkj5gr0/branch/master?svg=true)](https://ci.appveyor.com/project/jwhite-usgs/pestpp/branch/master)

<br>


## Documentation

The latest official report documenting PEST++ is available from the USGS:

[https://pubs.er.usgs.gov/publication/tm7C26](https://pubs.er.usgs.gov/publication/tm7C26)

##### Suggested Citation:

White, J.T., Hunt, R.J., Fienen, M.N., and Doherty, J.E., 2020, Approaches to Highly Parameterized Inversion: PEST++ Version 5, a Software Suite for Parameter Estimation, Uncertainty Analysis, Management Optimization and Sensitivity Analysis: U.S. Geological Survey Techniques and Methods 7C26, 52 p., [https://doi.org/10.3133/tm7C26](https://doi.org/10.3133/tm7C26).

##### User's Manual

The lastest PEST++ users manual is available [here](https://github.com/usgs/pestpp/tree/develop/documentation). Direct zip download [here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/usgs/pestpp/tree/develop/documentation)

<br>

## Links to latest binaries

As of version 4.3.11, PEST++ pre-compiled binaries for windows and mac are available as a github release.  For older version of PEST++, precompiled binaries are in the `bin` directory

<br>

## Compiling
The develop branch includes a Visual Studio solution, as well as CMake files for cross-compilation on all operating systems.

See details [here](documentation/cmake.md) to compile using CMake.

<br>

## Overview
The PEST++ software suite includes several stand-alone tools for model-independent (non-intrusive) computer model parameter estimation and uncertainty analysis.  Codes include:

* ``pestpp-glm``: deterministic GLM parameter estimation using "on-the-fly" subspace reparameterization, effectively reproducing the SVD-Assist methodology of PEST without any user intervention and FOSM-based parameter and (optional) forecast uncertainty estimation with support for generating posterior parameter realizations.

* ``pestpp-sen``: Global sensitivity analysis using either Morris or Sobol

* ``pestpp-swp``: a generic parallel run utility driven by a CSV file of parameter values

* ``pestpp-opt``: chance-constrainted linear programming

* ``pestpp-ies``: iterative ensemble smoother implementation of GLM (based on the work Chen and Oliver 2013) with support for generic localization (local analysis and/or covariance localization)

* ``pestpp-pso``: particle-swarm based inversion.

All members of the software suite can be compiled for PC, MAC, or Linux and have several run managers to support parallelization.  precompiled binaries are available in the "bin" folder.  Windows users with older OS versions should use the ``iwin`` binaries (starting "i", compiled with intel C++) to avoid the dreaded MSVC missing runtime DLL issue.

<br>

## Funding

Funding for PEST++ has been provided by the U.S. Geologial Survey. The New Zealand Strategic Science Investment Fund as part of GNS Science’s (https://www.gns.cri.nz/) Groundwater Research Programme has also funded contributions 2018-present.  Intera, Inc. also provides ongoing support for PEST++.

<br>

## Recent developements

PEST++ version 5 has been released!  Please see the users manual for current input instructions and options for all PEST++ tools.  Also, several new tools are in development, include PESTPP-DA (generalized data assimilation including iterative ensemble Kalman filter), PESTPP-MOU (single and multiple objective constrained optimization under uncertainty) and PESTPP-SQP (ensemble-based constrainted sequential quadratic programming under uncertainty).  If you would like to be an early adopter/beta tester, please let us know!

The PEST++ suite has been refactored to remove the fortran dependancy.  This includes the template and instruction file processing routines.

An updated control file format has been implemented.  All the existing PEST++ tools are still backward compatible with the standard control file format.  However, new tools that are in developement to support sequential data assimilation need the external file format.  Additional, new functionality is being developed that will require the updated control file format.  The new control file format is described in documentation and pyEMU can operate on this new format interchangeably.

<br>

## PEST++ References:

White, J.T., Hunt, R.J., Fienen, M.N., and Doherty, J.E., 2020, Approaches to Highly Parameterized Inversion: PEST++ Version 5, a Software Suite for Parameter Estimation, Uncertainty Analysis, Management Optimization and Sensitivity Analysis: U.S. Geological Survey Techniques and Methods 7C26, 52 p., https://doi.org/10.3133/tm7C26.

White, J. T., 2018, A model-independent iterative ensemble smoother for efficient history-matching and uncertainty quantification in very high dimensions. Environmental Modelling & Software. 109. 10.1016/j.envsoft.2018.06.009. <a ref="http://dx.doi.org/10.1016/j.envsoft.2018.06.009">http://dx.doi.org/10.1016/j.envsoft.2018.06.009</a>.

White, J. T., Fienen, M. N., Barlow, P. M., and Welter, D.E., 2017, A tool for efficient, model-independent management optimization under uncertainty. Environmental Modeling and Software.  <a ref="http://dx.doi.org/10.1016/j.envsoft.2017.11.019">http://dx.doi.org/10.1016/j.envsoft.2017.11.019</a>.

Welter, D.E., White, J.T., Hunt, R.J., and Doherty, J.E., 2015, Approaches in highly parameterized inversion— PEST++ Version 3, a Parameter ESTimation and uncertainty analysis software suite optimized for large environmental models: U.S. Geological Survey Techniques and Methods, book 7, chap. C12, 54 p., <a ref="http://dx.doi.org/10.3133/tm7C12">http://dx.doi.org/10.3133/tm7C12</a>.

Welter, D.E., Doherty, J.E., Hunt, R.J., Muffels, C.T., Tonkin, M.J., and Schreüder, W.A., 2012, Approaches in highly parameterized inversion—PEST++, a Parameter ESTimation code optimized for large environmental models: U.S. Geological Survey Techniques and Methods, book 7, section C5, 47 p., available only at <a ref="http://pubs.usgs.gov/tm/tm7c5">http://pubs.usgs.gov/tm/tm7c5</a>.

<br>

### Related Links:

* <a ref="https://www.usgs.gov/software/pest-parameter-estimation-code-optimized-large-environmental-models">https://www.usgs.gov/software/pest-parameter-estimation-code-optimized-large-environmental-models </a>
* <a ref="http://www.pesthomepage.org">http://www.pesthomepage.org </a>
* <a ref="https://github.com/pypest/pyemu">https://github.com/pypest/pyemu </a>

<br>

## Testing

The ``benchmarks`` folder contains a simple worked example and basic testing routines that are used for basic CI testing.  Many full-worked test problems of varying problem sizes are now located in separate repos:

* [pestpp-glm benchmarks](https://github.com/usgs/pestpp-glm_benchmarks)
* [pestpp-ies benchmarks](https://github.com/pestpp/pestpp-ies_benchmarks)
* [pestpp-opt benchmarks](https://github.com/pestpp/pestpp-opt_benchmarks)

<br>

## Dependencies

Much work has been done to avoid additional external dependencies in PEST++.  As currently designed, the project is fully self-contained.  

<br>

# optional ``++`` arguments

please see the PEST++ users manual in the ``documentation`` directory for a current and complete description of all ``++`` options

<br>

###### USGS disclaimer

This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use
