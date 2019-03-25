# PEST++

## Tools for non-intrusive and scalable parameter estimation and uncertainty quantification

PEST++ is a software suite aimed at supporting complex numerical models in the decision-support context.  Much focus has been devoted to supporting environmental models (groundwater, surface water, etc) but these tools are readily applicable to any computer model.

<br><br><br>

[![Travis Status](https://travis-ci.org/jtwhite79/pestpp.svg?branch=master)](https://travis-ci.org/jtwhite79/pestpp)
[![Build status](https://ci.appveyor.com/api/projects/status/jlyivnw4jhstp8l8?svg=true)](https://ci.appveyor.com/project/jtwhite79/pestpp)

## Documentation

The lastest PEST++ manual is available [here](https://github.com/jtwhite79/pestpp/tree/develop/documentation). Direct zip download [here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/jtwhite79/pestpp/tree/develop/documentation)

## Links to latest binaries

* [windows (users with current visual studio installed)](https://github.com/jtwhite79/pestpp/tree/master/bin/win).  Direct zip download [here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/jtwhite79/pestpp/tree/master/bin/win)
* [windows compiled with intel C++ (the 'i' prefix)](https://github.com/jtwhite79/pestpp/tree/master/bin/iwin).  Direct zip download [here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/jtwhite79/pestpp/tree/master/bin/iwin)
* [mac OS](https://github.com/jtwhite79/pestpp/tree/master/bin/mac).  Direct zip download [here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/jtwhite79/pestpp/tree/master/bin/mac)

## Compiling
The master branch includes a Visual Studio 2017 project, as well as makefiles for linux and mac. The suite has been succcessfully compiled with gcc (g++ and gfortran) 4,5,6 and 7 on ubuntu, fedora, centos, and on slurm/MPI clusters - to use gcc, you must C++11 support (if using gcc 4, only gcc4.9 has this) and you need to have both ``lapack`` and ``blas`` libraries available in the path.  Then, in the `src` directory:

`>>>make clean`

`>>>STATIC=no make install`


this should put the compiled binaries in the `bin/linux` directory.  See the `.travis.yml` file for an example

## Overview
The PEST++ software suite includes several stand-alone tools for model-independent (non-intrusive) computer model parameter estimation and uncertainty analysis.  Codes include:

* ``pestpp-glm``: deterministic GLM parameter estimation using "on-the-fly" subspace reparameterization, effectively reproducing the SVD-Assist methodology of PEST without any user intervention and FOSM-based parameter and (optional) forecast uncertainty estimation with support for generating posterior parameter realizations.

* ``pestpp-sen``: Global sensitivity analysis using either Morris or Sobol

* ``pestpp-swp``: a generic parallel run utility driven by a CSV file of parameter values

* ``pestpp-opt``: chance-constrainted linear programming

* ``pestpp-ies``: iterative ensemble smoother implementation of GLM (based on the work Chen and Oliver 2013) with support for generic localization (local analysis and/or covariance localization)

All members of the software suite can be compiled for PC, MAC, or Linux and have several run managers to support parallelization.  precompiled binaries are available in the "bin" folder.  Windows users with older OS versions should use the ``bin/iwin`` binaries (starting "i", compiled with intel C++) to avoid the dreaded MSVC missing runtime DLL issue

## Recent Updates

<b> update 20 March 2019 </b>: PEST++ version 4.2.4 is now available.  ``pestpp-ies`` now implements combined automatic adaptive localization with use of a localization matrix - please see the manual for a description of this functionality.

<b> update 26 February 2019 </b>: PEST++ version 4.2 is now available.  ``pestpp-ies`` now implements automatic adaptive localization, building on the methodology of Luo (2018).  All users need to supply to activate this functionality is set ``++ies_autoadaloc(true)`` (with an optional ``++ies_autoadaloc_sigma_dist()``) - please see the manual for a description of this functionality and these args.  Also, ``pestpp`` has been renamed ``pestpp-glm`` and ``pestpp-gsa`` has been renamed ``pestpp-sen`` to clarify what these tools do (and also in anticipation of some enhancements to the sensitivity analysis tool).  


<b> update 12 November 2018 </b>: ``pestpp-ies`` now supports localization available as combined local analysis/covariance localization.  This is controlled by a localization matrix which lists adjustable parameters and/or parameter groups as columns and non-zero weighted observations and/or observation groups as rows.  The format of this matrix can be PEST ascii, PEST binary or CSV.  

<b> update 3 November 2018 </b>: The FOSM calculations in ``pestpp`` now also support generation and evaluation of a posterior parameter ensemble.  If you add ``++num_reals(50)`` to the end of the pest control file, once the GLM parameter estimation process of ``pestpp`` is complete, during the FOSM parameter uncertainty calculations, 50 posterior parameter realizations are generated from the Schur-based (bayes linear) posterior parameter covariance matrix.  These realizations are saved to a CSV file and are also evaluated (in parallel is ``pestpp`` is running in parallel).

<b> update 3 November 2018 </b>: The repo structure has been significantly refactored in an effort to reduce the size the primary code and binaries.  The main culprit was the ``benchmarks`` directory, so the various benchmarks have now been split off into separate repos: [https://github.com/jtwhite79/pestpp_benchmarks](https://github.com/jtwhite79/pestpp_benchmarks), [https://github.com/jtwhite79/pestpp-ies_benchmarks](https://github.com/jtwhite79/pestpp-ies_benchmarks), and [https://github.com/jtwhite79/pestpp-opt_benchmarks](https://github.com/jtwhite79/pestpp-opt_benchmarks).  This makes the main repo (this repo) much smaller; the CI testing continues as before by using some tricks for automatically triggering downstream/dependent builds in travis - a commit against this repo will trigger a fake commit to each of those benchmark repos and subsuquently CI testing in travis and appveyor.

<b> update 15 September 2018 </b>: An official PEST++ V4 manual is now available in the Documentation directory.  It is a docx file so any contributions (revisions/extensions/clarifications/etc) are greatly appreciated!

<b> update 5 August 2018 </b>: Welcome to the active fork of ``pestpp-ies``, ``pestpp-opt``, and ``pestpp-swp``. I strive to actively support users, so please raise issues related to these codes as needed.  I will also try to support ``pestpp`` users are best I can although this code is a lower priority for me.

<b> update 30 July 2018 </b>: The PEST++ now supports parameter and observation names up to 200 characters in length. This allows for more descriptive naming and better support for problems with very large numbers of pars and obs.  Note that if a parameter name exceeds 12 chars or an observation name exceeds 20 chars, the resulting jacobian binary file will include truncated names (at 12 and 20 chars, respectively) and names exceeding these lengths are not backward compatible with PEST. 

<b> update 4 July 2018 </b>: PEST++ version 4.0.0 has been released to support the newly-developed ``pestpp-ies``. A manuscript documenting ``pestpp-ies`` is available here: [https://www.sciencedirect.com/science/article/pii/S1364815218302676](https://www.sciencedirect.com/science/article/pii/S1364815218302676).  Stay tuned for an actual manual to accompany version 4!

<b> update 2 May 2018 </b>: some refactoring is underway.  ``sweep`` has been renamed ``pestpp-swp`` and ``gsa`` has been renamed ``pestpp-gsa``.  Also, the initial version of the new iterative ensemble smoother is avaiable as ``pestpp-ies``.  The basic ``++`` options needed for fine-grained control of ``pestpp-ies`` are listed below.   

<b> update 09/20/2017</b>: the new optimization under uncertainty tool is ready!  A supporting publication is in the works and should be available soon (a link will be posted once it is accepted).  This new tool uses the same control file/template file/instruction file approach as other PEST(++) applications, so applying this tool to your problem should be seamless.  Optional "++" args for tool are available further done this page.

<b>update 01/25/2017</b>: intel C++ builds are avaiable for mac and for windows.  For mac users, these are statically-linked so they do not require compilers to be installed.  For windows users, the intel build circumvents the "missing VCOMP140.DLL" error.  Note the intel windows builds are currently in the ~~``intel_c_windows`` branch~~ ``bin/iwin/`` folder.

<b>update 11/25/2016</b>: PEST++ version 3.6 is now available. Some of the many enhancements available in 3.6 include:

* a new approach to implementing regularization. Rather than using the standard pest control file parameters such as ``phimlim``, ``fracphim``, etc, we now offer a single pest++ argument, ``++reg_frac()``, that allows users to specify what fraction  of the composite objective function should be regularization penalty. For example, ``++reg_frac(0.5)`` would result in equal parts data misfit and regularization penalty, which results in the *maximum a posteriori* (MAP) parameter estimate. Using ``++reg_frac()`` will result in substantial speed ups during the lambda calculation process

* a new program for sequential linear programming under uncertainty.  ``pestpp-opt`` is a new executable in the PEST++ suite that uses the standard PEST model independent interface to solve a (sequential) linear programming (LP) problem.  ``pestpp-opt`` relies on the COIN-OR Linear Programming (CLP) solver [https://projects.coin-or.org/Clp]((https://projects.coin-or.org/Clp)).  Also, users have the option to use FOSM-based uncertainty estimation in the evaluation of model-based constraints (such as water levels, stream flows, stream flow depletion, etc) so that a risk-based optimal solution can be found.  See below for the required and optional ``++`` arguments needed to apply ``pestpp-opt``.  Two example problems using the ``pestpp-opt`` tool have been added to the ``benchmarks`` dir.  A publication about this tool is in the works.

* global optimization with differential evolution.  We now have a fully-parallel global solver that implements the differential evolution algorithm (DE) integrated into the pest++ executable.  See below for required and optional ``++`` arguments needed to use the DE solver.

* a new randomization-based SVD solver using the implementation of [https://github.com/ntessore/redsvd-h](https://github.com/ntessore/redsvd-h).  This solver is activated using ``++svd_pack(redsvd)``.  Testing shows it to be very efficient for problems for a wide range of problem sizes, especially with judicious use of ``max_sing``.

* upgrade parameter covariance scaling.  Through the ``++parcov_scale_fac()``, pest++ can now scale the normal matrix (J^tQJ) by a user specified parameter covariance matrix.  If no ``++parcov_filename()`` is provided, pest++ will construct a diagonal parameter covariance matrix from the parameter bounds.  This is a relatively new option and needs more testing, but limited testing to date shows that upgrade vectors resulting from a covariance-scaled normal matrix are more in harmony with expert knowledge.

<b>Update 05/26/2016</b>: PEST++ V3 has been officially released.  It supports a number of really cool features, including global sensitivity analyses, and automatic Bayes linear (first-order, second-moment) parameter and forecast uncertainty estimates.  We also have a utility for fully-parallel parametric sweeps from csv-based parameter files, which is useful for Monte Carlo, design of experiments, surrogate construction, etc.  All of these tools are based on the model-independent communication framework of PEST, so if you have a problem already setup, these tools are ready for you!

<b>Update 10/1/2014</b>: recent stable versions of PEST++ implement dynamic regularization, full restart capabilities, additional options for formulating the normal equations, and an iterative SVD algorithm for very-large problems.  Additionally the YAMR run manager has been improved to use threaded workers so that the master can more easily load balance.  

## PEST++ References:

White, J. T., 2018, A model-independent iterative ensemble smoother for efficient history-matching and uncertainty quantification in very high dimensions. Environmental Modelling & Software. 109. 10.1016/j.envsoft.2018.06.009. <a ref="http://dx.doi.org/10.1016/j.envsoft.2018.06.009">http://dx.doi.org/10.1016/j.envsoft.2018.06.009</a>.

White, J. T., Fienen, M. N., Barlow, P. M., and Welter, D.E., 2017, A tool for efficient, model-independent management optimization under uncertainty. Environmental Modeling and Software.  <a ref="http://dx.doi.org/10.1016/j.envsoft.2017.11.019">http://dx.doi.org/10.1016/j.envsoft.2017.11.019</a>.

Welter, D.E., White, J.T., Hunt, R.J., and Doherty, J.E., 2015, Approaches in highly parameterized inversion— PEST++ Version 3, a Parameter ESTimation and uncertainty analysis software suite optimized for large environmental models: U.S. Geological Survey Techniques and Methods, book 7, chap. C12, 54 p., <a ref="http://dx.doi.org/10.3133/tm7C12">http://dx.doi.org/10.3133/tm7C12</a>.

Welter, D.E., Doherty, J.E., Hunt, R.J., Muffels, C.T., Tonkin, M.J., and Schreüder, W.A., 2012, Approaches in highly parameterized inversion—PEST++, a Parameter ESTimation code optimized for large environmental models: U.S. Geological Survey Techniques and Methods, book 7, section C5, 47 p., available only at <a ref="http://pubs.usgs.gov/tm/tm7c5">http://pubs.usgs.gov/tm/tm7c5</a>.

### Related Links:

* <a ref="http://www.pesthomepage.org">http://www.pesthomepage.org </a>
* <a ref="https://github.com/jtwhite79/pyemu">https://github.com/jtwhite79/pyemu </a>



## Testing
The ``benchmarks`` folder contains a simple worked example that is used for basic CI testing.  Many full-worked test problems of varying problem sizes are now located in separate repos:

* [https://github.com/jtwhite79/pestpp_benchmarks](https://github.com/jtwhite79/pestpp_benchmarks)
* [https://github.com/jtwhite79/pestpp-ies_benchmarks](https://github.com/jtwhite79/pestpp-ies_benchmarks)
* [https://github.com/jtwhite79/pestpp-opt_benchmarks](https://github.com/jtwhite79/pestpp-opt_benchmarks)

## Dependencies
Much work has been done to avoid additional external dependencies in PEST++.  As currently designed, the project is fully self-contained and statically linked. ``lapack`` and ``blas`` are also required - these are included with the intel fortran compiler; those using gcc will need to have these libraries available. 

# optional ``++`` arguments

## please see the PEST++ version 4 manual in the ``documentation`` directory for a current and complete description of all ``++`` options

### USGS disclaimer

This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use
