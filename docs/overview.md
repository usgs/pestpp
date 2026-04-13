# Overview

PEST++ is a parameter estimation and uncertainty analysis tool suite that includes:

- **pestpp-glm** — Gauss-Levenberg-Marquardt gradient-based parameter estimation
- **pestpp-ies** — Iterative Ensemble Smoother for uncertainty analysis
- **pestpp-da** — Data Assimilation
- **pestpp-opt** — Chance-constrained linear/non-linear optimization
- **pestpp-mou** — Multi-Objective Uncertainty analysis (evolutionary algorithm)
- **pestpp-sqp** — Sequential Quadratic Programming
- **pestpp-swp** — Parameter sweep utility
- **pestpp-gsa** — Global Sensitivity Analysis (Morris, Sobol, Tornado)

## Source Layout

```
src/
  libs/common/          # Core utilities, networking, error handling
  libs/pestpp_common/   # Shared PEST++ algorithms and data structures
  libs/run_managers/    # Distributed/parallel run management (Panther)
  libs/opt/             # LP/optimization helpers
  programs/             # Individual executable entry points
```

## Building

```bash
pixi run build-release
```

See the [users manual](https://github.com/usgs/pestpp/blob/master/documentation/pestpp_users_manual.md) for full usage documentation.
