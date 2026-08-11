# demofit

demofit is an R package for fitting parametric mortality curves and forecasting mortality rates with stochastic models.

## Features

- 14 parametric mortality curves (constrained and unconstrained parameterisations)
- 7 stochastic mortality forecasting models (LC, RH, APC, CBD, M6, M7, STAR)
- Lee–Carter model adjusted for sparse mortality data observed in irregular years
- Multiple optimisation methods for robust parameter estimation of mortality curves
- Ensemble point and interval mortality forecasting
- Multi-population mortality forecasting models (CFM, CAE)
- Unified interfaces through MC(), FCS(), and PFCS()
- S3 methods for obtaining fitting and forecasting results

## Installation

```r
install.packages("demofit")
remotes::install\_github("jackieli/demofit")
```

## Documentation

See the package manual and the package vignette for full details.
