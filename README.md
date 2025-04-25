# MVmixture
An R package to implement the methods in McGee & Antonelli, found here: https://doi.org/10.48550/arXiv.2504.17195

The package can be used to fit (possibly multivariate) distributed lag models using the `DLM=TRUE` argument, with clustering across outcomes and exposures. It can also apply to other exposure structures (like exposures measured via different biomarkers).

It can also be used to fit adaptive (possibly multivariate) multiple index models with unknown index structure using the `MIM=TRUE` argument.

## To install:
```
require(devtools)

install_github("glenmcgee/MVmixture")
```
