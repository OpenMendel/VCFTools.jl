# VCFTools

| **Documentation** | **Build Status** | **Code Coverage**  | **Citation**  |  
|-------------------|------------------|--------------------|--------------------|  
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/VCFTools.jl/dev/) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/VCFTools.jl/stable) | [![build Actions Status](https://github.com/OpenMendel/VCFTools.jl/workflows/CI/badge.svg)](https://github.com/OpenMendel/VCFTools.jl/actions) [![CI (Julia nightly)](https://github.com/openmendel/VCFTools.jl/workflows/JuliaNightly/badge.svg)](https://github.com/OpenMendel/VCFTools.jl/actions/workflows/JuliaNightly.yml) | [![codecov](https://codecov.io/gh/OpenMendel/VCFTools.jl/branch/master/graph/badge.svg?token=QtTQogesUk)](https://codecov.io/gh/OpenMendel/VCFTools.jl) | [![DOI](https://zenodo.org/badge/100287089.svg)](https://zenodo.org/badge/latestdoi/100287089) |

VCFTools.jl provide Julia utilities for handling VCF files.

## Installation


This package is registered in the default Julia package registry, and can be installed through standard package installation procedure: e.g., running the following code in Julia REPL.
```julia
using Pkg
pkg"add VCFTools"
```

This package supports Julia v1.5+.

## Citation

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. OPENMENDEL: a cooperative programming project for statistical genetics. Hum Genet. 2020 Jan;139(1):61-71. doi: 10.1007/s00439-019-02001-z. Epub 2019 Mar 26. PMID: 30915546; PMCID: [PMC6763373](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6763373/).*

## Acknowledgments

This project has been supported by the National Institutes of Health under awards R01GM053275, R01HG006139, R25GM103774, and 1R25HG011845.
