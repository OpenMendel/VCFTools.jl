# VCFTools.jl

*Julia utilities for handling VCF files*

`VCFTools.jl` implements some Julia utilities for handling [VCF](https://github.com/samtools/hts-specs) files. This package heavily uses the VCF parser developed in the [`BioJulia/GeneticVariation.jl`](https://github.com/BioJulia/GeneticVariation.jl) package.

## Package Features

- Calculate genotype statistics (minor allele frequencies, minor allele, count of missing genotypes, etc) from a VCF file. 
- Extract specific data fields from a VCF file. 
- Filtering/subsetting VCF files. 
- Match markers in two VCF files according to ID. 
- Calculation of Genetic Relationship Matrix.
- Rank SNPs by their ancestry informativeness.

## Installation

Use the Julia package manager to install `VCFTools.jl`. Press `]` and type:
```julia
(v1.5) pkg> add https://github.com/OpenMendel/VCFTools.jl
```
This package supports Julia `v1.5`+.

## Manual Outline

```@contents
Pages = [
    "man/summaryinfo.md",
    "man/filter.md",
    "man/convert.md",
    "man/conformgt.md",
    "man/grm.md"
    "man/aim.md"
    "man/api.md"
]
Depth = 2
```
