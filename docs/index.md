# VCFTools.jl

*Julia utilities for handling VCF files*

`VCFTools.jl` implements some Julia utilities for handling [VCF](https://github.com/samtools/hts-specs) files. This package heavily uses the VCF parser developed in the [`BioJulia/GeneticVariation.jl`](https://github.com/BioJulia/GeneticVariation.jl) package.

## Package Features

- Calculate genotype statistics (minor allele frequencies, minor allele, count of missing genotypes, etc) from a VCF file. 
- Extract specific data fields from a VCF file. 
- Filtering/subsetting VCF files. 
- Match markers in two VCF files according to ID. 

## Installation

Use the Julia package manager to install `VCFTools.jl`.
```julia
]add https://github.com/biona001/VCFTools.jl
```
This package supports Julia `1.0`+.

## Manual Outline

```@contents
Pages = [
    "summaryinfo.md",
    "filter.md",
    "convert.md",
    "conformgt.md",
    "api.md"
]
Depth = 2
```
