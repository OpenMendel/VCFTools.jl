# VCFTools.jl

*Julia utilities for handling VCF fiels*

`VCFTools.jl` implements some Julia utilities for handling [VCF](https://github.com/samtools/hts-specs) files. This package heavily uses the VCF parser developed in the [`BioJulia/GeneticVariation.jl`](https://github.com/BioJulia/GeneticVariation.jl) package.

## Package Features

- Calculate genotype statistics (minor allele frequencies, minor allele, count of missing genotypes, etc) from a VCF file.  
- Extract specific data fields from a VCF file.   
- Match GT records in two VCF files according to ID.  

## Installation

Use the Julia package manager to install VCFTools.jl.
```julia
Pkg.clone("https://github.com/OpenMendel/VCFTools.jl.git")
```
This package supports Julia `0.6`.

## Manual Outline

```@contents
Pages = [
    "summaryinfo.md",
    "conformgt.md"
]
Depth = 2
```
