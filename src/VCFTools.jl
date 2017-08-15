module VCFTools

using Distributions, GeneticVariation, GZip, ProgressMeter

export gtstats, conformgt


# package code goes here
include("gtstats.jl")
#include("conformgt")

end # module
