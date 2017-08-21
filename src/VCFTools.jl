module VCFTools

using CodecZlib, Distributions, GeneticVariation, ProgressMeter

export gtstats, hwe

# package code goes here
include("gtstats.jl")
#include("conformgt")

end # module
