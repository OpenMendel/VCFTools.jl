module VCFTools

using CodecZlib, Distributions, GeneticVariation, ProgressMeter

export gtstats, gtstats!, hwe

# package code goes here
include("gtstats.jl")
#include("conformgt")

end # module
