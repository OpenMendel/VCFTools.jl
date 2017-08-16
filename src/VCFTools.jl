module VCFTools

using Distributions, GeneticVariation, CodecZlib, ProgressMeter

export gtstats, hwe

# package code goes here
include("gtstats.jl")
#include("conformgt")

end # module
