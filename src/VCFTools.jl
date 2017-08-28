module VCFTools

using CodecZlib, Distributions, GeneticVariation.VCF, HypothesisTests, ProgressMeter

export conformgt_by_id, conformgt_by_pos, filter_genotype, gtstats, nrecords, nsamples, openvcf

# package code goes here
include("gtstats.jl")
include("conformgt.jl")

end # module
