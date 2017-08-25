module VCFTools

using CodecZlib, Distributions, GeneticVariation.VCF, HypothesisTests, ProgressMeter

export conformgt_by_id, filter_genotype, gtstats, hwe, flip_gt_allele, match_gt_allele

# package code goes here
include("gtstats.jl")
include("conformgt.jl")

end # module
