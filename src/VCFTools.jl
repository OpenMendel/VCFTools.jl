module VCFTools

using CodecZlib, Distributions, GeneticVariation.VCF, HypothesisTests,
    NullableArrays, ProgressMeter

export conformgt_by_id, conformgt_by_pos,
    convert_gt, copy_gt!,
    filter_genotype, gtstats,
    nrecords, nsamples, openvcf

# package code goes here
include("gtstats.jl")
include("conformgt.jl")
include("convert.jl")

end # module
