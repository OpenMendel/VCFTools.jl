module VCFTools

using CodecZlib
using Distributions
using HypothesisTests
using ProgressMeter
using DelimitedFiles
using Dates
import GeneticVariation.VCF

export conformgt_by_id, conformgt_by_pos,
    gtstats, geno_ismissing,
    nrecords, nsamples, openvcf,
    # convert functions
    convert_gt, copy_gt!,
    convert_ht, copy_ht!,
    convert_ds, copy_ds!, 
    copy_gt_trans!, copy_ht_trans!, copy_ds_trans!,
    # filter functions
    filter_genotype, 
    filter, filter_header, 
    mask_gt, find_duplicate_marker

# package code goes here
include("gtstats.jl")
include("conformgt.jl")
include("convert.jl")
include("filter.jl")

end # module
