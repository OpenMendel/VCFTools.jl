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
    save_snpinfo,
    convert_gt, copy_gt!,
    convert_ht, copy_ht!,
    convert_ds, copy_ds!, 
    copy_gt_trans!, copy_ht_trans!, copy_ds_trans!,
    # filter functions
    filter_genotype, 
    filter, filter_header, 
    filter_chr, filter_range,
    mask_gt, find_duplicate_marker,
    unsafe_convert_gt, unsafe_convert_gt2, unsafe_convert_gt3

# package code goes here
include("gtstats.jl")
include("conformgt.jl")
include("convert.jl")
include("filter.jl")
include("unsafe_convert.jl")

end # module
