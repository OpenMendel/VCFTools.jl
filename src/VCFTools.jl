module VCFTools

using CodecZlib
using Distributions
using HypothesisTests
using ProgressMeter
using DelimitedFiles
using Dates
using SpecialFunctions
using VariantCallFormat
import VariantCallFormat.findgenokey

export conformgt_by_id, conformgt_by_pos,
    gtstats, geno_ismissing,
    nrecords, nsamples, openvcf,
    sampleID,
    grm,
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
    # aim selection
    aim_select

# package code goes here
include("gtstats.jl")
include("conformgt.jl")
include("convert.jl")
include("filter.jl")
include("aim_select.jl")
include("grm.jl")

# test data directory
datadir(parts...) = joinpath(@__DIR__, "..", "test", parts...)    

end # module
