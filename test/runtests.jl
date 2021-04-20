using VCFTools
using Test
using VariantCallFormat
using CodecZlib
using DelimitedFiles

# packages needed only for testing
using Random
using CSV
using DataFrames
using StatsBase
using LinearAlgebra

include("gtstats_test.jl")
include("conformgt_test.jl")
include("convert_test.jl")
include("filter_test.jl")
include("aim_test.jl")
include("grm_test.jl")
