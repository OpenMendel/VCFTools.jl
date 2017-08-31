module TestVCFTools

using VCFTools
using Base.Test

include("gtstats_test.jl")
include("conformgt_test.jl")
include("convert_test.jl")

end
