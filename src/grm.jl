"""
    grm(vcffile::AbstractString, maf::AbstractVector, [method = :GRM])

Computes the genetic relationship matrix of a VCF file using specified method.
Missing data is handled according to:
`https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6220858/`
"""
function grm(
    vcffile::AbstractString,
    maf::AbstractVector,
    method::Symbol = :GRM,
    t::Type{T} = Float64
    ) where T <: Real

    if method == :GRM
        G = convert_gt(t, vcffile, model=ADDITIVE_MODEL, impute=true, 
            center=true, scale=true)
        Φ = G * transpose(G)
        Φ ./= 2n
    elseif method == :MoM
        G = convert_gt(t, vcffile, model=ADDITIVE_MODEL, impute=true)
        G .-= 1
        Φ = G * transpose(G)
        c = sum(x -> abs2(x) + abs2(1 - x), maf)
        shft, scal = n / 2 - c, 1 / (n - c)
        @inbounds @simd for i in eachindex(Φ)
            Φ[i] = (Φ[i] / 2 + shft) * scal
        end
    elseif method == :Robust
        G = convert_gt(t, vcffile, model=ADDITIVE_MODEL, center=true, impute=true)
        scal = sum(x -> 4x * (1 - x), maf)
        Φ = G * transpose(G)
        Φ ./= scal
    else
        throw(ArgumentError("method should be :GRM, :MoM, or :Robust"))
    end
    Φ
end