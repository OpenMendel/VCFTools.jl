"""
    grm(vcffile::AbstractString, method = :GRM, minmaf=0.01, colinds=nothing)

Computes the genetic relationship matrix of a VCF file using specified method.
Missing data is handled according to:
`https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6220858/`
"""
function grm(
    vcffile::AbstractString;
    method::Symbol = :Robust,
    minmaf::Real = 0.01,
    cinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    t::Type{T} = Float64,
    scale_missing::Bool = false
    ) where T <: Real

    # load genotype into memory. Each row is a sample
    x = convert_gt(t, vcffile, model=:additive, msg="importing genotypes")

    # compute summary statistics
    n, p = size(x)
    mis = zeros(n) # proportion of missingness SNPs in each sample
    maf = zeros(p) # minor allele frequency
    wts = zeros(p) # inverse standard deviation
    cts = zeros(p) # mean
    @inbounds for i in 1:p
        n00 = n01 = n11 = nmissing = 0
        for j in 1:n
            if ismissing(x[j, i])
                nmissing += 1
                mis[j] += 1
            elseif x[j, i] == 0
                n00 += 1
            elseif x[j, i] == 1
                n01 += 1
            elseif x[j, i] == 2
                n11 += 1
            else
                error("genotype is not 0, 1, 2, or missing!")
            end
        end
        n0 = 2n00 + n01
        n1 = 2n11 + n01
        altfreq = n1 / (n0 + n1)
        reffreq = n0 / (n0 + n1)
        maf[i] = altfreq < reffreq ? altfreq : reffreq
        cts[i] = 2altfreq
        wts[i] = altfreq == 0 ? 1.0 : 1.0 / √(2altfreq * (1 - altfreq))
    end

    mis ./= p
    colinds = something(cinds, maf .≥ minmaf)

    return @views grm!(x[:, colinds], method, maf[colinds], wts[colinds], 
        cts[colinds], mis, scale_missing)
end

"""
    grm!(G, method, maf, wts, cts)

Computes GRM of `G` where `Gᵢ ∈ {0, 1, 2, missing}`. `G` is modified. 

Missing data is handled according to:
`https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6220858/`
"""
function grm!(
    G::AbstractMatrix,
    method::Symbol,
    maf::AbstractVector,
    wts::AbstractVector,
    cts::AbstractVector,
    mis::AbstractVector,
    scale_missing::Bool = false
    )
    m, n = size(G)
    if method == :GRM
        _impute_center_scale!(G, wts, cts, impute=true, center=true, scale=true)
        Φ = G * transpose(G)
        Φ ./= 2n
    elseif method == :MoM
        _impute_center_scale!(G, wts, cts, impute=true)
        G .-= 1
        Φ = G * transpose(G)
        c = sum(x -> abs2(x) + abs2(1 - x), maf)
        shft, scal = n / 2 - c, 1 / (n - c)
        @inbounds @simd for i in eachindex(Φ)
            Φ[i] = (Φ[i] / 2 + shft) * scal
        end
    elseif method == :Robust
        _impute_center_scale!(G, wts, cts, impute=true, center=true)
        scal = sum(x -> 4x * (1 - x), maf)
        Φ = G * transpose(G)
        Φ ./= scal
    else
        throw(ArgumentError("method should be :GRM, :MoM, or :Robust"))
    end

    # scale missing data
    if scale_missing
        @inbounds for i in 1:m
            @simd for j in 1:m
                Φ[j, i] /= ((1 - mis[i]) * (1 - mis[j]))
            end
        end
    end

    return Φ
end

function _impute_center_scale!(
    x::AbstractMatrix,
    wts::AbstractVector,
    cts::AbstractVector;
    impute::Bool=false,
    center::Bool=false,
    scale::Bool=false
    )
    m, n = size(x)
    @inbounds for i in 1:n
        ct = cts[i]
        wt = wts[i]
        for j in 1:m
            impute && ismissing(x[j, i]) && (x[j, i] = ct)
            center && (x[j, i] -= ct)
            scale && wt > 0 && (x[j, i] *= wt)
        end
    end
    return nothing
end
