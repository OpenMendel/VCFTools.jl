abstract type VariantIterator end
abstract type Variant end

mutable struct VCFIterator <: VariantIterator
    vcffile::AbstractString
    vcf::VCF.Reader
end

function VCFIterator(vcffile::AbstractString)
    vcf = VCF.Reader(openvcf(vcffile, "r"))
    return VCFIterator(vcffile, vcf)
end

mutable struct VCFIndex <: Variant
    index::Int
end

mutable struct VCFRow <: Variant
    CHROM::String
    POS::Int64
    ID::String
    REF::String
    ALT::String
    QUAL::Float64
end

@inline function Base.eltype(::Type{<:VariantIterator})
    VCFRow
end

function Base.iterate(itr::VCFIterator, state=1)
    rows = nrecords(itr.vcffile)
    if state <= 0 || state > rows
        return nothing
    else
        count = 1
        for record in itr.vcf
            if count == state
                result = (view(record),state+1)
                break
            else
                count = count + 1 
            end
        end
        state = state + 1
        return result
    end
end

@inline function Base.length(itr::VCFIterator)
    return nrecords(itr.vcffile)
end

@inline function vcfiterator(vcffile::AbstractString)
    VCFIterator(vcffile)
end

function chrom(vcffile::AbstractString, state=1)::String
    return itr.vcf[state].CHROM
end

function pos(vcffile::AbstractString, state=1)::Int
    return itr.vcf[state].POS
end

function rsid(vcffile::AbstractString, state=1)::String
    return itr.vcf[state].ID
end

function alleles(vcffile::AbstractString, state=1)::Vector{String}
    return itr.vcf[state].ALT
end

function alt_allele(vcffile::AbstractString, state=1)::String
    return itr.vcf[state].ALT[1]
end

function ref_allele(vcffile::AbstractString, state=1)::String
    return itr.vcf[state].REF
end