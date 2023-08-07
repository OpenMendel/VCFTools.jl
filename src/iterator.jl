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
        reader = VCF.Reader(openvcf(itr.vcffile, "r"))
        vector = []
        chr = ""
        pos = 0
        ids = ""
        ref = ""
        alt = ""
        qual = 0.0

        count = 0

        for record in reader
            chr = VCF.chrom(record)
            pos = VCF.pos(record)
            ids = VCF.id(record) # VCF.id(record)[1]
            ref = VCF.ref(record)
            alt = VCF.alt(record) # VCF.alt(record)[1]
            qual = VCF.qual(record)
            push!(vector, (chr, pos, ids, ref, alt, qual))

            count += 1

            if count == state 
                break
            end
        end

        close(reader)

        if isempty(vector)
            return nothing
        else
            if state == count 
                println("Chromosome: $chr, Position: $pos, IDs: $ids, REF: $ref, ALT: $alt, QUAL: $qual")
                return (vector, state + 1)
            end
        end
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