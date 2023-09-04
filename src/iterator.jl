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

# return VCFRow item from the iterate function 

@inline function Base.eltype(::Type{<:VariantIterator}) # vector of string 
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
            vector = []
            chr = VCF.chrom(record)
            pos = VCF.pos(record)
            ids = VCF.id(record) # VCF.id(record)[1]
            ref = VCF.ref(record)
            alt = VCF.alt(record) # VCF.alt(record)[1]
            qual = VCF.qual(record)
            # return a tuple not an array 
            vcf_row = VCFRow(chr, pos, ids, ref, alt, qual)
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
                return (vcf_row, state + 1)
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
    iterator = vcfiterator(vcffile)
    result, state = iterate(iterator, state)
    return result[1]
end

function pos(vcffile::AbstractString, state=1)::Int
    iterator = vcfiterator(vcffile)
    result, state = iterate(iterator, state)
    return result[2]
end

function rsid(vcffile::AbstractString, state=1)::String
    iterator = vcfiterator(vcffile)
    result, state = iterate(iterator, state)
    return result[3]
end

function alleles(vcffile::AbstractString, state=1)::Vector{String}
    iterator = vcfiterator(vcffile)
    result, state = iterate(iterator, state)
    return result[4:5]
end

function alt_allele(vcffile::AbstractString, state=1)::String
    iterator = vcfiterator(vcffile)
    result, state = iterate(iterator, state)
    return result[5]
end

function ref_allele(vcffile::AbstractString, state=1)::String
    iterator = vcfiterator(vcffile)
    result, state = iterate(iterator, state)
    return result[4]
end