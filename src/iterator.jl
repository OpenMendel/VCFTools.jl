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

mutable struct VCFData <: GeneticData
    file_name::AbstractString
    io::IOStream 
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

function chrom(data::VCFData, row::VCFRow)::String
    return row.CHROM
end

function pos(data::VCFData, row::VCFRow)::Int
    return row.POS
end

function rsid(data::VCFData, row::VCFRow)::String
    return row.ID
end

function alleles(data::VCFData, row::VCFRow)::Vector{String}
    return (row.REF, row.ALT)
end

function alt_allele(data::VCFData, row::VCFRow)::String
    return row.ALT
end

function ref_allele(data::VCFData, row::VCFRow)::String
    return row.REF
end

# use GeneticData as argument for genetic variant base functions 
# GeneticData, VCFRow

function maf(data::VCFData, row::VCFRow)
    copy_gt!(out, row)
    ref_count = 0
    alt_count = 0

    for genotype in out
        alleles = split(genotype, '/')
        ref_count += count(x -> x == "0", alleles)
        alt_count += count(x -> x != "0", alleles)
    end

    # ./.
    # account for missing genotypes 
    
    total_alleles = ref_count + alt_count 
    
    if total_alleles == 0 
        return 0.0
    end 

    result = min(ref_count, total_alleles - ref_count) / total_alleles

    if result > 0.5
        return 1 - result
    else
        return result
    end

    # each iteration of the iterator gives you VCFRow 
    # for each iteration take the average over the vector
    # divide average over the samples by 2

end

#copyto! function
#copygt! copydt
#reads in genotypes 0 1 2 average divided by two allele frequency of alternate allele 
#dosages for each snp 0-2

function hwepval(data::VCFData, row::VCFRow)
    copy_gt!(out, row)
    ref_count = 0
    alt_count = 0

    for genotype in out
        alleles = split(genotype, '/')
        if "0" in alleles
            ref_count += count(x -> x == "0", alleles)
        elseif "1" in alleles
            alt_count += count(x -> x == "1", alleles)
        else
            continue  # Skip missing genotypes
        end
    end
    
    total_samples = ref_count + alt_count

    p_obs = ref_count / (2 * total_samples)
    q_obs = alt_count / (2 * total_samples)

    p_exp = p_obs^2
    q_exp = 2 * p_obs * q_obs
    q_exp2 = q_obs^2

    chi_squared = sum((genotype_counts[geno] - total_samples * expected_frequency)^2 / (total_samples * expected_frequency) for (geno, expected_frequency) in [("0/0", p_exp), ("0/1", q_exp), ("1/1", r_exp)])

    p_value = 1.0 - cdf(Chisq(2), chi_squared)

    return p_value
end

#inside snparrays hwe