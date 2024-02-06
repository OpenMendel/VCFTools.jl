abstract type VariantIterator end
abstract type Variant end
abstract type GeneticData end

mutable struct VCFIterator <: VariantIterator
    vcffile::AbstractString
    vcf::VCF.Reader
    sample_num::Int
end

function VCFIterator(vcffile::AbstractString; sample_num::Int=nsamples(VCF.Reader(openvcf(vcffile, "r"))))
    vcf = VCF.Reader(openvcf(vcffile, "r")) 
    return VCFIterator(vcffile, vcf, sample_num) #then returns an iostream 
end

mutable struct VCFData <: GeneticData # feed VCFData VCFIterator function output output of openvcf function 
    file_name::AbstractString
end
# placeholder because other genotype files need that 

mutable struct VCFRow <: Variant
    CHROM::String
    POS::Int64
    ID::Vector{String}
    REF::String
    ALT::Vector{String}
    QUAL::Float64
    GENOTYPE::Vector{Union{Float64, Missing}}
    DOSAGES::Vector{Union{Float64, Missing}}
    INDEX::Int64
end

# for dosages there is a key DS 
# copy_gt! and copy_ds!
# for genotypes taking values from GT

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
        geno_out = zeros(Union{Missing, Float64}, nsamples(itr.vcffile))
        ds_out = zeros(Union{Missing, Float64}, nsamples(itr.vcffile))
        geno = copy_gt!(geno_out, reader)
        ds = copy_ds!(ds_out,reader)
        # dos = copy_ds

        count = 0

        for record in reader
            vector = []
            count += 1
            chr = VCF.chrom(record)
            pos = VCF.pos(record)
            ids = VCF.id(record) # VCF.id(record)[1]
            ref = VCF.ref(record)
            alt = VCF.alt(record) # VCF.alt(record)[1]
            qual = VCF.qual(record)
            # out = zeros(Union{Missing, Float64}, nsamples(itr.vcffile))
            # geno = copy_gt!(out, reader)
            # return a tuple not an array 
            index = count 
            vcf_row = VCFRow(chr, pos, ids, ref, alt, qual, geno, ds, index)

            if count == state 
                close(reader)
                return (vcf_row, state + 1)
            end
        end

        close(reader)
        return nothing 

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

function maf(data::VCFData)

    records, samples, lines, missing_by_sample, missings_by_record, maf_by_record, minor_allele_by_record, hwe_by_record = gtstats(data.file_name)
    maf = maf_by_record[row.INDEX]
    return maf

    # gt_row = row.GENOTYPE
    
    # ref_count = 0
    # alt_count = 0


    # s = 0.0 # total of the alternate allele count 
    # count = 0 # count of non missing samples 
    # for i in 1:sample_num
    #      if gt_row[i] !== missing
    #         s += gt_row[i]
    #         count += 1
    #     end
    # end
    # alt_allele_frequency = s / count
   
    # if alt_allele_frequency > 0.5
    #     return 1 - alt_allele_frequency
    # else
    #     return alt_allele_frequency
    # end

    # ./.
    # account for missing genotypes 

    # each iteration of the iterator gives you VCFRow 
    # for each iteration take the average over the vector
    # divide average over the samples by 2

end

#copyto! function
#copygt! copydt
#reads in genotypes 0 1 2 average divided by two allele frequency of alternate allele 
#dosages for each snp 0-2

function hwepval(data::VCFData, row::VCFRow)
    
    records, samples, lines, missing_by_sample, missings_by_record, maf_by_record, minor_allele_by_record, hwe_by_record = gtstats(data.file_name)
    p_value = hwe_by_record[row.INDEX]
    return p_value 

    # ref_count = 0
    # alt_count = 0
    # gt_row = row.GENOTYPE

    # for genotype in out
    #     alleles = split(genotype, '/')
    #     alleles = string.(alleles)
    #     if "0" in alleles
    #         ref_count += count(x -> x == "0", alleles)
    #     elseif "1" in alleles
    #         alt_count += count(x -> x == "1", alleles)
    #     else
    #         continue  # Skip missing genotypes
    #     end
    # end
    
    # total_samples = ref_count + alt_count

    # p_obs = ref_count / (2 * total_samples)
    # q_obs = alt_count / (2 * total_samples)

    # p_exp = p_obs^2
    # q_exp = 2 * p_obs * q_obs
    # q_exp2 = q_obs^2

    # chi_squared = sum((genotype_counts[geno] - total_samples * expected_frequency)^2 / (total_samples * expected_frequency) for (geno, expected_frequency) in [("0/0", p_exp), ("0/1", q_exp), ("1/1", r_exp)])

    # # Using Distributions.jl 

    # p_value = 1.0 - cdf(Chisq(2), chi_squared)

end

function alt_dosages!(arr::AbstractArray{T}, data::VCFData, row::VCFRow; use_genotype::Bool=false) where T <: Real
    if use_genotype
        genotypes = row.GENOTYPE
        copyto!(arr, genotypes)
    else
        dosages = row.DOSAGES
        for i in 1:length(dosages)
            arr[i] = dosages[i]
        end
    end
   
    return arr
end


function alt_genotypes!(arr::AbstractArray{T}, data::VCFData, row::VCFRow) where T <: Real
    # Extract genotype information from the VCFRow
    genotypes = row.GENOTYPE

    for i in 1:length(genotypes)
        if genotypes[i] != 0.0  # Not homozygous for the reference allele
            # Add the genotype to the array
            arr[i] = genotypes[i]
        else
            # Set to missing if not an alternate genotype
            arr[i] = missing
        end
    end

    return arr
end


# filter.jl line 188

#inside snparrays hwe