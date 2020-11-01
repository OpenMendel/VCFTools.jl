"""
    aim_select(vcffile::AbstractString, sampleID_to_population::Dict{String, String}, 
        [maf_threshold::Float64=0.01])

Ranks SNPs by their ancestry information content. All people should 
be assigned an ancestry origin and be fully typed. P-values are computed 
via a likelihood ratio heterogeneity test. Sex chromosome is ignored (for now). 

# Inputs
- `vcffile`: VCF file, ending with `.vcf` or `.vcf.gz`
- `sampleID_to_population`: A dictionary mapping every sample ID to a population
    origin. 

# Optional Inputs
- `maf_threshold`: Minimum minor allele frequency for SNPs considered (default 0.01)

# Note:
Method is adapted from MendelAimSelection.jl
https://github.com/OpenMendel/MendelAimSelection.jl
"""
function aim_select(
    vcffile::AbstractString, 
    sampleID_to_population::Dict{String, String},
    maf_threshold::Float64 = 0.01,
    )
    0.0 < maf_threshold < 0.5 || error("maf_threshold should be between 0 and 0.5.")

    # data information
    reader = VCF.Reader(openvcf(vcffile, "r"))
    sampleID = VCF.header(reader).sampleID
    # ethnics = sort!(ethnic(sampleID, sampleID_to_population)) # to match with MendelAimSelection
    ethnics = ethnic(sampleID, sampleID_to_population)
    populations = unique(ethnics)

    # preallocated vectors
    alleles = Vector{Int}(undef, length(populations))
    genes = Vector{Int}(undef, length(populations))
    records = nrecords(vcffile)
    pvalues = zeros(records)

    # loop over SNPs
    pbar = ProgressMeter.Progress(records, 5)
    for (i, record) in enumerate(reader)
        pvalues[i] = aim_select(record, ethnics, populations, alleles, genes, maf_threshold)
        next!(pbar)
    end
    close(reader)

    return pvalues
end

function aim_select(
    record::VCF.Record,
    ethnics::Vector{String},
    populations::Vector{String} = unique(ethnics),
    alleles::Vector{Int} = Vector{Int}(undef, length(populations)),
    genes::Vector{Int} = Vector{Int}(undef, length(populations)),
    maf_threshold::Float64 = 0.01
    )
    fill!(alleles, 0)
    fill!(genes, 0)

    # summarize current record
    people = length(ethnics)
    n00, n01, n11, n0, n1, altfreq, reffreq, missings,
        minorallele, maf, hwepval = gtstats(record)

    # skip SNPs with maf too small
    maf ≤ maf_threshold && return 1.0

    # Tally reference alleles and genes in each population
    @inbounds for i in 1:people
        # which population this person belongs
        j = something(findfirst(isequal(ethnics[i]), populations))
        # get genotype: "0" (REF) => 0x30, "1" (ALT) => 0x31
        gtkey = VCF.findgenokey(record, "GT")
        gtkey === nothing && return 1.0 # if no "GT" field, skip this record
        geno = record.genotype[i]
        if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
            return 1.0 # if genotype missing, skip this record
        end
        a1 = record.data[geno[gtkey][1]] == 0x31
        a2 = record.data[geno[gtkey][3]] == 0x31
        # increment counters
        alleles[j] += a1 + a2
        genes[j] += 2
    end

    # Add the maximum loglikelihoods for the different populations
    lrt = 0.0
    @inbounds for j in 1:length(populations)
        p = 0.0
        if genes[j] > 0
            p = alleles[j] / genes[j]
        end
        if 0.0 < p < 1.0
            n = genes[j]
            x = alleles[j]
            lrt += x * log(p) + (n - x)*log(1.0 - p)
        end
    end

    # Subtract the maximum loglikelihood for the entire sample.
    # Based on this, compute the likelihood ratio p-value.
    p = sum(alleles) / sum(genes)
    if 0.0 < p < 1.0
        n = sum(genes)
        x = sum(alleles)
        lrt = 2.0*(lrt - x * log(p) - (n - x)*log(1.0 - p))
        return ccdf(Chisq(1), lrt) 
    else
        return 1.0
    end
end

"""
    ethnic(sampleID::Vector{String}, sampleID_to_population::Dict{String, String})

Computes the population origin for each sample in `sampleID`, returns a dictionary

# Inputs
- `sampleID`: Vector of sample IDs
- `sampleID_to_population`: a `Dict` where sample IDs are keys and populations
    are values.
"""
function ethnic(
    sampleID::Vector{String}, 
    sampleID_to_population::Dict{String, String}
    )
    ethnics = Vector{String}(undef, length(sampleID))
    for (i, id) in enumerate(sampleID)
        @inbounds ethnics[i] = sampleID_to_population[id]
    end
    return ethnics
end
