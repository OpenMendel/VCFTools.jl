"""
    aim_select(vcffile::AbstractString)

Ranks SNPs by their ancestry information content. All people should 
be assigned ancestry fractions and be fully typed. Ranks are assigned
by a likelihood ratio heterogeneity test. Sex chromosome is ignored. 

# Inputs
- `vcffile`: VCF file, ending with .vcf or .vcf.gz
- `sampleID_to_population`: A dictionary mapping every sample ID to a population
    origin. 

# Note:
Method is adapted from MendelAimSelection.jl
https://github.com/OpenMendel/MendelAimSelection.jl
"""
function aim_select(
    vcffile::AbstractString, 
    sampleID_to_population::Dict{String, String}
    )
    reader = VCF.Reader(openvcf(vcffile, "r"))
    sampleID = VCF.header(reader).sampleID
    ethnics = ethnic(sampleID, sampleID_to_population)
    populations = unique(ethnics)

    # preallocated vectors
    alleles = Vector{Int}(undef, length(ethnics))
    genes = Vector{Int}(undef, length(ethnics))
    records = nrecords(vcffile)
    pvalues = zeros(records)

    # loop over SNPs
    pbar = ProgressMeter.Progress(records, 5)
    for (i, record) in enumerate(reader)
        pvalues[i] = aim_select(record, ethnics, alleles, genes)
        next!(pbar)
    end
    close(out); close(reader)

    return pvalues
end

function aim_select(
    record::VCF.Record,
    ethnics::Vector{String},
    alleles::Vector{Int},
    genes::Vector{Int}
    )
    fill!(alleles, 0)
    fill!(genes, 0)

    # summarize current record
    people = length(ethnics)
    n00, n01, n11, n0, n1, altfreq, reffreq, missings,
        minorallele, maf, hwepval = gtstats(record)

    # skip SNPs with maf too small
    maf ≤ 0.01 && return 1.0

    # Tally reference alleles and genes in each population
    for i in 1:people
        # which population this person belongs
        ethnic = ethnics[i]
        j = something(findfirst(x -> x == ethnic, ethnics))
        # get genotype: "0" (REF) => 0x30, "1" (ALT) => 0x31
        geno = record.genotype[i]
        gtkey = VCF.findgenokey(record, "GT")
        gtkey === nothing && return 1.0 # if no "GT" field, skip this record
        a1 = record.data[geno[gtkey][1]] == 0x31
        a2 = record.data[geno[gtkey][3]] == 0x31
        # increment counters
        alleles[j] += a1 + a2
        genes[j] += 2
    end

    # Add the maximum loglikelihoods for the different populations
    lrt = 0.0
    for j = 1:length(ethnics)
        p = 0.0
        if genes[j] > 0
            p = alleles[j] / genes[j]
        end
        if 0.0 < p < 1.0
            n = genes[j]
            x = alleles[j]
            lrt = lrt + x * log(p) + (n - x)*log(1.0 - p)
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

Computes the the population origin for each sample in `sampleID`. 
`sampleID_to_population` is a `Dict` where sample IDs are keys and populations
are values. 
"""
function ethnic(
    sampleID::Vector{String}, 
    sampleID_to_population::Dict{String, String}
    )
    populations = Vector{String}(undef, length(sampleID))
    for (i, id) in enumerate(sampleID)
        populations[i] = sampleID_to_population[id]
    end
    return populations
end
