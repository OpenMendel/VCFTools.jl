"""
    aim_select(vcffile::AbstractString)

Ranks SNPs by their ancestry information content. All people should 
be assigned ancestry fractions and be fully typed. Ranks are assigned
by a likelihood ratio heterogeneity test. Sex chromosome is ignored. 

# Note:
Method is adapted from MendelAimSelection.jl
https://github.com/OpenMendel/MendelAimSelection.jl
"""
function aim_select(
    vcffile::AbstractString, 
    sampleID_to_population::Dict{String, String}
    )
    pvalues = zeros(nrecords(vcffile))
    populations = unique_populations(sampleID_to_population)
    alleles = Vector{Int}(undef, length(populations)),
    genes = Vector{Int}(undef, length(populations))

    # loop over SNPs
    reader = VCF.Reader(openvcf(vcffile, "r"))
    for (i, record) in enumerate(reader)
        # if no "GT" field, skip this record
        if VCF.findgenokey(record, "GT") === nothing
            pvalues[i] = 1.0
            continue
        end
        pvalues[i] = aim_select(record, sampleID_to_population, populations,
            alleles, genes)
    end
    close(out); close(reader)

    return pvalues
end

function aim_select(
    record::VCF.Record,
    sampleID_to_population::Dict{String, String},
    populations::Vector{String},
    alleles::Vector{Int} = Vector{Int}(undef, length(populations)),
    genes::Vector{Int} = Vector{Int}(undef, length(populations))
    )
    # summarize current record
    n00, n01, n11, n0, n1, altfreq, reffreq, missings,
        minorallele, maf, hwepval = gtstats(record)

    # skip SNPs with maf too small
    maf ≤ 0.01 && return 1.0

    # Tally reference alleles and genes in each population
    people = length(sampleID_to_population)
    for (person, ethnic) in sampleID_to_population
        # which population this person belongs
        j = something(findfirst(x -> x == ethnic, populations))
        # get genotype: "0" (REF) => 0x30, "1" (ALT) => 0x31
        geno = record.genotype[i]
        a1 = record.data[geno[gtkey][1]] == 0x31
        a2 = record.data[geno[gtkey][3]] == 0x31
        # increment counters
        alleles[j] += a1 + a2
        genes[j] += 2
    end

    # Add the maximum loglikelihoods for the different populations
    lrt = 0.0
    for j = 1:length(populations)
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
    unique_populations(x::Dict{String, String})

Computes the unique list of populations, preserving order. `x` is a `Dict`
where sample IDs are keys and populations are values. 
"""
function unique_populations(x::Dict{String, String})
    populations = String[]
    for (key, val) in x
        val ∉ populations && push!(populations, val)
    end
    return populations
end