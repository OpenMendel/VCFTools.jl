"""
    gtstats(vcffile, [out=devnull])

Calculate genotype statistics for each marker with GT field in a VCF file.

# Input
- `vcffile`: VCF file, ending with .vcf or .vcf.gz
- `out`: output file name or IOStream. Default is `out=devnull` (no output).
One line with 15 tab-delimited fiels is written per marker to `out`:
    - 1-8)  VCF fixed fields (CHROM, POS, ID, REF, ALT, QUAL, FILT, INFO)
    -   9)  Missing genotype count
    -  10)  Missing genotype frequency
    -  11)  ALT allele count
    -  12)  ALT allele frequency
    -  13)  Minor allele count             (REF allele vs ALT alleles)
    -  14)  Minor allele frequency         (REF allele vs ALT alleles)
    -  15)  HWE P-value                    (REF allele vs ALT alleles)

# Output
- `records`: number of records in the input VCF file
- `samples`: number of individuals in the input VCF file
- `lines`  : number of lines written to `out`; equivalently number of markers
    with GT field
- `missings_by_sample`: number of missing genotypes in each sample
- `missings_by_record`: number of missing genotypes in each record (marker)
- `maf_by_record`: minor allele frequency in each record (marker)
- `minorallele_by_record`: a Boolean vector indicating the minor allele in each
    record (marker). `minorallele_by_record[i]=true` means the minor allele is
    the REF allele for marker `i`; `minorallele_by_record[i]=false` means the
    minor allele is the ALT allele for marker `i`
"""
function gtstats(
    vcffile::AbstractString, 
    out::IO=devnull, 
    pval_method::Symbol = :Pearson
    )
    # open VCF file
    reader = VCF.Reader(openvcf(vcffile, "r"))
    # set up progress bar
    records = nrecords(vcffile)
    out == stdout || (pbar = ProgressMeter.Progress(records, 1))
    # allocate ouput arrays
    samples = nsamples(reader)
    missings_by_sample = zeros(Int, samples)
    missings_by_record = Int[]
    maf_by_record = Float64[]
    minorallele_by_record = Bool[]
    sizehint!(maf_by_record, records)
    sizehint!(missings_by_record, records)
    sizehint!(minorallele_by_record, records)
    # loop over records
    records = lines = 0
    for record in reader
        records += 1
        # if no "GT" field, skip this record
        VCF.findgenokey(record, "GT") === nothing && continue
        # calcuate summary statistics
        lines += 1
        n00, n01, n11, n0, n1, altfreq, reffreq, missings, minorallele, maf, 
            hwepval = gtstats(record, missings_by_sample, pval_method)
        missfreq = missings / (n0 + n1)
        altfreq  = n1 / (n0 + n1)
        minoralleles = minorallele ? n0 : n1
        push!(missings_by_record, missings)
        push!(maf_by_record, maf)
        push!(minorallele_by_record, minorallele)
        # output
        nbytes = record.format[1][1] - 2
        unsafe_write(out, pointer(record.data), nbytes)
        print(out, '\t', missings, '\t', missfreq, '\t', n0, '\t',
        altfreq, '\t', minoralleles, '\t', maf, '\t', hwepval, '\n')
        # update progress bar
        out == stdout || ProgressMeter.update!(pbar, records)
    end
    close(out); close(reader)
    return records, samples, lines, missings_by_sample, missings_by_record,
        maf_by_record, minorallele_by_record
end

function gtstats(vcffile::AbstractString, out::AbstractString)
    ofile = endswith(out, ".gz") ? GzipCompressorStream(open(out, "w")) : open(out, "w")
    gtstats(vcffile, ofile)
end

"""
    gtstats(record, [missings_by_sample=nothing])

Calculate genotype statistics for a VCF record with GT field.

# Input
- `record`: a VCF record
- `missings_by_sample`: accumulator of missings by sample, `missings_by_sample[i]` is incremented by 1 if `i`-th individual has missing genotype in this record

# Output
- `n00`: number of homozygote REF/REF or REF|REF
- `n01`: number of heterozygote REF/ALT or REF|ALT
- `n11`: number of homozygote ALT/ALT or ALT|ALT
- `n0`: number of REF alleles
- `n1`: number of ALT alleles
- `altfreq`: proportion of ALT alleles
- `reffreq`: proportion of REF alleles
- `missings`: number of missing genotypes
- `minorallele`: minor allele: `false` (ALT is minor allele) or `true` (REF is minor allele)
- `maf`: minor allele frequency
- `hwepval`: Hardy-Weinberg p-value
"""
function gtstats(
    record::VCF.Record,
    missings_by_sample::Union{Vector,Nothing}=nothing,
    pval_method::Symbol=:Pearson
    )
    # n11: number of homozygote ALT/ALT or ALT|ALT
    # n00: number of homozygote REF/REF or REF|REF
    # n01: number of heterozygote REF/ALT or REF|ALT
    missings = n00 = n10 = n11 = 0
    gtkey = VCF.findgenokey(record, "GT")
    for i in 1:lastindex(record.genotype)
        geno = record.genotype[i]
        # dropped field or "." => 0x2e
        if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
            missings += 1
            missings_by_sample === nothing || (missings_by_sample[i] += 1)
        else
            # "0" => 0x30, "1" => 0x31
            if record.data[geno[gtkey][1]] == 0x30
                n00 += record.data[geno[gtkey][3]] == 0x30 ? 1 : 0
            elseif record.data[geno[gtkey][1]] == 0x31
                n11 += record.data[geno[gtkey][3]] == 0x31 ? 1 : 0
            end
        end
    end
    n01         = length(record.genotype) - missings - n00 - n11
    n0          = n01 + 2n00
    n1          = n01 + 2n11
    altfreq     = n1 / (n0 + n1)
    reffreq     = n0 / (n0 + n1)
    minorallele = n0 < n1 # true if REF is minor 
    maf         = minorallele ? reffreq : altfreq
    hwepval     = hwe(n00, n01, n11, pval_method)
    return n00, n01, n11, n0, n1, altfreq, reffreq, missings,
        minorallele, maf, hwepval
end

"""
    hwe(n00, n01, n11, [method=:Pearson])

Hardy-Weinberg equilibrium test using Pearson Chi squared test or Fisher's
exact test. 
"""
function hwe(n00::Integer, n01::Integer, n11::Integer, method::Symbol=:Pearson)
    n = n00 + n01 + n11
    n == 0 && return 1.0
    p0 = (n01 + 2n00) / 2n
    (p0 ≤ 0.0 || p0 ≥ 1.0) && return 1.0
    p1 = 1 - p0
    if method == :Pearson
        # Pearson's Chi-squared test
        e00 = n * p0 * p0
        e01 = 2n * p0 * p1
        e11 = n * p1 * p1
        ts = (n00 - e00)^2 / e00 + (n01 - e01)^2 / e01 + (n11 - e11)^2 / e11
        pval = ccdf(Chi(1), ts)
    elseif method == :Fisher
        # Fisher exact test
        n0 = n01 + 2n00
        pval, _, _ = hwe_fisher(n01, n, n0)
    end
    return pval
end

"""
    fisher_exact(n01::Integer, N::Integer, n0::Integer)

Computes probability of observing n01 heterozygote genotypes. 

# Arguments:
- `n01` = heterozygous counts
- `N` = total samples with nonmissing genotypes
- `n0` = number of 0 alleles = `2*n00 + n01`
"""
function fisher_exact(n01::Integer, N::Integer, n0::Integer)
    n1  = 2N - n0
    n00 = (n0 - n01) >> 1 # fast division by 2
    n11 = (n1 - n01) >> 1
    logp = n01*log(2) + logfactorial(N) - logfactorial(n00) - 
        logfactorial(n01) - logfactorial(n11) + logfactorial(n0) + 
        logfactorial(n1) - logfactorial(2N)
    return exp(logp)
end

"""
    hwe_fisher(o01::Integer, N::Integer, n0::Integer)

Computes p-value for deviation from Hardy weinburg equilibrium using Fisher's
exact test. Also computes tail probabilities P(n01 ≤ o01) and P(n01 ≥ o01) 
where n01 is the number of heterozygotes and o01 is number of observed 
heterozygotes. 

# Arguments:
- `o01` = observed heterozygous counts
- `N` = total samples with non-missing genotypes
- `n0` = number of 0 alleles = 2*n00 + n01
"""
function hwe_fisher(o01::Integer, N::Integer, n0::Integer)
    n1 = 2N - n0

    # P(number of heterozygous ≤ n01)
    ts_low = 0.0
    for n01 in isodd(n0):2:o01
        ts_low += fisher_exact(n01, N, n0)
    end

    # P(number of heterozygous ≥ n01)
    ts_high = 0.0
    for n01 in o01:2:min(n1, n0)
        ts_high += fisher_exact(n01, N, n0)
    end

    # P_hwe according to https://pubmed.ncbi.nlm.nih.gov/15789306/
    hwepval = 0.0
    for n01 in isodd(n0):2:min(n1, n0)
        nAB     = fisher_exact(o01, N, n0)
        nABstar = fisher_exact(n01, N, n0)
        if nAB ≥ nABstar
            hwepval += nABstar
        end
    end

    return hwepval, ts_low, ts_high
end

"""
    countvcflines(vcffile)

Count the number of lines in a VCF (`.vcf` or `.vcf.gz`) file.
"""
function countvcflines(vcffile::AbstractString)
    countlines(openvcf(vcffile, "r"))
end

"""
    openvcf(vcffile, [mode = "r"])

Open VCF file (`.vcf` or `.vcf.gz`) and return an IO stream.
"""
function openvcf(vcffile::AbstractString, mode::AbstractString="r")
    if endswith(vcffile, ".vcf")
        return open(vcffile, mode)
    elseif endswith(vcffile, ".vcf.gz") && mode == "r"
        return GzipDecompressorStream(open(vcffile, mode))
    elseif endswith(vcffile, ".vcf.gz") && mode ∈ ["w", "a"]
        return GzipCompressorStream(open(vcffile, mode))
    elseif endswith(vcffile, ".vcf.gz") && mode ∉ ["r", "w", "a"]
        throw(ArgumentError("mode can only be r, w, or a for vcf.gz file"))
    else
        throw(ArgumentError("vcffile name should end with vcf or vcf.gz"))
    end
end

"""
    nrecords(vcffile)

Number of records (markers) in a VCF file. Each record is a row. 
"""
function nrecords(vcffile::AbstractString)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    records = countvcflines(vcffile) - length(VCF.header(reader)) - 1
    close(reader)
    records
end

function nsamples(reader::VCF.Reader)
    length(VCF.header(reader).sampleID)
end

"""
    nsamples(vcffile)

Number of samples (individuals) in a VCF file.
"""
function nsamples(vcffile::AbstractString)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    samples = length(VCF.header(reader).sampleID)
    close(reader)
    samples
end
