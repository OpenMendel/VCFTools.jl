""" 
Convert a two-bit genotype (of ALT alleles) to a real number 
of type `t` according to specified SNP model. 
""" 
function convert_gt(
    t::Type{T},
    a::NTuple{2, Bool},
    model::Symbol = :additive
    ) where T <: Real
    if model == :additive
        return convert(T, a[1] + a[2])
    elseif model == :dominant
        return convert(T, a[1] | a[2])
    elseif model == :recessive
        return convert(T, a[1] & a[2])
    else
        throw(ArgumentError("un-recognized SNP model: $model"))
    end
end

"""
    copy_gt!(A, reader; ...)

Fill the rows of matrix `A` by the GT data from VCF records in `reader`.
Records with missing genotypes are converted to `missing`.

# Input
- `A`: a matrix or vector such that `eltype(A) <: Union{Missing, Real}`. Each
    row is a person's genotype. 
- `reader`: a VCF reader

# Optional argument
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`
- `msg`: A message that will be printed to indicate progress. Defaults to
    not printing. 
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 

# Output
- `A`: `ismissing(A[i, j]) == true` indicates missing genotype. If `impute=true`,
    `ismissing(A[i, j]) == false` for all entries.
"""
function copy_gt!(
    A::Union{AbstractMatrix{Union{Missing, T}}, AbstractVector{Union{Missing, T}}},
    reader::VCF.Reader;
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false,
    msg::String = "", 
    sampleID::Union{AbstractVector, Nothing} = nothing,
    record_chr::Union{AbstractVector, Nothing} = nothing,
    record_pos::Union{AbstractVector, Nothing} = nothing,
    record_ids::Union{AbstractVector, Nothing} = nothing,
    record_ref::Union{AbstractVector, Nothing} = nothing,
    record_alt::Union{AbstractVector, Nothing} = nothing
    ) where T <: Real
    msg != "" && (pmeter = Progress(size(A, 2), 5, msg))
    sampleID === nothing || (sampleID .= header(reader).sampleID)
    record = VCF.Record()
    for j in 1:size(A, 2)
        if eof(reader)
            @warn("Reached end of reader; columns $j-$(size(A, 2)) are " *
                "set to missing values")
            fill!(view(A, :, j:size(A, 2)), missing)
            break
        else
            read!(reader, record)
        end
        gtkey = findgenokey(record, "GT")
        # if no GT field, fill by missing values
        if gtkey === nothing
            @inbounds @simd for i in 1:size(A, 1)
                A[i, j] = missing
            end
        end
        # save record's information
        record_chr === nothing || (record_chr[j] = VCF.chrom(record))
        record_pos === nothing || (record_pos[j] = VCF.pos(record))
        record_ids === nothing || (record_ids[j] = try VCF.id(record) catch; ["."] end)
        record_ref === nothing || (record_ref[j] = VCF.ref(record))
        record_alt === nothing || (record_alt[j] = VCF.alt(record))
        # second pass: impute, convert, center, scale
        _, _, _, _, _, alt_freq, _, _, _, _, _ = gtstats(record, nothing)
        ct = 2alt_freq
        wt = alt_freq == 0 ? 1.0 : 1.0 / √(2alt_freq * (1 - alt_freq))
        for i in 1:size(A, 1)
            geno = record.genotype[i]
            # Missing genotype: dropped field or when either haplotype contains "."
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                if impute
                    a1, a2 = rand() ≤ alt_freq, rand() ≤ alt_freq
                    A[i, j] = convert_gt(T, (a1, a2), model)
                else
                    A[i, j] = missing
                end
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[i, j] = convert_gt(T, (a1, a2), model)
            end
            # center and scale if asked
            center && !ismissing(A[i, j]) && (A[i, j] -= ct)
            scale && !ismissing(A[i, j]) && (A[i, j] *= wt)
        end

        # update progress
        msg != "" && next!(pmeter)
    end
    return A
end

"""
    copy_gt_trans!(A, reader; ...)

Fill the rows of matrix `A` by the GT data from VCF records in `reader` where
the ALT allele in each record is interpreted as a `1`. Records with missing
genotypes are converted to `missing`.

# Input
- `A`: a matrix or vector such that `eltype(A) <: Union{Missing, Real}`. Each
    column is one person's genotype. 
- `reader`: a VCF reader

# Optional argument
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`
- `msg`: A message that will be printed to indicate progress. Defaults to
    not printing. 
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 

# Output
- `A`: `ismissing(A[i, j]) == true` indicates missing genotype. If `impute=true`,
    `ismissing(A[i, j]) == false` for all entries.
"""
function copy_gt_trans!(
    A::Union{AbstractMatrix{Union{Missing, T}}, AbstractVector{Union{Missing, T}}},
    reader::VCF.Reader;
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false,
    msg::String = "",
    sampleID::Union{AbstractVector, Nothing} = nothing,
    record_chr::Union{AbstractVector, Nothing} = nothing,
    record_pos::Union{AbstractVector, Nothing} = nothing,
    record_ids::Union{AbstractVector, Nothing} = nothing,
    record_ref::Union{AbstractVector, Nothing} = nothing,
    record_alt::Union{AbstractVector, Nothing} = nothing
    ) where T <: Real
    msg != "" && (pmeter = Progress(size(A, 1), 5, msg))
    sampleID === nothing || (sampleID .= header(reader).sampleID)
    record = VCF.Record()

    for j in 1:size(A, 1)
        if eof(reader)
            @warn("Reached end of reader; rows $j-$(size(A, 1)) are set to " * 
                "missing values")
            fill!(view(A, j:size(A, 1), :), missing)
            break
        else
            read!(reader, record)
        end
        gtkey = findgenokey(record, "GT")
        # if no GT field, fill by missing values
        if gtkey === nothing
            fill!(view(A, j, :), missing)
        end
        # save record's information
        record_chr === nothing || (record_chr[j] = VCF.chrom(record))
        record_pos === nothing || (record_pos[j] = VCF.pos(record))
        record_ids === nothing || (record_ids[j] = try VCF.id(record) catch; ["."] end)
        record_ref === nothing || (record_ref[j] = VCF.ref(record))
        record_alt === nothing || (record_alt[j] = VCF.alt(record))
        # second pass: impute, convert, center, scale
        _, _, _, _, _, alt_freq, _, _, _, _, _ = gtstats(record, nothing)
        ct = 2alt_freq
        wt = alt_freq == 0 ? 1.0 : 1.0 / √(2alt_freq * (1 - alt_freq))
        for i in 1:size(A, 2)
            geno = record.genotype[i]
            # Missing genotype: dropped field or when either haplotype contains "."
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                if impute
                    a1, a2 = rand() ≤ alt_freq, rand() ≤ alt_freq
                    A[j, i] = convert_gt(T, (a1, a2), model)
                else
                    A[j, i] = missing
                end
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[j, i] = convert_gt(T, (a1, a2), model)
            end
            # center and scale if asked
            center && !ismissing(A[j, i]) && (A[j, i] -= ct)
            scale && !ismissing(A[j, i]) && (A[j, i] *= wt)
        end
        # update progress
        msg != "" && next!(pmeter)
    end
    return nothing
end

function geno_ismissing(record::VCF.Record, range::UnitRange{Int})
    return record.data[first(range)] == UInt8('.') || 
           record.data[last(range)] == UInt8('.')
end

"""
    convert_gt(t, vcffile; [impute=false], [center=false], [scale=false])

Convert the GT data from a VCF file to a matrix of type `t`, allowing for 
missing values. Records with missing genotypes are converted to `missing`.

# Input
- `t`: a type `t <: Real`
- `vcffile`: VCF file path

# Optional argument
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`
- `trans`: whether to import data transposed so that each column is 1 genotype,
    default `false`
- `msg`: A message that will be printed to indicate progress. Defaults to
    not printing. 
- `save_snp_info`: Boolean. If true, will also output sample ID, chrom, pos,
    snp id, ref, and alt vectors

# Output
- `out`: The genotype matrix stored in column or rows depending on `trans`.
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 
"""
function convert_gt(
    t::Type{T},
    vcffile::AbstractString;
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false,
    trans::Bool = false,
    msg::String = "", 
    save_snp_info::Bool = false
    ) where T <: Real
    reader = VCF.Reader(openvcf(vcffile, "r"))
    records = nrecords(vcffile)
    samples = nsamples(vcffile)
    sampleID   = (save_snp_info ? Vector{String}(undef, samples) : nothing)
    record_chr = (save_snp_info ? Vector{String}(undef, records) : nothing)
    record_pos = (save_snp_info ? zeros(Int, records) : nothing)
    record_ids = (save_snp_info ? Vector{Vector{String}}(undef, records) : nothing)
    record_ref = (save_snp_info ? Vector{String}(undef, records) : nothing)
    record_alt = (save_snp_info ? Vector{Vector{String}}(undef, records) : nothing)

    if trans
        out = Matrix{Union{t, Missing}}(undef, records, samples)
        copy_gt_trans!(out, reader; model = model, impute = impute,
            center = center, scale = scale, msg = msg, sampleID=sampleID,
            record_chr=record_chr, record_pos=record_pos, record_ids=record_ids,
            record_ref=record_ref, record_alt=record_alt)
    else
        out = Matrix{Union{t, Missing}}(undef, samples, records)
        copy_gt!(out, reader; model = model, impute = impute,
            center = center, scale = scale, msg = msg, sampleID=sampleID,
            record_chr=record_chr, record_pos=record_pos, record_ids=record_ids,
            record_ref=record_ref, record_alt=record_alt)
    end
    close(reader)
    if save_snp_info
        return out, sampleID, record_chr, record_pos, record_ids, record_ref, 
            record_alt
    else
        return out
    end
end

"""
    convert_ht(t, vcffile)

Converts the GT data from a VCF file to a haplotype matrix of type `t`. One 
VCF record will become 2 column/row in the resulting matrix. Records cannot
have missing genotypes and the target file must be phased. 

# Input
- `t`: a type `t <: Real`. If `t` is `Bool`, output will be a BitMatrix.
- `vcffile`: a phased VCF file
- `trans`: whether to import data transposed so that each column is a
    haplotype, default `false`.
- `msg`: A message that will be printed to indicate import progress. Defaults to
    not printing. 
- `save_snp_info`: Boolean. If true, will also output sample ID, chrom, pos,
    snp id, ref, and alt vectors

# Output
- `out`: The haplotype matrix stored in column or rows depending on `trans`.
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 
"""
function convert_ht(
    t::Type{T},
    vcffile::AbstractString;
    trans::Bool = false,
    msg::String = "",
    save_snp_info::Bool = false
    ) where T <: Real
    reader = VCF.Reader(openvcf(vcffile, "r"))
    records = nrecords(vcffile)
    samples = nsamples(vcffile)
    sampleID   = (save_snp_info ? Vector{String}(undef, samples) : nothing)
    record_chr = (save_snp_info ? Vector{String}(undef, records) : nothing)
    record_pos = (save_snp_info ? zeros(Int, records) : nothing)
    record_ids = (save_snp_info ? Vector{Vector{String}}(undef, records) : nothing)
    record_ref = (save_snp_info ? Vector{String}(undef, records) : nothing)
    record_alt = (save_snp_info ? Vector{Vector{String}}(undef, records) : nothing)

    M = (t == Bool ? BitArray{2} : Matrix{t})
    if trans
        out = M(undef, records, 2samples)
        copy_ht_trans!(out, reader, msg = msg, sampleID=sampleID, 
            record_chr=record_chr, record_pos=record_pos, record_ids=record_ids,
            record_ref=record_ref, record_alt=record_alt)
    else
        out = M(undef, 2samples, records)
        copy_ht!(out, reader, msg = msg, sampleID=sampleID,
            record_chr=record_chr, record_pos=record_pos, record_ids=record_ids,
            record_ref=record_ref, record_alt=record_alt)
    end
    close(reader)
    if save_snp_info
        return out, sampleID, record_chr, record_pos, record_ids, record_ref, 
            record_alt
    else
        return out
    end
end

"""
    copy_ht!(A, reader)

Fill 2 columns of `A` by the GT data from VCF records in `reader`, each record
filling 2 columns. Missing data is not allowed and genotypes must be phased. 

# Input
- `A`: matrix where rows are haplotypes. Person `i`'s haplotype are filled in
    rows `2i - 1` and `2i`. 
- `reader`: a VCF reader

# Optional inputs
- `msg`: A message that will be printed to indicate progress. Defaults to
    not printing.
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 
"""
function copy_ht!(
    A::Union{AbstractMatrix{T}, AbstractVector{T}},
    reader::VCF.Reader;
    msg::String = "",
    sampleID::Union{AbstractVector, Nothing} = nothing,
    record_chr::Union{AbstractVector, Nothing} = nothing,
    record_pos::Union{AbstractVector, Nothing} = nothing,
    record_ids::Union{AbstractVector, Nothing} = nothing,
    record_ref::Union{AbstractVector, Nothing} = nothing,
    record_alt::Union{AbstractVector, Nothing} = nothing
    ) where T <: Real

    n, p = size(A)
    nn   = n >>> 1
    sampleID === nothing || (sampleID .= header(reader).sampleID)
    msg != "" && (pmeter = Progress(p, 5, msg)) # update every 5 seconds
    record = VCF.Record()

    for j in 1:p
        if eof(reader)
            @warn("Reached end of record! Columns $j-$p are filled " * 
                "with $(zero(T))s and are NOT haplotypes!")
            fill!(view(A, :, j:size(A, 2)), zero(T))
            break
        else
            read!(reader, record)
        end
        gtkey = findgenokey(record, "GT")

        # haplotype reference files must have GT field
        if gtkey === nothing
            error("Missing GT field for record $j. Reference panels cannot "
                * "have missing data!")
        end

        # save record's information
        record_chr === nothing || (record_chr[j] = VCF.chrom(record))
        record_pos === nothing || (record_pos[j] = VCF.pos(record))
        record_ids === nothing || (record_ids[j] = try VCF.id(record) catch; ["."] end)
        record_ref === nothing || (record_ref[j] = VCF.ref(record))
        record_alt === nothing || (record_alt[j] = VCF.alt(record))

        # second pass: convert
        for i in 1:nn
            geno = record.genotype[i]
            # Missing genotype: dropped field or when either haplotype contains "."
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                error("Missing GT field for record $j entry $(2i - 1). " * 
                    "Reference panels cannot have missing data!")
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[2i - 1, j] = convert(T, a1)
                A[2i    , j] = convert(T, a2)
            end
        end

        # update progress
        msg != "" && next!(pmeter)
    end
    return A
end

"""
    copy_ht_trans!(A, reader)

Fill 2 rows of `A` by the GT data from VCF records in `reader`, each record
filling 2 columns. Missing data is not allowed and genotypes must be phased. 

# Input
- `A`: matrix where columns are haplotypes. Person `i`'s haplotype are filled in
    rows `2i - 1` and `2i`.
- `reader`: a VCF reader

# Optional inputs
- `msg`: A message that will be printed to indicate progress. Defaults to
    not printing.
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 
"""
function copy_ht_trans!(
    A::Union{AbstractMatrix{T}, AbstractVector{T}},
    reader::VCF.Reader;
    msg::String = "",
    sampleID::Union{AbstractVector, Nothing} = nothing,
    record_chr::Union{AbstractVector, Nothing} = nothing,
    record_pos::Union{AbstractVector, Nothing} = nothing,
    record_ids::Union{AbstractVector, Nothing} = nothing,
    record_ref::Union{AbstractVector, Nothing} = nothing,
    record_alt::Union{AbstractVector, Nothing} = nothing
    ) where T <: Real

    p, n = size(A)
    nn   = n >>> 1
    msg != "" && (pmeter = Progress(p, 5, msg)) # update every 5 seconds
    sampleID === nothing || (sampleID .= header(reader).sampleID)
    record = VCF.Record()

    for j in 1:p
        if eof(reader)
            @warn("Reached end of record! Rows $j-$p are filled with "
                * "$(zero(T))s and are NOT haplotypes!")
            fill!(view(A, j:size(A, 2)), zero(T))
            break
        else
            read!(reader, record)
        end
        gtkey = findgenokey(record, "GT")

        # save record's information
        record_chr === nothing || (record_chr[j] = VCF.chrom(record))
        record_pos === nothing || (record_pos[j] = VCF.pos(record))
        record_ids === nothing || (record_ids[j] = try VCF.id(record) catch; ["."] end)
        record_ref === nothing || (record_ref[j] = VCF.ref(record))
        record_alt === nothing || (record_alt[j] = VCF.alt(record))

        # haplotype reference files must have GT field
        if gtkey === nothing
            error("Missing GT field for record $j. Reference panels cannot" *
                " have missing data!")
        end

        # second pass: convert
        for i in 1:nn
            geno = record.genotype[i]
            # Missing genotype: dropped field or when either haplotype contains "."
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                error("Missing GT field for record $j entry $(2i - 1). " *
                    "Reference panels cannot have missing data!")
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[j, 2i - 1] = convert(T, a1)
                A[j, 2i] = convert(T, a2)
            end
        end

        # update progress
        msg != "" && next!(pmeter)
    end
    return nothing
end

"""
    convert_ds(t, vcffile; [key = "DS"], ...)

Converts dosage data from a VCF file to a numeric matrix of type `t`. Here `key`
specifies the FORMAT field of the VCF file that encodes the dosage
(default = "DS"). 

# Arguments
- `t`: type of output matrix. 
- `vcffile`: VCF file path

# Optional argument
- `key`: The FIELD name of the VCF file that encodes dosages, defaults to `"DS"`
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype with 2 times the ALT allele frequency or
    not. Default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`
- `trans`: whether to import data transposed so that each column is 1 genotype,
    default `false`
- `msg`: A message that will be printed to indicate progress. Defaults to
    not printing. 

# Output
- `out`: Dosage matrix 
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 
"""
function convert_ds(
    t::Type{T},
    vcffile::AbstractString;
    key::String = "DS",
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false,
    trans::Bool = false,
    save_snp_info::Bool=false,
    msg::String = ""
    ) where T <: Real
    reader = VCF.Reader(openvcf(vcffile, "r"))
    records = nrecords(vcffile)
    samples = nsamples(vcffile)
    sampleID   = (save_snp_info ? Vector{String}(undef, samples) : nothing)
    record_chr = (save_snp_info ? Vector{String}(undef, records) : nothing)
    record_pos = (save_snp_info ? zeros(Int, records) : nothing)
    record_ids = (save_snp_info ? Vector{Vector{String}}(undef, records) : nothing)
    record_ref = (save_snp_info ? Vector{String}(undef, records) : nothing)
    record_alt = (save_snp_info ? Vector{Vector{String}}(undef, records) : nothing)

    if trans
        out = Matrix{Union{Missing, t}}(undef, records, samples)
        copy_ds_trans!(out, reader, key = key, impute = impute, center = center, 
            scale = scale, msg = msg, sampleID=sampleID, record_chr=record_chr, 
            record_pos=record_pos, record_ids=record_ids, record_ref=record_ref, 
            record_alt=record_alt)
    else
        out = Matrix{Union{Missing, t}}(undef, samples, records)
        copy_ds!(out, reader, key = key, impute = impute, center = center, 
            scale = scale, msg = msg, sampleID=sampleID, record_chr=record_chr, 
            record_pos=record_pos, record_ids=record_ids, record_ref=record_ref, 
            record_alt=record_alt)
    end
    close(reader)
    if save_snp_info
        return out, sampleID, record_chr, record_pos, record_ids, record_ref,
            record_alt
    else
        return out
    end
end

"""
    copy_ds!(A, reader; [key="DS"], ...)

Fill the columns of matrix `A` by the dosage data from VCF records in `reader`.
Records with missing fields (i.e. `key` field is `.`) are converted to
`missing`.

# Input
- `A`: matrix where rows are dosages. 
- `reader`: a VCF reader

# Optional argument
- `key`: The field key for accessing dosage values, default `"DS"`
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`
- `msg`: A message that will be printed to indicate progress. Defaults to
    not printing.
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 
"""
function copy_ds!(
    A::Union{AbstractMatrix{Union{Missing, T}}, AbstractVector{Union{Missing, T}}},
    reader::VCF.Reader;
    key::String = "DS",
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false,
    msg::String = "",
    sampleID::Union{AbstractVector, Nothing} = nothing,
    record_chr::Union{AbstractVector, Nothing} = nothing,
    record_pos::Union{AbstractVector, Nothing} = nothing,
    record_ids::Union{AbstractVector, Nothing} = nothing,
    record_ref::Union{AbstractVector, Nothing} = nothing,
    record_alt::Union{AbstractVector, Nothing} = nothing
    ) where T <: Real
    msg != "" && (pmeter = Progress(size(A, 2), 5, msg)) # update every 5 seconds
    sampleID === nothing || (sampleID .= header(reader).sampleID)
    record = VCF.Record()
    for j in 1:size(A, 2)
        if eof(reader)
            @warn("Reached end of reader; columns $j-$(size(A, 2)) are set " * 
                "to missing values")
            fill!(view(A, :, j:size(A, 2)), missing)
            break
        else
            read!(reader, record)
        end
        dskey = findgenokey(record, key)

        # if no dosage field, fill by missing values
        if dskey === nothing
            @inbounds @simd for i in 1:size(A, 1)
                A[i, j] = missing
            end
        end

        # save record's information
        record_chr === nothing || (record_chr[j] = VCF.chrom(record))
        record_pos === nothing || (record_pos[j] = VCF.pos(record))
        record_ids === nothing || (record_ids[j] = try VCF.id(record) catch; ["."] end)
        record_ref === nothing || (record_ref[j] = VCF.ref(record))
        record_alt === nothing || (record_alt[j] = VCF.alt(record))

        # loop over every marker in record
        # _, _, _, _, _, alt_freq, _, _, _, _, _ = gtstats(record, nothing) 
        # ct = 2alt_freq
        # wt = alt_freq == 0 ? 1.0 : 1.0 / √(2alt_freq * (1 - alt_freq))
        for i in 1:size(A, 1)
            if dskey === nothing
                break
            end
            geno = record.genotype[i]
            # Missing genotype: dropped field or "." => 0x2e
            if dskey > lastindex(geno) || record.data[geno[dskey]] == [0x2e]
                A[i, j] = missing # (impute ? ct : missing)
            else # not missing
                A[i, j] = parse(T, String(record.data[geno[dskey]]))
            end

        end
        if center || scale || impute
            total_ds = zero(T)
            cnt = 0
            for i in 1:size(A, 1)
                if A[i, j] !== missing
                    total_ds += A[i, j]
                    cnt += 1
                end
            end
            ct = total_ds / cnt
            wt = ct == 0 ? 1.0 : 1.0 / √(ct * (1 - ct/2))
            # center and scale if asked
            center && !ismissing(A[i, j]) && (A[i, j] -= ct)
            scale && !ismissing(A[i, j]) && (A[i, j] *= wt)
            if impute
                for i in 1:size(A, 1)
                    if A[i, j] === missing
                        A[i, j] = ct
                    end
                end
            end
        end
        # update progress
        msg != "" && next!(pmeter)
    end
end

"""
    copy_ds_trans!(A, reader; [key="DS"], ...)

Fill the rows of matrix `A` by the dosage data from VCF records in `reader`.
Records with missing fields (i.e. `key` field is `.`) are converted to `missing`.

# Input
- `A`: matrix where columns are dosages. 
- `reader`: a VCF reader

# Optional argument
- `key`: The field key for accessing dosage values, default `"DS"`
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`
- `msg`: A message that will be printed to indicate progress. Defaults to
    not printing. 
- `sampleID`: The sample ID for each person
- `record_chr`: The chromosome number for each record
- `record_pos`: The position for each record
- `record_ids`: The record ID for each record. Only the first one is
    recorded when there are multiple IDs. 
- `record_ref`: The ref allele for each record
- `record_alt`: The alternate allele for each record. Only the first one is
    recorded when there are multiple alternate alleles. 
"""
function copy_ds_trans!(
    A::Union{AbstractMatrix{Union{Missing, T}}, AbstractVector{Union{Missing, T}}},
    reader::VCF.Reader;
    key::String = "DS",
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false,
    msg::String = "",
    sampleID::Union{AbstractVector, Nothing} = nothing,
    record_chr::Union{AbstractVector, Nothing} = nothing,
    record_pos::Union{AbstractVector, Nothing} = nothing,
    record_ids::Union{AbstractVector, Nothing} = nothing,
    record_ref::Union{AbstractVector, Nothing} = nothing,
    record_alt::Union{AbstractVector, Nothing} = nothing
    ) where T <: Real
    msg != "" && (pmeter = Progress(size(A, 1), 5, msg)) # update every 5 seconds
    sampleID === nothing || (sampleID .= header(reader).sampleID)
    record = VCF.Record()
    for j in 1:size(A, 1)
        if eof(reader)
            @warn("Reached end of reader; rows $j-$(size(A, 1)) " * 
                "are set to missing values")
            fill!(view(A, j:size(A, 2), :), missing)
            break
        else
            read!(reader, record)
        end
        dskey = findgenokey(record, key)

        # if no dosage field, fill by missing values
        if dskey === nothing
            @inbounds @simd for i in 1:size(A, 2)
                A[j, i] = missing
            end
        end

        # save record's information
        record_chr === nothing || (record_chr[j] = VCF.chrom(record))
        record_pos === nothing || (record_pos[j] = VCF.pos(record))
        record_ids === nothing || (record_ids[j] = try VCF.id(record) catch; ["."] end)
        record_ref === nothing || (record_ref[j] = VCF.ref(record))
        record_alt === nothing || (record_alt[j] = VCF.alt(record))
        
        # loop over every marker in record
        # _, _, _, _, _, alt_freq, _, _, _, _, _ = gtstats(record, nothing) 
        # ct = 2alt_freq
        # wt = alt_freq == 0 ? 1.0 : 1.0 / √(2alt_freq * (1 - alt_freq))
        for i in 1:size(A, 2)
            geno = record.genotype[i]
            # Missing genotype: dropped field or "." => 0x2e
            if dskey > lastindex(geno) || record.data[geno[dskey]] == [0x2e]
                A[j, i] = missing #(impute ? ct : missing)
            else # not missing
                A[j, i] = parse(T, String(record.data[geno[dskey]]))
            end
        end
        if center || scale || impute
            total_ds = zero(T)
            cnt = 0
            for i in 1:size(A, 2)
                if A[j, i] !== missing
                    total_ds += A[j, i]
                    cnt += 1
                end
            end
            ct = total_ds / cnt
            wt = ct == 0 ? 1.0 : 1.0 / √(ct * (1 - ct/2))
            # center and scale if asked
            center && !ismissing(A[j, i]) && (A[j, i] -= ct)
            scale && !ismissing(A[j, i]) && (A[j, i] *= wt)
            if impute
                for i in 1:size(A, 2)
                    if A[j, i] === missing
                        A[j, i] = ct
                    end
                end
            end
        end

        #update progress
        msg != "" && next!(pmeter)
    end
end
