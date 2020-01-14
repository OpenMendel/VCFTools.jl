"""
Convert a two-bit genotype to a real number (minor allele count) of type
`t` according to specified SNP model. Missing genotype is converted
to null. `minor_allele==true` indicates `REF` is the minor allele;
`minor_allele==false` indicates `ALT` is the minor allele.
"""
function convert_gt(
    t::Type{T},
    a::NTuple{2, Bool},
    minor_allele::Bool,
    model::Symbol = :additive
    ) where T <: Real
    if minor_allele # REF is the minor allele
        if model == :additive
            return convert(T, a[1] + a[2])
        elseif model == :dominant
            return convert(T, a[1] | a[2])
        elseif model == :recessive
            return convert(T, a[1] & a[2])
        else
            throw(ArgumentError("un-recognized SNP model: $model"))
        end
    else # ALT is the minor allele
        if model == :additive
            return convert(T, !a[1] + !a[2])
        elseif model == :dominant
            return convert(T, !a[1] | !a[2])
        elseif model == :recessive
            return convert(T, !a[1] & !a[2])
        else
            throw(ArgumentError("un-recognized SNP model: $model"))
        end
    end
end

"""
    copy_gt!(A, reader; [model=:additive], [impute=false], [center=false], [scale=false])

Fill the columns of matrix `A` by the GT data from VCF records in `reader` where the minor 
allele for each record is computed on-the-fly. Each column of `A` corresponds to one record. 
Record without GT field is converted to `missing`.

# Input
- `A`: a matrix or vector such that `eltype(A) <: Union{Missing, Real}`
- `reader`: a VCF reader

# Optional argument
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`

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
    scale::Bool = false
    ) where T <: Real
    for j in 1:size(A, 2)
        if eof(reader)
            @warn("Only $j records left in reader; columns $(j+1)-$(size(A, 2)) are set to missing values")
            fill!(view(A, :, (j + 1):size(A, 2)), missing)
            break
        else
            record = read(reader)
        end
        gtkey = VCF.findgenokey(record, "GT")
        # if no GT field, fill by missing values
        if gtkey == nothing
            @inbounds @simd for i in 1:size(A, 1)
                A[i, j] = missing
            end
        end
        # convert GT field to numbers according to specified genetic model
        _, _, _, _, _, _, _, _, minor_allele, maf, _ = gtstats(record, nothing)
        # second pass: impute, convert, center, scale
        ct = 2maf
        wt = maf == 0 ? 1.0 : 1.0 / √(2maf * (1 - maf))
        for i in 1:size(A, 1)
            geno = record.genotype[i]
            # Missing genotype: dropped field or when either haplotype contains "."
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                if impute
                    if minor_allele # REF is the minor allele
                        a1, a2 = rand() ≤ maf, rand() ≤ maf
                    else # ALT is the minor allele
                        a1, a2 = rand() > maf, rand() > maf
                    end
                    A[i, j] = convert_gt(T, (a1, a2), minor_allele, model)
                else
                    A[i, j] = missing
                end
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[i, j] = convert_gt(T, (a1, a2), minor_allele, model)
            end
            # center and scale if asked
            center && !ismissing(A[i, j]) && (A[i, j] -= ct)
            scale && !ismissing(A[i, j]) && (A[i, j] *= wt)
        end
    end
    A
end

"""
    copy_gt_as_is!(A, reader; [model=:additive], [impute=false], [center=false], [scale=false])

Fill the columns of matrix `A` by the GT data from VCF records in `reader` where the ALT allele
in each record is interpreted as a `1`. Each column of `A` corresponds to one record. 
Record without GT field is converted to `missing`.

# Input
- `A`: a matrix or vector such that `eltype(A) <: Union{Missing, Real}`
- `reader`: a VCF reader

# Optional argument
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`

# Output
- `A`: `ismissing(A[i, j]) == true` indicates missing genotype. If `impute=true`,
    `ismissing(A[i, j]) == false` for all entries.
"""
function copy_gt_as_is!(
    A::Union{AbstractMatrix{Union{Missing, T}}, AbstractVector{Union{Missing, T}}},
    reader::VCF.Reader;
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false
    ) where T <: Real
    for j in 1:size(A, 2)
        if eof(reader)
            @warn("Only $j records left in reader; columns $(j+1)-$(size(A, 2)) are set to missing values")
            fill!(view(A, :, (j + 1):size(A, 2)), missing)
            break
        else
            record = read(reader)
        end
        gtkey = VCF.findgenokey(record, "GT")
        # if no GT field, fill by missing values
        if gtkey == nothing
            @inbounds @simd for i in 1:size(A, 1)
                A[i, j] = missing
            end
        end
        # second pass: impute, convert, center, scale
        _, _, _, _, _, maf, _, _, _, _, _ = gtstats(record, nothing) #reffreq in gtstats actually count number of "1"
        ct = 2maf
        wt = maf == 0 ? 1.0 : 1.0 / √(2maf * (1 - maf))
        for i in 1:size(A, 1)
            geno = record.genotype[i]
            # Missing genotype: dropped field or when either haplotype contains "."
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                if impute
                    a1, a2 = rand() ≤ maf, rand() ≤ maf
                    A[i, j] = convert_gt(T, (a1, a2), true, model)
                else
                    A[i, j] = missing
                end
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[i, j] = convert(T, a1 + a2)
            end
            # center and scale if asked
            center && !ismissing(A[i, j]) && (A[i, j] -= ct)
            scale && !ismissing(A[i, j]) && (A[i, j] *= wt)
        end
    end
    return A
end

function geno_ismissing(record::VCF.Record, range::UnitRange{Int})
    return record.data[first(range)] == UInt8('.') || record.data[last(range)] == UInt8('.')
end

"""
    convert_gt!(t, vcffile; [impute=false], [center=false], [scale=false])

Convert the GT data from a VCF file to a matrix of type `t`, allowing for 
missing values. Each column of the matrix corresponds to one VCF record. 
Record without GT field is converted to equivalent of missing genotypes.

# Input
- `t`: a type `t <: Real`
- `vcffile`: VCF file path

# Optional argument
- `as_minorallele`: convert VCF data (1 indicating ALT allele) to minor allele count or not. If `false`, 0 and 1 will be read as stored in VCF file, default `false`. 
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`

# Output
- `A`: matrix where `eltype(A) <: Union{missing, Real}`. `ismissing(A[i, j]) == true`
    indicates missing genotype.
"""
function convert_gt(
    t::Type{T},
    vcffile::AbstractString;
    as_minorallele::Bool = false,
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false,
    ) where T <: Real
    out = Matrix{Union{t, Missing}}(undef, nsamples(vcffile), nrecords(vcffile))
    reader = VCF.Reader(openvcf(vcffile, "r"))
    if as_minorallele
        copy_gt!(out, reader; model = model, impute = impute,
            center = center, scale = scale)
    else
        copy_gt_as_is!(out, reader; model = model, impute = impute,
            center = center, scale = scale)
    end
    close(reader)
    out
end

"""
    convert_ht(t, vcffile; [as_minorallele=false])

Converts the GT data from a VCF file to a haplotype matrix of type `t`. One 
column of the VCF record will become 2 column in the resulting matrix. If haplotypes
are not phased, 1 will be on the left column (i.e. `1/0`).

# Input
- `t`: a type `t <: Real`
- `vcffile`: VCF file path

# Optional argument
- `as_minorallele`: convert VCF data (1 indicating ALT allele) to minor allele count or not. If `false`, 0 and 1 will be read as stored in VCF file, default `false`. 
"""
function convert_ht(
    t::Type{T},
    vcffile::AbstractString;
    as_minorallele::Bool = false
    ) where T <: Real
    out = Matrix{t}(undef, 2nsamples(vcffile), nrecords(vcffile))
    reader = VCF.Reader(openvcf(vcffile, "r"))
    if as_minorallele
        copy_ht!(out, reader)
    else
        copy_ht_as_is!(out, reader)
    end
    close(reader)
    return out
end

"""
Convert a one-bit haplotype to a real number (minor allele count) of type
`t`. Here `minor_allele==true` indicates `REF` is the minor allele;
`minor_allele==false` indicates `ALT` is the minor allele.
"""
function convert_ht(
    t::Type{T},
    a::Bool,
    minor_allele::Bool
    ) where T <: Real
    if minor_allele # REF is the minor allele
        return convert(T, a)
    else # ALT is the minor allele
        return convert(T, !a)
    end
end

"""
    copy_ht!(A, reader)

Fill 2 columns of BitArray `A` by the GT data from VCF records in `reader`, each record filling 2 columns.
The minor allele for each record is computed on-the-fly. Record without GT field is converted to `missing`.
Missing GT field is NOT allowed. 

# Input
- `A`: a BitArray{2} or BitArray{1}. 
- `reader`: a VCF reader
"""
function copy_ht!(
    A::Union{AbstractMatrix{T}, AbstractVector{T}},
    reader::VCF.Reader
    ) where T <: Real

    n, p = size(A)
    nn   = Int(n / 2)

    # loop over every record, each filling 2 columns of A
    for j in 1:p
        if eof(reader)
            @warn("Reached end! Check if output is correct!")
            A = A[:, 1:2j]
            break
        else
            record = read(reader)
        end
        gtkey = VCF.findgenokey(record, "GT")

        # haplotype reference files must have GT field
        if gtkey == nothing
            error("Missing GT field for record $j. Reference panels cannot have missing data!")
        end

        # check whether allele is REF or ALT
        _, _, _, _, _, _, _, _, minor_allele, _, _ = gtstats(record, nothing)

        # loop over every marker in record
        for i in 1:nn
            geno = record.genotype[i]
            # Missing genotype: dropped field or when either haplotype contains "."
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                error("Missing GT field for record $j entry $i. Reference panels cannot have missing data!")
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[2i - 1, j] = convert_ht(T, a1, minor_allele)
                A[2i    , j] = convert_ht(T, a2, minor_allele)
            end
        end
    end
end

"""
    copy_ht_as_is!(A, reader)

Fill 2 columns of BitArray `A` by the GT data from VCF records in `reader`, each record filling 2 columns.
The minor allele rfor each ecord is the ALT allele (i.e. will read 0s and 1s of the vcffile as-is). Record 
without GT field is converted to `missing`. Missing GT field is NOT allowed. 

# Input
- `A`: a BitArray{2} or BitArray{1}. 
- `reader`: a VCF reader
"""
function copy_ht_as_is!(
    A::Union{AbstractMatrix{T}, AbstractVector{T}},
    reader::VCF.Reader;
    ) where T <: Real

    n, p = size(A)
    nn   = Int(n / 2)

    for j in 1:size(A, 2)
        if eof(reader)
            @warn("Reached end! Check if output is correct!")
            A = A[:, 1:2j]
            break
        else
            record = read(reader)
        end
        gtkey = VCF.findgenokey(record, "GT")

        # haplotype reference files must have GT field
        if gtkey == nothing
            error("Missing GT field for record $j. Reference panels cannot have missing data!")
        end

        # second pass: convert
        for i in 1:nn
            geno = record.genotype[i]
            # Missing genotype: dropped field or when either haplotype contains "."
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                error("Missing GT field for record $j entry $i. Reference panels cannot have missing data!")
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[2i - 1, j] = convert(T, a1)
                A[2i    , j] = convert(T, a2)
            end
        end
    end
    return A
end

"""
    convert_ds(t, vcffile; [key = "DS"], [impute=false], [center=false], [scale=false])

Converts dosage data from a VCF file to a numeric matrix of type `t`. Here `key` specifies
the FORMAT field of the VCF file that encodes the dosage (default = "DS"). If `impute=true`,
we fill missing entries with 2 times the ALT allele frequency. 

# Arguments
- `t`: type of output matrix. 

# Optional argument
- `as_minorallele`: convert VCF data (1 indicating ALT allele) to minor allele count or not. If `false`, 0 and 1 will be read as stored in VCF file, default `false`. 
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`
"""
function convert_ds(
    t::Type{T},
    vcffile::AbstractString;
    key::String = "DS",
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false
    ) where T <: Real
    out = Matrix{Union{Missing, t}}(undef, nsamples(vcffile), nrecords(vcffile))
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_ds!(out, reader, key = key, impute = impute, center = center, scale = scale)
    close(reader)
    return out
end

"""
    copy_ds!(A, reader; [key="DS"], [model=:additive], [impute=false], [center=false], [scale=false])

Fill the columns of matrix `A` by the dosage data from VCF records in `reader`. Record without GT field 
is converted to `missing`. 

# Input
- `A`: a matrix or vector such that `eltype(A) <: Union{Missing, Real}`
- `reader`: a VCF reader

# Optional argument
- `key`: The field key for accessing dosage values, default `"DS"`
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`

# Output
- `A`: `ismissing(A[i, j]) == true` indicates missing genotype. If `impute=true`,
    `ismissing(A[i, j]) == false` for all entries.
"""
function copy_ds!(
    A::Union{AbstractMatrix{Union{Missing, T}}, AbstractVector{Union{Missing, T}}},
    reader::VCF.Reader;
    key::String = "DS",
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false
    ) where T <: Real
    for j in 1:size(A, 2)
        if eof(reader)
            @warn("Only $j records left in reader; columns $(j+1)-$(size(A, 2)) are set to missing values")
            fill!(view(A, :, (j + 1):size(A, 2)), missing)
            break
        else
            record = read(reader)
        end
        dskey = VCF.findgenokey(record, key)

        # if no dosage field, fill by missing values
        if dskey == nothing
            @inbounds @simd for i in 1:size(A, 1)
                A[i, j] = missing
            end
        end

        # loop over every marker in record
        _, _, _, _, _, maf, _, _, _, _, _ = gtstats(record, nothing) #reffreq in gtstats actually count number of "1"
        ct = 2maf
        wt = maf == 0 ? 1.0 : 1.0 / √(2maf * (1 - maf))
        for i in 1:size(A, 1)
            geno = record.genotype[i]
            # Missing genotype: dropped field or "." => 0x2e
            if dskey > lastindex(geno) || record.data[geno[dskey]] == [0x2e]
                A[i, j] = (impute ? ct : missing)
            else # not missing
                A[i, j] = parse(T, String(record.data[geno[dskey]]))
            end
            # center and scale if asked
            center && !ismissing(A[i, j]) && (A[i, j] -= ct)
            scale && !ismissing(A[i, j]) && (A[i, j] *= wt)
        end
    end
end
