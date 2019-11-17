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

Fill the columns of a nullable matrix `A` by the GT data from VCF records in
`reader`. Each column of `A` corresponds to one record. Record without GT field
is converted to `NaN`.

# Input
- `A`: a nullable matrix or nullable vector
- `reader`: a VCF reader

# Optional argument
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`

# Output
- `A`: `isnull(A[i, j]) == true` indicates missing genotype. If `impute=true`,
    `isnull(A[i, j]) == false` for all entries.
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
        if gtkey == 0
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
            # Missing genotype: dropped field or "." => 0x2e
            if gtkey > lastindex(geno) || record.data[geno[gtkey]] == [0x2e]
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
                # "0" (ALT) => 0x30, "1" (REF) => 0x31
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
    convert_gt!(t, vcffile; [impute=false], [center=false], [scale=false])

Convert the GT data from a VCF file to a nullable matrix of type `t`. Each
column of the matrix corresponds to one VCF record. Record without GT field
is converted to equivalent of missing genotypes.

# Input
- `t`: a type `t <: Real`
- `vcffile`: VCF file path

# Optional argument
- `model`: genetic model `:additive` (default), `:dominant`, or `:recessive`
- `impute`: impute missing genotype or not, default `false`
- `center`: center gentoype by 2maf or not, default `false`
- `scale`: scale genotype by 1/√2maf(1-maf) or not, default `false`

# Output
- `A`: a nulalble matrix of type `NullableMatrix{T}`. `isnull(A[i, j]) == true`
    indicates missing genotype, even when `A.values[i, j]` may hold the imputed
    genotype
"""
function convert_gt(
    t::Type{T},
    vcffile::AbstractString;
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false
    ) where T <: Real
    out = Matrix{Union{t, Missing}}(undef, nsamples(vcffile), nrecords(vcffile))
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_gt!(out, reader; model = model, impute = impute,
        center = center, scale = scale)
    close(reader)
    out
end

"""
    convert_ht(t, vcffile)

Converts the GT data from a VCF file to a haplotype matrix of type `t`. One 
column of the VCF record will become 2 column in the resulting matrix. If haplotypes
are not phased, monomorphic alleles will have a 1 on the left column.
"""
function convert_ht(
    t::Type{T},
    vcffile::AbstractString
    ) where T <: Real
    out = Matrix{Union{t, Missing}}(undef, nsamples(vcffile), 2nrecords(vcffile))
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_ht!(out, reader)
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

function copy_ht!(
    A::Union{AbstractMatrix{Union{Missing, T}}, AbstractVector{Union{Missing, T}}},
    reader::VCF.Reader
    ) where T <: Real

    n, p = size(A)
    pp   = Int(p / 2)

    # loop over every record, each filling 2 columns of A
    for j in 1:pp
        if eof(reader)
            @warn("Reached end! Check if output is correct!")
            A = A[:, 1:2j]
        else
            record = read(reader)
        end
        gtkey = VCF.findgenokey(record, "GT")

        # haplotype reference files must have GT field
        if gtkey == 0
            error("Missing GT field for record $j. Reference panels cannot have missing data!")
        end

        # check whether allele is REF or ALT
        _, _, _, _, _, _, _, _, minor_allele, _, _ = gtstats(record, nothing)

        # loop over every marker in record
        for i in 1:n
            geno = record.genotype[i]
            # Missing genotype: dropped field or "." => 0x2e
            if gtkey > lastindex(geno) || record.data[geno[gtkey]] == [0x2e]
                error("Missing GT field for record $j entry $i. Reference panels cannot have missing data!")
            else # not missing
                # "0" (ALT) => 0x30, "1" (REF) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                A[i, 2j - 1] = convert_ht(T, a1, minor_allele)
                A[i, 2j]     = convert_ht(T, a2, minor_allele)
            end
        end
    end
end
