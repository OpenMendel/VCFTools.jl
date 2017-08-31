using NullableArrays
import SnpArrays.isnan, SnpArrays.convert, SnpArrays.randgeno

"""
    copy_gt!(A, reader; [impute=false], [center=false], [scale=false])

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
- `A`: `isnull(A[i, j]) == true` indicates missing genotype, even when
    `A.values[i, j]` may hold the imputed genotype
"""
function copy_gt!(
    A::Union{NullableMatrix{T}, NullableVector{T}},
    reader::VCF.Reader;
    model::Symbol = :additive,
    impute::Bool = false,
    center::Bool = false,
    scale::Bool = false
    ) where T <: Real
    for j in 1:size(A, 2)
        if eof(reader)
            warning("Only $j records left in reader; $(j+1)-th to last column are set to missing values")
            fill!(view(A, :, (j + 1):size(A, 2)), Nullable(zero(T), false))
            break
        else
            record = read(reader)
        end
        gtkey = VCF.findgenokey(record, "GT")
        # if no GT field, fill by missing values
        if gtkey == 0
            @inbounds @simd for i in 1:size(A, 1)
                A[i, j] = Nullable(zero(T), false)
            end
        end
        # convert GT field to numbers according to specified genetic model
        _, _, _, _, _, _, _, _, minor_allele, maf, _ = gtstats(record, nothing)
        # second pass: impute, convert, center, scale
        ct = 2maf
        wt = maf == 0 ? 1.0 : 1.0 / √(2maf * (1 - maf))
        @simd for i in 1:size(A, 1)
            geno = record.genotype[i]
            # dropped field or "." => 0x2e
            if gtkey > endof(geno) || record.data[geno[gtkey]] == [0x2e]
                A[i, j] = Nullable(zero(T), false)
            else
                # "0" (ALT) => 0x30, "1" (REF) => 0x31
                # In SnpArray: A1 = ALT, A2 = REF
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                if a1 && !a2 # (true, false) in SnpArrays is preserved for missing genotype
                    a1, a2 = (false, true)
                end
            end
            # impute if asked
            if isnan(a1, a2) && impute
                a1, a2 = randgeno(maf, minor_allele)
            end
            A[i, j] = Nullable(convert(T, (a1, a2), !minor_allele, model), true)
            # center and scale if asked
            center && (A.values[i, j] -= ct)
            scale && (A.values[i, j] *= wt)
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
    out = NullableArray(t, nsamples(vcffile), nrecords(vcffile))
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_gt!(out, reader; model = model, impute = impute,
        center = center, scale = scale)
    out
end
