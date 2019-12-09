"""
    filter(src, record_index, sample_index; des = "filtered." * src)

Filter vcf files (.vcf or .vcf.gz) according to row indices `record_index` 
and column indices `sample_index` and write to a new set of vcf files `des`.

# Input:
- `src`: full vcf file name
- `record_index`: row indices to keep.
- `sample_index`: column indices to keep.

# Optional arguments:
- `des`: output vcf file name; default it `"filtered." * src`.
"""
function filter(
    src::AbstractString, 
    record_index::AbstractVector{<:Integer}, 
    sample_index::AbstractVector{<:Integer}; 
    des::AbstractString = "filtered." * src 
    )
    # create record (row) and sample (column) masks
    records, samples = nrecords(src), nsamples(src)
    if eltype(record_index) == Bool
        record_mask = record_index
    else
        record_mask = falses(records)
        record_mask[record_index] .= true
    end
    if eltype(sample_index) == Bool
        sample_mask = sample_index
    else
        sample_mask = falses(samples)
        sample_mask[sample_index] .= true
    end

    # create input and output VCF files
    reader = VCF.Reader(openvcf(src, "r"))
    writer = VCF.Writer(openvcf(des, "w"), filter_header(reader, sample_mask))

    # write to des
    for (i, record) in enumerate(reader)
        if record_index[i] 
            filter_record!(record, sample_mask)
            VCF.write(writer, record)
        end
    end

    close(reader)
    flush(writer)
    close(VCF.BioCore.IO.stream(writer))
end

function filter_header(
    reader::VCF.Reader,
    sample_mask::BitVector
    )
    # save meta information
    fileformat = VCF.header(reader).metainfo[1]
    filedate   = VCF.MetaInfo("##filedate=" * string(Dates.today()))
    signature  = VCF.MetaInfo("##source=VCFTools.jl")
    metainfo = vcat(fileformat, filedate, signature, VCF.header(reader).metainfo[2:end])

    # filter sampleID
    sampleID = VCF.header(reader).sampleID[sample_mask]
    return VCF.Header(metainfo, sampleID)
end

"""
    filter_record!(record, sample_mask)
"""
function filter_record!(
    record::VCF.Record,
    sample_mask::BitVector
    )

    # quick return
    if all(sample_mask)
        return nothing
    end

    p = length(record.genotype)
    new_data = UInt8[]
    new_genotype = Vector{UnitRange{Int64}}[]
    sizehint!(new_data, p)
    sizehint!(new_genotype, p)

    # copy chrom, pos, id, ref, alt, qual, filter, information, format into new_data
    for i in 1:(record.genotype[1][1][1] - 1)
        push!(new_data, record.data[i])
    end

    # copy genotype indices and data 
    for i in 1:p
        if sample_mask[i]
            old_geno = record.genotype[i]
            new_geno = UnitRange{Int64}[]
            for g in old_geno
                # save genotypes
                new_start = length(new_data) + 1
                new_end   = new_start + length(g) - 1
                push!(new_geno, new_start:new_end)

                # save strings to data
                old_start = g[1]
                old_end   = g[end]
                [push!(new_data, record.data[i]) for i in old_start:old_end]
                push!(new_data, 0x3a) # 0x3a = byte equivalent of char ':'
            end
            new_data[end] = 0x09 # turn last ':' into '\t'
            push!(new_genotype, new_geno)
        end
    end
    resize!(new_data, length(new_data) - 1) # get rid of last '\t'

    # update pointers to data and genotype indices
    record.data = new_data
    record.genotype = new_genotype
    record.filled = 1:length(new_data)
end

"""
    mask_gt(src, masks; des = "masked." * src, separator = '/')

Creates a new VCF file `des` where genotype entry (i, j) of `src` 
is missing if `masks[i, j]` is true. `src` is unchanged.

# Arguments
- `src`: Input VCF file name.
- `masks`: Bit matrix. `masks[i, j] = true` means mask entry (i, j).
- `des`: output VCF file name.
- `separator`: Separator of VCF genotypes. Can be '/' (default) or '|'.  

# Useful byte reprensetations
- '\t': String([0x09])
- ':': String([0x3a])
- '0': String([0x30])
- '1': String([0x31])
- '.': String([0x2e])
- '|': String([0x7c])
- '/': String([0x2f])

- 'A': String([0x41])
- 'T': String([0x54])
- 'C': String([0x43])
- 'G': String([0x47])
- 'N': String([0x4e])
"""
function mask_gt(
    src::AbstractString, 
    masks::BitArray{2}; 
    des::AbstractString = "masked." * src,
    separator::Char = '/'
    )
    records, samples = nrecords(src), nsamples(src)
    p, n = size(masks)
    if !(records == p && samples == n)
        throw(DimensionMismatch("size(src) = ($records, $samples) ≠ size(masks) = ($p, $n)."))
    end

    # define byte representation of separator
    if separator == '/' 
        separator = 0x2f
    elseif separator == '|'
        separator = 0x7c
    else
        throw(ArgumentError("separator must be / or | but got $separator."))
    end

    # create input and output VCF files
    reader = VCF.Reader(openvcf(src, "r"))
    writer = VCF.Writer(openvcf(des, "w"), VCF.header(reader))

    # loop over each record
    for (i, record) in enumerate(reader)
        gtkey = VCF.findgenokey(record, "GT")
        if gtkey == nothing 
            # loop over genotypes
            for (j, geno) in enumerate(record.genotype)
                if masks[i, j]
                    record.data[geno[gtkey][1]] = 0x2e
                    record.data[geno[gtkey][2]] = separator
                    record.data[geno[gtkey][3]] = 0x2e
                end
            end
        end
        write(writer, record)
    end
    flush(writer); close(reader); close(writer)
end
