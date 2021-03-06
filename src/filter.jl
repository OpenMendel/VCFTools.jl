# Useful byte reprensetations
# - '\t': String([0x09])
# - '\n': String([0x0a])
# - '#': String([0x23])
# - '\\': String([0x5c])
# - '\r' String([0x0d])
# - ':': String([0x3a])
# - '0': String([0x30])
# - '1': String([0x31])
# - '9': String([0x39])
# - '.': String([0x2e])
# - '|': String([0x7c])
# - '/': String([0x2f])

# - 'A': String([0x41])
# - 'T': String([0x54])
# - 'C': String([0x43])
# - 'G': String([0x47])
# - 'N': String([0x4e])
"""
    filter(src, record_index, sample_index; des = "filtered." * src)

Filter vcf files (.vcf or .vcf.gz) according to row indices `record_index` 
and column indices `sample_index` and write to a new set of vcf files `des`.

# Input:
- `src`: full vcf file name
- `record_index`: row indices to keep.
- `sample_index`: column indices to keep.

# Optional arguments:
- `des`: output vcf file name; default `"filtered." * src`.
- `allow_multiallelic`: If `false`, multi-allelic SNPs will be filtered out. default `false` 
"""
function filter(
    src::AbstractString, 
    record_index::AbstractVector{<:Integer}, 
    sample_index::AbstractVector{<:Integer}; 
    des::AbstractString = "filtered." * src,
    allow_multiallelic::Bool = false
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
    pmeter = Progress(records, 5, "filtering vcf file...")

    # write to des
    for (i, record) in enumerate(reader)
        if record_mask[i] # keep this record
            if !allow_multiallelic && length(VCF.alt(record)) > 1 
                next!(pmeter)
                continue
            end
            VCFTools.filter_record!(record, sample_mask)
            VCF.write(writer, record)
        end
        next!(pmeter)
    end

    close(reader)
    flush(writer)
    close(writer)
    return nothing
end

"""
    filter_header(reader, sample_mask)

Filters the header info of a VCF reader. Will not keep sample IDs for entry `i` if `sample_mask[i] = false`.

# Inpute:
- `reader`: a VCF reader object.
- `sample_mask`: a BitVector. `sample_mask[i] = true` means keep sample `i`. 

TODO: Create VCFTool.jl signature in meta info. If filedate/signature exist already, delete them before writing new ones
"""
function filter_header(
    reader::VCF.Reader,
    sample_mask::BitVector
    )
    # save meta information
    # fileformat = header(reader).metainfo[1]
    # filedate   = VCF.MetaInfo("##filedate=" * string(Dates.today()))
    # signature  = VCF.MetaInfo("##source=VCFTools.jl")
    # metainfo = vcat(fileformat, filedate, signature, header(reader).metainfo[2:end])

    # filter sampleID
    metainfo = header(reader).metainfo
    sampleID = header(reader).sampleID[sample_mask]
    return VCF.Header(metainfo, sampleID)
end

"""
    filter_record!(record, sample_mask)

Filters a VCF record by samples. 

# Inputs:
- `record`: a VCF record
- `sample_mask`: a BitVector. `sample_mask[i] = true` means keep sample `i`. 
"""
function filter_record!(
    record::VCF.Record,
    sample_mask::BitVector,
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
    return nothing
end

"""
    mask_gt(src, masks; [des = "masked." * src], [separator = '/'])

Creates a new VCF file `des` where genotype entry (i, j) of `src` 
is missing if `masks[i, j]` is true. `src` is unchanged.

# Arguments
- `src`: Input VCF file name.
- `masks`: Bit matrix. `masks[i, j] = true` means mask entry (i, j).
- `des`: output VCF file name.
- `unphase`: If `true`, all heterozygous genotypes will be `1/0` and homozygous separator will be `/`. If `false`, all unmasked genotypes will be retained.  
"""
function mask_gt(
    src::AbstractString, 
    masks::BitArray{2}; 
    des::AbstractString = "masked." * src,
    unphase::Bool = false
    )
    records, samples = nrecords(src), nsamples(src)
    p, n = size(masks)
    if !(records == p && samples == n)
        throw(DimensionMismatch("size(src) = ($records, $samples) ≠ size(masks) = ($p, $n)."))
    end
    pmeter = Progress(records, 5, "masking vcf file...")

    # define byte representation of separator
    separator = (unphase ? 0x2f : 0x7c) # '/' or '|'

    # create input and output VCF files
    reader = VCF.Reader(openvcf(src, "r"))
    writer = VCF.Writer(openvcf(des, "w"), header(reader))
    record = read(reader)

    # loop over each record
    i = 1 # record (row) counter
    while true
        gtkey = findgenokey(record, "GT")
        if gtkey != nothing 
            # loop over genotypes
            for (j, geno) in enumerate(record.genotype)
                if masks[i, j]
                    # change entry to "./." or ".|."
                    record.data[geno[gtkey][1]] = 0x2e 
                    record.data[geno[gtkey][2]] = separator
                    record.data[geno[gtkey][3]] = 0x2e
                end
                if unphase
                    a1 = record.data[geno[gtkey][1]] == 0x31 # "1" (ALT) => 0x31
                    a2 = record.data[geno[gtkey][3]] == 0x31
                    if a1 + a2 == 1
                        record.data[geno[gtkey][1]] = 0x31
                        record.data[geno[gtkey][2]] = separator
                        record.data[geno[gtkey][3]] = 0x30
                    else
                        record.data[geno[gtkey][2]] = separator
                    end
                end
            end
        end

        write(writer, record)

        if eof(reader) 
            next!(pmeter); break
        else
            read!(reader, record)
            i += 1
        end
        next!(pmeter)
    end
    flush(writer); close(reader); close(writer)
    return nothing
end

"""
    find_duplicate_marker(vcffile::String)

Loops over a vcf file and outputs a `BitVector` indicating duplicate records (SNPs)
by checking marker positions. `true` means the marker has been seen before. 

The first occurance will be false and any subsequent occurance will be true. 
"""
function find_duplicate_marker(vcffile::String)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    records = nrecords(vcffile)
    duplicates = falses(records)
    seen = BitSet()
    pmeter = Progress(records, 5, "finding duplicate markers...")
    @inbounds for (i, record) in enumerate(reader)
        curr = VCF.pos(record)
        if curr in seen
            duplicates[i] = true
        else
            push!(seen, curr)
        end
        next!(pmeter)
    end
    close(reader)
    return duplicates
end

"""
    filter_chr(src, chr; des="filtered.(chr)." * src)

Leave only the variants from the chromosome `chr`. Returns the number of records remaining.
"""
function filter_chr(
    src::AbstractString, 
    chr::Union{AbstractString, Int}; 
    des::AbstractString = "filtered.chr$chr." * src
    )

    # create input and output VCF files
    reader = VCF.Reader(openvcf(src, "r"))
    metainfo = header(reader).metainfo
    sampleID = header(reader).sampleID
    writer = VCF.Writer(openvcf(des, "w"), VCF.Header(metainfo, sampleID))
    cnt = 0
    for record in reader
        if VCF.chrom(record) == string(chr)
            cnt += 1
            VCF.write(writer, record)
        end
    end

    close(reader)
    flush(writer)
    close(writer)
    return cnt
end

"""
    filter_range(src, chr, from, to; des="filtered.chr(chr):(from)-(to)." * src)

Leave only the variants from the range `chr:from-to`. Returns the number of records remaining.
"""
function filter_range(
    src::AbstractString,
    chr::Union{AbstractString, Int},
    from::Int, to::Int;
    des::AbstractString = "filtered.chr$chr:$from-$to." * src
    )
    reader = VCF.Reader(openvcf(src, "r"))
    metainfo = header(reader).metainfo
    sampleID = header(reader).sampleID
    writer = VCF.Writer(openvcf(des, "w"), VCF.Header(metainfo, sampleID))
    cnt = 0
    for record in reader
        if VCF.chrom(record) == string(chr)
            if from <= VCF.pos(record) <= to
                VCF.write(writer, record)
                cnt += 1
            end
        end
    end
    close(reader)
    flush(writer)
    close(writer)
    return cnt
end