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
        record_index[i] && VCF.write(writer, filter_record(record, sample_mask))
    end

    close(reader)
    flush(writer)
    close(VCF.BioCore.IO.stream(writer))
end

function filter_header(
    reader::VCF.Reader,
    sample_mask::BitVector
    )
    #save meta information
    fileformat = VCF.header(reader).metainfo[1]
    filedate   = VCF.MetaInfo("##filedate=" * string(Dates.today()))
    signature  = VCF.MetaInfo("##source=VCFTools.jl")
    metainfo = vcat(fileformat, filedate, signature, VCF.header(reader).metainfo[2:end])

    #filter sampleID
    sampleID = VCF.header(reader).sampleID[sample_mask]
    return VCF.Header(metainfo, sampleID)
end

"""
    filter_record(record, sample_mask)

TODO: make this efficient
"""
function filter_record(
    record::VCF.Record,
    sample_mask::BitVector
    )
    # filter only genotype data
    new_record = copy(record)
    new_record.genotype = record.genotype[sample_mask]
    return new_record
end
