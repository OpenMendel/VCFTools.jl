function filter(
	src::AbstractString, 
	record_index::AbstractVector{<:Integer}, 
	sample_index::AbstractVector{<:Integer}; 
	des::AbstractString = "filtered." * src 
	)

    # create record (row) and sample (column) masks
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
    writer = VCF.Writer(openvcf(des, "w"), filtered_header(reader, sample_mask))
    
    # needed constants
    records = nrecords(src)
    samples = nsamples(src)

    # write to des
    # lines = 0
    # for record in reader
    # 	lines += 1
    # 	VCF.write(writer, record)
    # end

    close(reader)
    flush(writer)
    close(VCF.BioCore.IO.stream(writer))
end

function filtered_header(
	reader::VCF.Reader,
	sample_mask::BitVector
	)
	metainfo = VCF.header(reader).metainfo
	sampleID = VCF.header(reader).sampleID[sample_mask]
	return VCF.Header(metainfo, sampleID)
end

