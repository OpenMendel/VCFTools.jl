@testset "filter by sample/record index" begin
    # filter samples only specified via integer indices 
    vcffile = "test.08Jun17.d8b.vcf.gz"
    samples = nsamples(vcffile)
    records = nrecords(vcffile)
    record_mask = collect(1:records)
    sample_mask = collect(2:(samples - 1))
    des = "filtered.test.08Jun17.d8b.vcf.gz"
    @time VCFTools.filter(vcffile, record_mask, sample_mask, des=des)
    # @code_warntype VCFTools.filter(vcffile, record_mask, sample_mask, des=des)
    @test nsamples(des) == samples - 2
    @test nrecords(des) == records
    X = convert_gt(Float32, vcffile, as_minorallele=false)
    X_filter = convert_gt(Float32, des, as_minorallele=false)
    @test all(X[sample_mask, :] .== X_filter)

    # filter records only specified by bitmasks
    record_mask = trues(records)
    sample_mask = trues(samples)
    record_mask[1] = record_mask[end] = false
    des = "filtered.test.08Jun17.d8b.vcf.gz"
    @time VCFTools.filter(vcffile, record_mask, sample_mask, des=des)
    @test nsamples(des) == samples
    @test nrecords(des) == records - 2
    X = convert_gt(Float32, vcffile, as_minorallele=false)
    X_filter = convert_gt(Float32, des, as_minorallele=false)
    @test all(X[:, record_mask] .== X_filter)

    # filter record and samples together
    record_mask = trues(records)
    sample_mask = trues(samples)
    record_mask[1] = record_mask[end] = false
    sample_mask[1] = sample_mask[end] = false
    des = "filtered.test.08Jun17.d8b.vcf.gz"
    @time VCFTools.filter(vcffile, record_mask, sample_mask, des=des)
    @test nsamples(des) == sum(sample_mask)
    @test nrecords(des) == sum(record_mask)
    X = convert_gt(Float32, vcffile, as_minorallele=false)
    X_filter = convert_gt(Float32, des, as_minorallele=false)
    @test all(X[sample_mask, record_mask] .== X_filter)
end

@testset "mask_gt" begin
    vcffile = "test.08Jun17.d8b.vcf.gz"
    outfile = "masked.test.08Jun17.d8b.vcf.gz"
    samples = nsamples(vcffile)
    records = nrecords(vcffile)
    masks   = falses(records, samples)
    [masks[i, i] = true for i in 1:10]
    @time mask_gt(vcffile, masks, des = outfile)
    # @code_warntype mask_gt(vcffile, masks, des = outfile)
    @test nsamples(outfile) == samples
    @test nrecords(outfile) == records
    X = convert_gt(Float32, vcffile, as_minorallele=false)
    X = copy(X') # need to transpose back to dimension in VCF file
    Xmasked = convert_gt(Float32, outfile, as_minorallele=false)
    Xmasked = copy(Xmasked') 
    @test eltype(Xmasked) == Union{Missing, eltype(X)}
    @test all(Xmasked[.!masks] .== X[.!masks])
    @test all(Xmasked[masks] .=== missing)
end
