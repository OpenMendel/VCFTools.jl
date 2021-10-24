@testset "convert_gt(vcfile)" begin
    vcffile = "test.08Jun17.d8b.vcf.gz"
    isfile(vcffile) || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        abspath(joinpath(dirname(pathof(VCFTools)), "..", "test/$vcffile"))) 
    #@code_warntype convert_gt(UInt8, vcffile; impute = false, center = false, scale = false)
    # @inferred convert_gt(UInt8, vcffile; impute = false, center = false, scale = false)
    # convert to a matrix of UInt8, no impute/center/scale (default)
    @time A = convert_gt(UInt8, vcffile)
    @test size(A) == (191, 1356)
    @test eltype(A) == Union{UInt8, Missing}
    # convert to a matrix of Float64, impute = center = scale = true
    @time B = convert_gt(Float64, vcffile; impute = true, center = true, scale = true)
    @test size(B) == (191, 1356)
    @test eltype(B) == Union{Float64, Missing}
    # convert to a matrix of Float32
    # @code_warntype convert_gt(Float32, vcffile)
    @time C = convert_gt(Float32, vcffile)
    @test size(C) == (191, 1356)
    @test eltype(C) == Union{Float32, Missing}
    # read the first record to a UInt8 vector, no impute/center/scale (default)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    A1 = Vector{Union{UInt8, Missing}}(undef, nsamples(vcffile))
    @time copy_gt!(A1, reader)
    @test all(A1[findall(!ismissing, A1)] .== A[findall(!ismissing, A[:, 1])])
    @test all(findall(ismissing, A1) .== findall(ismissing, A[:, 1]))
    # convert next 5 records to a Float64 matrix, impute = center = scale = true
    B5 = Matrix{Union{Float64, Missing}}(undef, nsamples(vcffile), 5)
    @time copy_gt!(B5, reader; impute = true, center = true, scale = true)
    B_copy = copy(B[:, 2:6])
    @test all(B5[findall(!ismissing, B5)] .== B_copy[findall(!ismissing, B[:, 2:6])])
    @test all(B5[findall(ismissing, B5)] .== B_copy[findall(ismissing, B_copy)])
    # convert next 5 records to a Float32 matrix
    C5 = Matrix{Union{Float64, Missing}}(undef, nsamples(vcffile), 5)
    @time copy_gt!(C5, reader)
    @test findall(x -> x === missing, C5) == []
    close(reader)
end

@testset "convert_ht(vcfile)" begin
    vcffile = "test.08Jun17.d8b.vcf.gz"
    # convert to a reference haplotype panel of Float32
    @time H1 = convert_ht(Float32, vcffile)
    @test size(H1) == (382, 1356)
    @test typeof(H1) == Matrix{Float32}
    # convert to a reference haplotype panel of Float64
    @time H2 = convert_ht(Float64, vcffile)
    @test size(H2) == (382, 1356)
    @test typeof(H2) == Matrix{Float64}
    # convert to bitarray
    Hb = convert_ht(Bool, vcffile)
    @test all(H1 .== convert(Matrix{Float32}, Hb))
    # convert first record into 2 haplotype vectors and check their sum is the genotype vector
    reader = VCF.Reader(openvcf(vcffile, "r"))
    h1h2 = Matrix{Float64}(undef, 2, nrecords(vcffile))
    copy_ht!(h1h2, reader)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    g1 = Matrix{Union{Float64, Missing}}(undef, 1, nrecords(vcffile))
    copy_gt!(g1, reader)
    @test all(sum(h1h2, dims=1) .== g1)
end

@testset "convert_ds" begin
    vcffile = "test.08Jun17.d8b.vcf.gz"
    D = convert_ds(Float64, vcffile, key = "DS", impute=false, center=false, scale=false)
    @test D[1, 5] == 1.0
    @test D[1, 1345] == 0.05
    @test D[2, 13] == 1.75
    @test eltype(D) == Union{Missing, Float64}
    @test size(D) == (191, 1356)
end

@testset "convert to transpose functions" begin
    vcffile = "test.08Jun17.d8b.vcf"

    # convert_gt_trans!
    @time A  = convert_gt(Float64, vcffile)
    @time At = convert_gt(Float64, vcffile, trans=true)
    @test all(At .== A')

    # test if eof(reader) is working
    out = Matrix{Union{Float64, Missing}}(undef, 191, 1400)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_gt!(out, reader)
    @test all(ismissing.(out[:, 1357:end]))

    out = Matrix{Union{Float64, Missing}}(undef, 1400, 191)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_gt_trans!(out, reader)
    @test all(ismissing.(out[1357:end, :]))

    # test impute, center, scale (these use mean and var functions, currently not imported)
    # A = Matrix{Union{Float64, Missing}}(undef, 191, 1400)
    # reader = VCF.Reader(openvcf(vcffile, "r"))
    # copy_gt!(A, reader, impute=true, center=true, scale=true)
    # @test isapprox(mean(A[:, 5]), 0, atol=10)
    # @test isapprox(var(A[:, 5]), 1, atol=10)
    # At = Matrix{Union{Float64, Missing}}(undef, 1400, 191)
    # reader = VCF.Reader(openvcf(vcffile, "r"))
    # copy_gt_trans!(At, reader, impute=true, center=true, scale=true)
    # @test isapprox(mean(At[5, :]), 0, atol=10)
    # @test isapprox(var(At[5, :]), 1, atol=10)

    # convert_ht_trans!
    @time H  = convert_ht(Float64, vcffile)
    @time Ht = convert_ht(Float64, vcffile, trans=true)
    @test all(Ht .== H')

    # test if eof(reader) is working
    out = Matrix{Float64}(undef, 382, 1400)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_ht!(out, reader)
    @test size(out) == (382, 1400)
    @test all(out[:, 1358:end] .== 0.0)

    out = Matrix{Float64}(undef, 1400, 382)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_ht_trans!(out, reader)
    @test size(out) == (1400, 382)
    @test all(out[1358:end, :] .== 0.0)

    # convert_ds_trans!
    @time D  = convert_ds(Float64, vcffile)
    @time Dt = convert_ds(Float64, vcffile, trans=true)
    @test all(Dt .== D')

    # test if eof(reader) is working
    out = Matrix{Union{Missing, Float64}}(undef, 191, 1400)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_ds!(out, reader)
    @test size(out) == (191, 1400)
    @test all(ismissing.(out[:, 1358:end]))

    out = Matrix{Union{Missing, Float64}}(undef, 1400, 191)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    copy_ds_trans!(out, reader)
    @test size(out) == (1400, 191)
    @test all(ismissing.(out[1358:end, :]))
end

@testset "record info" begin
    vcffile = "test.08Jun17.d8b.vcf"
    X, sampleID, chr, pos, SNPid, ref, alt = convert_gt(Float64, vcffile, 
        trans=false, save_snp_info=true, msg = "Importing genotype file...")
    Xt, sampleIDt, chrt, post, SNPidt, reft, altt = convert_gt(Float64, vcffile, 
        trans=true, save_snp_info=true, msg = "Importing genotype file...")

    @test all(X .== Xt')
    @test all(sampleID .== sampleIDt)
    @test all(chr .== chrt)
    @test all(pos .== post)
    @test all(SNPid .== SNPidt)
    @test all(ref .== reft)
    @test all(alt .== altt)

    D, sampleID, chr, pos, SNPid, ref, alt = convert_ds(Float64, vcffile, 
        trans=false, save_snp_info=true, msg = "Importing genotype file...")
    Dt, sampleIDt, chrt, post, SNPidt, reft, altt = convert_ds(Float64, vcffile, 
        trans=true, save_snp_info=true, msg = "Importing genotype file...")

    @test all(D .== Dt')
    @test all(sampleID .== sampleIDt)
    @test all(chr .== chrt)
    @test all(pos .== post)
    @test all(SNPid .== SNPidt)
    @test all(ref .== reft)
    @test all(alt .== altt)
end

@testset "write VCF" begin
    # write haplotypes
    H1 = bitrand(100, 200)
    H2 = bitrand(100, 200)
    write_vcf("test.write.vcf.gz", H1, H2)
    H = convert_ht(Bool, "test.write.vcf.gz")
    @test all(H1 .== view(H, 1:2:size(H, 1), :))
    @test all(H2 .== view(H, 2:2:size(H, 1), :))

    # write genotypes
    X = H1 + H2
    write_vcf("test.write.vcf.gz", X)
    G = convert_gt(Float64, "test.write.vcf.gz")
    @test all(G .== X)
end
