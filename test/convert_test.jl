using GeneticVariation

@testset "convert_gt(vcfile)" begin
    vcffile = "test.08Jun17.d8b.vcf.gz"
    isfile(vcffile) || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        abspath(joinpath(dirname(pathof(VCFTools)), "..", "test/$vcffile"))) 
    #@code_warntype convert_gt(UInt8, vcffile; impute = false, center = false, scale = false)
    @inferred convert_gt(UInt8, vcffile; impute = false, center = false, scale = false)
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
