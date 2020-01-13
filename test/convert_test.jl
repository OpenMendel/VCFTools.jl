using GeneticVariation

@testset "convert_gt" begin
    #@code_warntype convert_gt(Float64, (false, false), true, :additive)
    @inferred convert_gt(Float64, (false, false), true, :additive)
    # REF is the minor allele
    @test convert_gt(Float64, (false, false), true, :additive) == 0.0
    @test convert_gt(Float64, (false, true), true, :additive) == 1.0
    @test convert_gt(Float64, (true, false), true, :additive) == 1.0
    @test convert_gt(Float64, (true, true), true, :additive) == 2.0
    @test convert_gt(Float64, (false, false), true, :dominant) == 0.0
    @test convert_gt(Float64, (false, true), true, :dominant) == 1.0
    @test convert_gt(Float64, (true, false), true, :dominant) == 1.0
    @test convert_gt(Float64, (true, true), true, :dominant) == 1.0
    @test convert_gt(Float64, (false, false), true, :recessive) == 0.0
    @test convert_gt(Float64, (false, true), true, :recessive) == 0.0
    @test convert_gt(Float64, (true, false), true, :recessive) == 0.0
    @test convert_gt(Float64, (true, true), true, :recessive) == 1.0
    # ALT is the minor allele
    @test convert_gt(Float64, (false, false), false, :additive) == 2.0
    @test convert_gt(Float64, (false, true), false, :additive) == 1.0
    @test convert_gt(Float64, (true, false), false, :additive) == 1.0
    @test convert_gt(Float64, (true, true), false, :additive) == 0.0
    @test convert_gt(Float64, (false, false), false, :dominant) == 1.0
    @test convert_gt(Float64, (false, true), false, :dominant) == 1.0
    @test convert_gt(Float64, (true, false), false, :dominant) == 1.0
    @test convert_gt(Float64, (true, true), false, :dominant) == 0.0
    @test convert_gt(Float64, (false, false), false, :recessive) == 1.0
    @test convert_gt(Float64, (true, false), false, :recessive) == 0.0
    @test convert_gt(Float64, (false, true), false, :recessive) == 0.0
    @test convert_gt(Float64, (true, true), false, :recessive) == 0.0
end

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
    # convert to a matrix of Float32 without checking minor alleles (reading 0/1 as-is)
    # @code_warntype convert_gt(Float32, vcffile; as_minorallele = false)
    @time C = convert_gt(Float32, vcffile; as_minorallele = false)
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
    # convert next 5 records to a Float32 matrix without checking which is minor alleles
    C5 = Matrix{Union{Float64, Missing}}(undef, nsamples(vcffile), 5)
    @time copy_gt_as_is!(C5, reader)
    @test findall(x -> x === missing, C5) == []
    close(reader)
end

@testest "convert_ht" begin
    # REF is the minor allele
    @test convert_ht(Float64, true, true) == 1.0
    @test convert_ht(Float64, false, true) == 0.0
    # ALT is the minor allele
    @test convert_ht(Float64, true, false) == 0.0
    @test convert_ht(Float64, false, false) == 1.0
end

@testset "convert_ht(vcfile)" begin
    vcffile = "test.08Jun17.d8b.vcf.gz"
    # convert to a reference haplotype panel of Float32, no checking for minor allele
    @time H1 = convert_ht(Float32, vcffile, as_minorallele = false)
    @test size(H1) == (382, 1356)
    @test typeof(H1) == Matrix{Float32}
    # convert to a reference haplotype panel of Float64, checking for minor allele
    @time H2 = convert_ht(Float64, vcffile, as_minorallele = true)
    @test size(H2) == (382, 1356)
    @test typeof(H2) == Matrix{Float64}
    # convert first record into 2 haplotype vectors and check their sum is the genotype vector
    reader = VCF.Reader(openvcf(vcffile, "r"))
    h1h2 = Matrix{UInt8}(undef, 2, nrecords(vcffile))
    # g1   = Matrix{Union{UInt8, Missing}}(undef, 1, nrecords(vcffile))
    # copy_ht_as_is!(h1h2, reader)
    # copy_gt_as_is!(g1, reader)
end
