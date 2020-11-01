@testset "grm" begin
    vcf = "test.08Jun17.d8b.vcf.gz"
    isfile(vcf) || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        abspath(joinpath(dirname(pathof(VCFTools)), "..", "test/$vcf"))) 
    # standard GRM method
    Φgrm = grm(vcf, method=:GRM)
    n = nsamples(vcf)
    p = nrecords(vcf)
    @test size(Φgrm) == (n, n)
    @test issymmetric(Φgrm)
    @test eigmin(Φgrm) > -1e-8
    # MoM estimates
    Φmom = grm(vcf, method=:MoM)
    @test size(Φmom) == (n, n)
    @test issymmetric(Φmom)
    # Robust GRM
    Φrbs = grm(vcf, method=:Robust)
    @test size(Φrbs) == (n, n)
    @test issymmetric(Φrbs)
    @test eigmin(Φrbs) > -1e-8
    # test scale_missing
    Random.seed!(2020)
    masked_vcf = "masked.test.08Jun17.d8b.vcf.gz"
    mask_gt(vcf, bitrand(p, n), des=masked_vcf) # create vcf file with missings
    Φrbs_sm = grm(masked_vcf, method=:Robust, scale_missing=true)
    @test issymmetric(Φrbs_sm)
    @test eigmin(Φrbs_sm) > -1e-8
    @test size(Φrbs_sm) == (n, n)
end