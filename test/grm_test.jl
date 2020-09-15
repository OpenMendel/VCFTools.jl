using Revise
using VCFTools
using Test
using LinearAlgebra

@testset "grm" begin
    cd("/Users/biona001/Benjamin_Folder/UCLA/research/BIG summer/2020 summer/data")
    vcf = "odis_ddrad.merged.het_filtered.hard_filtered.arthropod.minDP3.vcf"
    # standard GRM method
    Φgrm = grm(vcf, method=:GRM)
    n = nsamples(vcf)
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

    # difference
    Φrbs - Φmom

    [Φgrm[:, 1] Φrbs[:, 1] Φmom[:, 1]]

end