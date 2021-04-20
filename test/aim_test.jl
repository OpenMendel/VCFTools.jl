@testset "AIM selection" begin
    # download data
    vcffile = "chr22.1kg.phase3.v5a.vcf.gz"
    isfile(vcffile) || 
        download("http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz", 
        joinpath(pwd(), vcffile))

    # cd to test folder and read population origin
    joinpath(pathof(VCFTools), "test")
    df = CSV.read("1000genomes.population.txt", DataFrame)

    # create dictionary with key = ID, value = population 
    sampleID_to_population = Dict{String, String}()
    for (id, population) in eachrow(df)
        sampleID_to_population[id] = population
    end
    sampleID_to_population

    pvals = VCFTools.aim_select(vcffile, sampleID_to_population)
    @test all(0.0 .≤ pvals .≤ 1.0)
    @test length(pvals) == 424147
    @test pvals[1] == 1.0
    @test pvals[11] ≈ 5.684678238417218e-120
end
