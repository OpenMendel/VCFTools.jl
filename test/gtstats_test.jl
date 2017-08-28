using GeneticVariation, CodecZlib

@testset "H-W test" begin
    # <https://en.wikipedia.org/wiki/Hardy–Weinberg_principle#Example_.7F.27.22.60UNIQ--postMath-0000001E-QINU.60.22.27.7F_test_for_deviation>
    #@code_warntype hwe(1469, 138, 5)
    @test VCFTools.hwe(1469, 138, 5) ≈ 0.406002899
end

@testset "openvcf" begin
    # download test file if not exist
    isfile("test.08Jun17.d8b.vcf.gz") || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(Pkg.dir("VCFTools"), "test/test.08Jun17.d8b.vcf.gz"))
    @test_throws ArgumentError openvcf("test.08Jun17.d8b.vcf.gz", "r+")
    @test_throws ArgumentError openvcf("test.08Jun17.d8b.vcf.gz", "w+")
    @test_throws ArgumentError openvcf("test.08Jun17.d8b.vcf.gz", "c")
    # VCF file name has to end with vcf or vcf.gz
    @test_throws ArgumentError openvcf("test.08Jun17.d8b")
end


@testset "nrecords(vcf)" begin
    # download test file if not exist
    isfile("test.08Jun17.d8b.vcf.gz") || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(Pkg.dir("VCFTools"), "test/test.08Jun17.d8b.vcf.gz"))
    @test nrecords("test.08Jun17.d8b.vcf.gz") == 1356
end

@testset "nsamples(vcf)" begin
    # download test file if not exist
    isfile("test.08Jun17.d8b.vcf.gz") || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(Pkg.dir("VCFTools"), "test/test.08Jun17.d8b.vcf.gz"))
    @test nsamples("test.08Jun17.d8b.vcf.gz") == 191
end

@testset "gtstats(record)" begin
    record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
    #@code_warntype gtstats(record)
    @test gtstats(record)[1:end-1] == (1, 1, 0, 3, 1, 0.75, 0.25, 0, true, 0.25)
    @test gtstats(record)[end] ≈ 0.82414089574
    # with missing data
    record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t.\t1|0:48:8:51,51")
    @test gtstats(record)[1:end-1] == (0, 1, 0, 1, 1, 0.5, 0.5, 1, false, 0.5)
    @test gtstats(record)[end] ≈ 0.31731050786291415
end

@testset "gtstats(vcf)" begin
    #@code_warntype gtstats("test.08Jun17.d8b.vcf.gz", STDOUT)
    #@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz")
    #@time gtstats("test.08Jun17.d8b.vcf")
    # output tuple: (records, samples, lines, missings_by_sample, missings_by_record, maf_by_record, minorallele_by_record)

    # download test file and unzip
    isfile("test.08Jun17.d8b.vcf.gz") || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(Pkg.dir("VCFTools"), "test/test.08Jun17.d8b.vcf.gz"))
    write("test.08Jun17.d8b.vcf", read(GzipDecompressionStream(open("test.08Jun17.d8b.vcf.gz", "r"))))

    @testset "input: text file, output: text file" begin
        @time out = gtstats("test.08Jun17.d8b.vcf", "gtstats.out.txt")
        @test out[1:3] == (1356, 191, 1356)
        outtxt = readdlm("gtstats.out.txt")
        @test all(outtxt[:, 9]  .== out[5]) # missings_by_record
        @test all(outtxt[:, 14] .== out[6]) # maf_by_record
    end

    @testset "input: text file, output: gz file" begin
        @time out = gtstats("test.08Jun17.d8b.vcf", "gtstats.out.gz")
        @test out[1:3] == (1356, 191, 1356)
        outtxt = readdlm(GzipDecompressionStream(open("gtstats.out.gz", "r")))
        @test all(outtxt[:, 9]  .== out[5]) # missings_by_record
        @test all(outtxt[:, 14] .== out[6]) # maf_by_record
    end

    @testset "input: gz file, output: text file" begin
        @time out = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.txt")
        @test out[1:3] == (1356, 191, 1356)
        outtxt = readdlm("gtstats.out.txt")
        @test all(outtxt[:, 9]  .== out[5]) # missings_by_record
        @test all(outtxt[:, 14] .== out[6]) # maf_by_record
    end

    @testset "input: gz file, output: gz file" begin
        @time out = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.gz")
        @test out[1:3] == (1356, 191, 1356)
        outtxt = readdlm(GzipDecompressionStream(open("gtstats.out.gz", "r")))
        @test all(outtxt[:, 9]  .== out[5]) # missings_by_record
        @test all(outtxt[:, 14] .== out[6]) # maf_by_record
    end

    #@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.gz")
    #@time records, samples, lines = gtstats("chr22.1kg.ref.phase1_release_v3.20101123.vcf.gz", "gtstats.out.gz") # about 180 seconds
    #@show maf_by_record
end
