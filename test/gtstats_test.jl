@testset "H-W test" begin
    # <https://en.wikipedia.org/wiki/Hardy–Weinberg_principle#Example_.7F.27.22.60UNIQ--postMath-0000001E-QINU.60.22.27.7F_test_for_deviation>
    #@code_warntype hwe(1469, 138, 5)
    @test VCFTools.hwe(1469, 138, 5) ≈ 0.406002899
end

@testset "openvcf" begin
    # download and extract test file if not exist
    isfile("test.08Jun17.d8b.vcf.gz") || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz"))
    @test_throws ArgumentError openvcf("test.08Jun17.d8b.vcf.gz", "r+")
    @test_throws ArgumentError openvcf("test.08Jun17.d8b.vcf.gz", "w+")
    @test_throws ArgumentError openvcf("test.08Jun17.d8b.vcf.gz", "c")
    # VCF file name has to end with vcf or vcf.gz
    @test_throws ArgumentError openvcf("test.08Jun17.d8b")
end

@testset "nrecords(vcf)" begin
    # download and extract test file if not exist
    isfile("test.08Jun17.d8b.vcf.gz") ||   download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz"))
    @test nrecords("test.08Jun17.d8b.vcf.gz") == 1356
end

@testset "nsamples(vcf)" begin
    # download test file if not exist
    isfile("test.08Jun17.d8b.vcf.gz") ||  download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz"))
    @test nsamples("test.08Jun17.d8b.vcf.gz") == 191
end

@testset "gtstats(record)" begin
    record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
    #@code_warntype gtstats(record)
    @test gtstats(record)[1:end-1] == (1, 1, 0, 3, 1, 0.25, 0.75, 0, false, 0.25)
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

    # download and extract test file and unzip
    isfile("test.08Jun17.d8b.vcf.gz") || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz"))
    write("test.08Jun17.d8b.vcf", read(GzipDecompressorStream(open("test.08Jun17.d8b.vcf.gz", "r"))))

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
        outtxt = readdlm(GzipDecompressorStream(open("gtstats.out.gz", "r")))
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
        outtxt = readdlm(GzipDecompressorStream(open("gtstats.out.gz", "r")))
        @test all(outtxt[:, 9]  .== out[5]) # missings_by_record
        @test all(outtxt[:, 14] .== out[6]) # maf_by_record
    end

    #@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.gz")
    #@time records, samples, lines = gtstats("chr22.1kg.ref.phase1_release_v3.20101123.vcf.gz", "gtstats.out.gz") # about 180 seconds
    #@show maf_by_record
end

@testset "Fishers exact HWE test" begin
    # See p20-21 in http://courses.washington.edu/b516/lectures_2009/HWE_Lecture.pdf
    n0 = 179
    n1 = 21
    N = (n0 + n1) >> 1
    n01 = 1
    @test VCFTools.fisher_exact(n01, N, n0) < .000001
    n01 = 3
    @test VCFTools.fisher_exact(n01, N, n0) < .000001
    n01 = 5
    @test VCFTools.fisher_exact(n01, N, n0) < .000001
    n01 = 7
    @test round(VCFTools.fisher_exact(n01, N, n0), digits=6) == .000001
    n01 = 9
    @test round(VCFTools.fisher_exact(n01, N, n0), digits=6) == .000047
    n01 = 11
    @test round(VCFTools.fisher_exact(n01, N, n0), digits=6) == .000870
    n01 = 13
    @test round(VCFTools.fisher_exact(n01, N, n0), digits=6) == .009375
    n01 = 15
    @test round(VCFTools.fisher_exact(n01, N, n0), digits=6) == .059283
    n01 = 17
    @test round(VCFTools.fisher_exact(n01, N, n0), digits=6) == .214465
    n01 = 19
    @test round(VCFTools.fisher_exact(n01, N, n0), digits=6) == .406355
    n01 = 21
    @test round(VCFTools.fisher_exact(n01, N, n0), digits=6) == .309604

    # Picked a few values in Table 1 of 
    # AJHG Volume 76, Issue 5, May 2005, Pages 887-893
    # https://www.sciencedirect.com/science/article/pii/S0002929707607356
    N  = 100
    n0 = 21
    o01 = 11
    hwepval, ts_low, ts_high = VCFTools.hwe_fisher(o01, N, n0)
    @test round(hwepval, digits=6) == .000919
    @test round(ts_low, digits=6)  == .000919
    @test round(ts_high, digits=6) == .999952
    o01 = 17
    hwepval, ts_low, ts_high = VCFTools.hwe_fisher(o01, N, n0)
    @test round(hwepval, digits=6) == .284042
    @test round(ts_low, digits=6)  == .284042
    @test round(ts_high, digits=6) == .930424
    o01 = 21
    hwepval, ts_low, ts_high = VCFTools.hwe_fisher(o01, N, n0)
    @test round(hwepval, digits=6) == .593645
    @test round(ts_low, digits=6)  == 1.0
    @test round(ts_high, digits=6) == .309604
end