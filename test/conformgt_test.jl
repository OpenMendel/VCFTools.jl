@testset "filter_genotype(record)" begin
    record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
    #@code_warntype filter_genotype(record, ["GT"])
    # GT field
    record_out = filter_genotype(record, ["GT"])
    @test VCF.format(record_out) == ["GT"]
    @test VCF.genotype(record_out) == [["0|0"], ["1|0"]]
    # DP field
    record_out = filter_genotype(record, ["DP"])
    @test VCF.format(record_out) == ["DP"]
    @test VCF.genotype(record_out) == [["1"], ["8"]]
    # GT and DP fields
    record_out = filter_genotype(record, ["GT", "DP"])
    @test VCF.format(record_out) == ["DP", "GT"]
    @test VCF.genotype(record_out) == [["1", "0|0"], ["8", "1|0"]]
    # Only one matching format
    record_out = filter_genotype(record, ["GT", "DS"])
    @test VCF.format(record_out) == ["GT"]
    @test VCF.genotype(record_out) == [["0|0"], ["1|0"]]
    # no matching field
    record_out = filter_genotype(record, ["DS"])
    @test VCF.hasformat(record_out) == false
end

@testset "filter_genotype(file)" begin
    # retrieve (compressed) GT data 
    filter_genotype("test.08Jun17.d8b.vcf.gz")
    reader_in = VCF.Reader(openvcf("test.08Jun17.d8b.vcf.gz", "r"))
    reader_out = VCF.Reader(openvcf("filtered.vcf.gz", "r"))
    iter_state_in, iter_state_out = iterate(reader_in), iterate(reader_out)
    while iter_state_in !== nothing
        record_in, state_in = iter_state_in
        record_out, state_out = iter_state_out
        @test VCF.format(record_out) == ["GT"]
        @test VCF.genotype(record_in, :, "GT") == VCF.genotype(record_out, :, "GT")
        iter_state_in = iterate(reader_in, state_in)
        iter_state_out = iterate(reader_out, state_out)
    end
    close(reader_in); close(reader_out)
    # no matching formats
    filter_genotype("test.08Jun17.d8b.vcf.gz", "filtered.vcf.gz", ["GQ"])
    reader_out = VCF.Reader(openvcf("filtered.vcf.gz", "r"))
    iter_state_out = iterate(reader_out)
    while iter_state_out !== nothing
        record_out, state_out = iter_state_out
        @test VCF.hasformat(record_out) == false
        iter_state_out = iterate(reader_out, state_out)
    end
    close(reader_out)
end

@testset "flip_gt_allele" begin
    record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
    #@code_warntype flip_gt_allele(record)
    record_out = VCFTools.flip_gt_allele(record)
    #@show record_out
    @test VCF.ref(record_out) == "A"
    @test VCF.alt(record_out) == ["G"]
    @test VCF.format(record_out) == ["GT"]
    @test VCF.genotype(record_out) == [["1|1"], ["0|1"]]
end

@testset "match_gt_allele" begin
    record1 = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
    record2 = VCF.Record("20\t14370\trs6054257\tA\tG\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|1:48:1:51,51\t0|1:48:8:51,51")
    #@code_warntype match_gt_allele(record1, record2)
    #@inferred match_gt_allele(record1, record2)
    record1_out, record2_out = VCFTools.match_gt_allele(record1, record2)
    @test VCF.ref(record1_out) == VCF.ref(record2_out)
    @test VCF.alt(record1_out) == VCF.alt(record2_out)
    @test VCF.genotype(record1_out) == VCF.genotype(record2_out)
    # ambiguous case
    record3 = VCF.Record("20\t14370\trs6054257\tC\tT\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|1:48:1:51,51\t0|1:48:8:51,51")
    @test_throws ArgumentError VCFTools.match_gt_allele(record1, record3)
end

@testset "binomial_proportion_test" begin
    # large sample test
    # <http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/binotest.htm>
    pval = VCFTools.binomial_proportion_test(32, 6, 39, 5)
    @test pval ≈ 0.55760 atol=1e-4
    # Fisher exact test
    # Lady tasting tea: <https://en.wikipedia.org/wiki/Fisher%27s_exact_test>
    pval = VCFTools.binomial_proportion_test(1, 9, 11, 3)
    @test pval ≈ 0.002759456 atol=1e-6
end

@testset "conformgt_by_id" begin
    refvcf = "chr22.1kg.phase3.v5a.vcf.gz"    
    tgtvcf = "test.08Jun17.d8b.vcf.gz"
    outvcf = "conformgt.matched"
    isfile(refvcf) || Downloads.download("https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz", 
        abspath(joinpath(dirname(pathof(VCFTools)), "..", "test/$refvcf")))
    isfile(tgtvcf) || Downloads.download("https://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz", 
        abspath(joinpath(dirname(pathof(VCFTools)), "..", "test/$tgtvcf")))
    #@code_warntype conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    #@test @inferred conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    @time lines = conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    @test lines == 823
    reader_ref = VCF.Reader(openvcf(join([outvcf, ".ref.vcf.gz"]), "r"))
    reader_tgt = VCF.Reader(openvcf(join([outvcf, ".tgt.vcf.gz"]), "r"))
    iter_state_ref, iter_state_tgt = iterate(reader_ref), iterate(reader_tgt)
    while iter_state_ref !== nothing
        record_ref, state_ref = iter_state_ref
        record_tgt, state_tgt = iter_state_tgt
        @test VCF.id(record_ref) == VCF.id(record_tgt)
        @test VCF.ref(record_ref) == VCF.ref(record_tgt)
        # @test VCF.alt(record_ref) == VCF.alt(record_tgt)
        iter_state_ref = iterate(reader_ref, state_ref)
        iter_state_tgt = iterate(reader_tgt, state_tgt)
    end
    @time lines = conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, 0.05)
    @test lines == 488
    @time lines = conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, 1)
    @test lines == 0
    # Profile.clear()
    # @profile conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    # Profile.print(format=:flat)
end

@testset "conformgt_by_pos" begin
    refvcf = "chr22.1kg.phase3.v5a.vcf.gz"    
    tgtvcf = "test.08Jun17.d8b.vcf.gz"
    outvcf = "conformgt.matched"
    isfile(refvcf) || Downloads.download("https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz", 
        abspath(joinpath(dirname(pathof(VCFTools)), "..", "test/$refvcf"))) 
    isfile(tgtvcf) || Downloads.download("https://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz", 
        abspath(joinpath(dirname(pathof(VCFTools)), "..", "test/$tgtvcf")))
    #@code_warntype conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    #@test @inferred conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    @time lines = conformgt_by_pos(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    @test lines == 833
    reader_ref = VCF.Reader(openvcf(join([outvcf, ".ref.vcf.gz"]), "r"))
    reader_tgt = VCF.Reader(openvcf(join([outvcf, ".tgt.vcf.gz"]), "r"))
    iter_state_ref, iter_state_tgt = iterate(reader_ref), iterate(reader_tgt)
    while iter_state_ref !== nothing
        record_ref, state_ref = iter_state_ref
        record_tgt, state_tgt = iter_state_tgt
        @test VCF.chrom(record_ref) == VCF.chrom(record_tgt)
        @test VCF.pos(record_ref) == VCF.pos(record_tgt)
        @test VCF.ref(record_ref) == VCF.ref(record_tgt)
        # @test VCF.alt(record_ref) == VCF.alt(record_tgt)
        iter_state_ref = iterate(reader_ref, state_ref)
        iter_state_tgt = iterate(reader_tgt, state_tgt)
    end
    @time lines = conformgt_by_pos(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, 0.05)
    @test lines == 493
    @time lines = conformgt_by_pos(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, 1)
    @test lines == 0
    # Profile.clear()
    # @profile conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    # Profile.print(format=:flat)
end
