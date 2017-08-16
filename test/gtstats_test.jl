using GeneticVariation

info("Hardy-Weinberg equilibrium test")
# <https://en.wikipedia.org/wiki/Hardy–Weinberg_principle#Example_.7F.27.22.60UNIQ--postMath-0000001E-QINU.60.22.27.7F_test_for_deviation>
#@code_warntype hwe(1469, 138, 5)
@test hwe(1469, 138, 5) ≈ 0.406002899

info("gtstat(record)")
record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
#@code_warntype gtstats(record)
@test gtstats(record)[1:end-1] == (1, 1, 0, 3, 1, 0.75, 0.25, 0, 1, 0.25)
@test gtstats(record)[end] ≈ 0.82414089574

info("gtstats(vcfile, ofile)")
#@code_warntype gtstats("test.08Jun17.d8b.vcf.gz", STDOUT)
#@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz")
#@time gtstats("test.08Jun17.d8b.vcf")
@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.txt")
#@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.gz")
#@time records, samples, lines = gtstats("chr22.1kg.ref.phase1_release_v3.20101123.vcf.gz", "gtstats.out.gz") # about 180 seconds
@show records, samples, lines
