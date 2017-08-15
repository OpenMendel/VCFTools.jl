#@code_warntype gtstats("test.08Jun17.d8b.vcf.gz", STDOUT)
#@time gtstats("test.08Jun17.d8b.vcf.gz")
#@time gtstats("test.08Jun17.d8b.vcf")
#@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.txt")
@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.gz")
#@time records, samples, lines = gtstats("chr22.1kg.ref.phase1_release_v3.20101123.vcf.gz", "gtstats.out.gz") # about 180 seconds
@show records, samples, lines
