#@code_warntype gtstats("test.08Jun17.d8b.vcf.gz", STDOUT)
#@time gtstats("test.08Jun17.d8b.vcf.gz")
#@time gtstats("test.08Jun17.d8b.vcf")
@time records, samples, lines = gtstats("test.08Jun17.d8b.vcf.gz", "gtstats.out.txt")
@show records, samples, lines
