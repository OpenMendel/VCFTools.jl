
# Write

Sometimes one wishes to save a numeric matrix to VCF file. This can be accomplished by the [write_vcf](https://openmendel.github.io/VCFTools.jl/dev/man/api/#VCFTools.write_vcf) function.


```julia
using VCFTools
using Random
using Test
```

## Phased genotypes

Let us simulate haplotypes `H1` and `H2`, where each row is a haplotype and each column is a SNP. So Sample 1's gegnotype is `H1[1, :] + H2[1, :]`. We will save these haplotypes into `0|0, 0|1, 1|0, 1|1` accordingly. 


```julia
H1 = bitrand(100, 200) # simulated haplotype 1
H2 = bitrand(100, 200) # simulated haplotype 2
write_vcf("test.write.vcf.gz", H1, H2) # write routine
```

The writen file can of course be re-imported, and we can check that they are the same as H1 and H2:


```julia
H = convert_ht(Bool, "test.write.vcf.gz")
@test all(H1 .== view(H, 1:2:size(H, 1), :))
@test all(H2 .== view(H, 2:2:size(H, 1), :))
```




    [32m[1mTest Passed[22m[39m



## Unphased genotypes

Similarly, `write_vcf` can also save unphased genotypes stored in a numeric matrix `X`. In this case, all heterozygous genotypes will be saved as `1/0`. 


```julia
X = H1 + H2
write_vcf("test.write.vcf.gz", X)
```

We can also check that the written VCF file is the same as the original after re-importing it


```julia
G = convert_gt(Float64, "test.write.vcf.gz")
@test all(G .== X)
```




    [32m[1mTest Passed[22m[39m




```julia
# clean up
rm("test.write.vcf.gz", force=true)
```
