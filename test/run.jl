#trying out VCFTools.jl and GeneticVariation.jl
using Revise
using GeneticVariation
using Random
using VCFTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
reader = VCF.Reader(openvcf(vcffile, "r"))
record = read(reader)

record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
record = VCF.Record("20\t14370\t.\tG\tA\t29\tPASS;fdsa\t.\tDS\t1.99\t1.01")
record = VCF.Record("20\t14370\t.\tG\tA\t29\tPASS\t.\tGT\t1\\0\t0\\0")
record.genotype

#test convert_ht
using Revise
using GeneticVariation
using VCFTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
H = convert_ht(Float32, vcffile)

# byte representation mapping
String([0x09]) # '\t'
String([0x3a]) # ':'
String([0x30]) # '0'
String([0x31]) # '1'
String([0x2e]) # '.'
String([0x7c]) # '|'
String([0x2f]) # '/'


# test filter
using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
samples = nsamples(vcffile)
records = nrecords(vcffile)

vcffile = "test.08Jun17.d8b.vcf"
reader = VCF.Reader(openvcf(vcffile, "r"))
record = read(reader)


Random.seed!(123)
record_index = bitrand(records)
# record_index = trues(records)
sample_index = trues(samples)
sample_index[1] = false
sample_index[end] = false

VCFTools.filter(vcffile, record_index, sample_index)
X = convert_gt(Float32, vcffile)
X_filter = convert_gt(Float32, "filtered." * vcffile)

@benchmark VCFTools.filter(vcffile, record_index, sample_index)
#new: 149.763 ms, 214.48 MiB
#old: 134.655 ms, 177.42 MiB

des = "filter." * vcffile
reader = VCF.Reader(openvcf(vcffile, "r"))
writer = VCF.Writer(openvcf(des, "w"), filter_header(reader, sample_index))
    
reader = VCF.Reader(openvcf(vcffile, "r"))
record = read(reader)


# test masking

using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
samples = nsamples(vcffile)
records = nrecords(vcffile)
masks = bitrand(records, samples)
mask_gt(vcffile, masks)
@benchmark mask_gt(vcffile, masks) #82.221 ms, 108.69 MiB, 800187 alloc


# test convert_ds
using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools


cd("/Users/biona001/Benjamin_Folder/UCLA/research/2nd_project/benchmarks/AFRPed/minimac4")
vcffile = "minimac4_result.dose.vcf"
dosage = convert_ds(Float64, vcffile)


