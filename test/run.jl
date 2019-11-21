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
record.genotype

#test convert_ht
using Revise
using GeneticVariation
using VCFTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
H = convert_ht(Float32, vcffile)



# test filter
using Revise
using GeneticVariation
using Random
using VCFTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf.gz"
samples = nsamples(vcffile)
records = nrecords(vcffile)

Random.seed!(123)
record_index = bitrand(records)
sample_index = bitrand(samples)
VCFTools.filter(vcffile, record_index, sample_index)



des = "filter." * vcffile
reader = VCF.Reader(openvcf(vcffile, "r"))
writer = VCF.Writer(openvcf(des, "w"), filter_header(reader, sample_index))
    

reader = VCF.Reader(openvcf(vcffile, "r"))
record = read(reader)


