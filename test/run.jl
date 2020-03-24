#trying out VCFTools.jl and GeneticVariation.jl
using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
reader = VCF.Reader(openvcf(vcffile, "r"))
record = read(reader)

record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
record = VCF.Record("20\t14370\t.\tG\tA\t29\tPASS;fdsa\t.\tDS\t1.99\t1.01")
record = VCF.Record("20\t14370\t.\tG\tA\t29\tPASS\t.\tGT\t1/0\t0/0")
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

#import data
cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf.gz"
samples = nsamples(vcffile)
records = nrecords(vcffile)
reader = VCF.Reader(openvcf(vcffile, "r"))
record = read(reader)

# define masks
Random.seed!(123)
record_mask = trues(records)
# sample_mask = bitrand(samples)
sample_mask = trues(samples)
sample_mask[1] = false
sample_mask[end] = false

des = "filtered.test.08Jun17.d8b.vcf.gz"
VCFTools.filter(vcffile, record_mask, sample_mask, des=des)
reader = VCF.Reader(openvcf(des, "r"))
record = read(reader)
record.genotype

nsamples(des)

X = convert_gt(Float32, vcffile)
X_filter = convert_gt(Float32, "filtered." * vcffile)

@benchmark VCFTools.filter(vcffile, record_mask, sample_mask)
#new: 149.763 ms, 214.48 MiB
#old: 134.655 ms, 177.42 MiB

des = "filter." * vcffile
reader = VCF.Reader(openvcf(vcffile, "r"))
writer = VCF.Writer(openvcf(des, "w"), filter_header(reader, sample_mask))
    
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






# test overwrite_alt_ref_allele
using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools

function overwrite_alt_ref_allele(
    src::AbstractString,
    des::AbstractString,
    ref::String = "A",
    alt::String = "T"
    )

    # create input and output VCF files
    reader = VCF.Reader(openvcf(src, "r"))
    writer = VCF.Writer(openvcf(des, "w"), VCF.header(reader))
    
    # define ref/alt allele
    dict = Dict([("A", 0x41), ("T", 0x54), ("C", 0x43), ("G", 0x47), ("N", 0x4e)])
    haskey(dict, ref) || error("ref allele must be A, T, C, G, or N.")
    haskey(dict, alt) || error("ref allele must be A, T, C, G, or N.")
    ref, alt = dict[ref], dict[alt]
    for record in reader
        record.data[record.alt[1]] .= alt
        record.data[record.ref[1]] = ref
        write(writer, record)
    end
    flush(writer); close(reader); close(writer)
end

src = "full_test.vcf"
reader = VCF.Reader(openvcf(src, "r"))
record = read(reader)
# overwrite_alt_ref_allele("test.vcf")
overwrite_alt_ref_allele("test.vcf", "test_result.vcf")



# Simulate genotype matrix using real haplotypes 
# 1000 Genomes Phase 1, chr 22 data downloaded from: https://genome.sph.umich.edu/wiki/Minimac4
using Revise
using MendelImpute
using VCFTools
using BenchmarkTools

cd("/Users/biona001/.julia/dev/MendelImpute/simulation/compare4/VCF_Files")

#2184 haplotypes each with 365644 SNPs -> 99MB of RAM
@time H = convert_ht(Bool, "ALL.chr22.phase1_v3.snps_indels_svs.genotypes.all.noSingleton.vcf.gz")

#simulate 100 samples each with 365644 SNPs
X = simulate_genotypes(H', 100)



# test convert_gt transpose

using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools
using Test

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"

A  = convert_gt(Float64, vcffile)
At = convert_gt(Float64, vcffile, trans=true)
@test all(At .== A')
# size(A) = (191, 1356)

@benchmark convert_gt(Float64, vcffile) # 63.926 ms, 104.59 MiB, 1070395 alloc
@benchmark convert_gt(Float64, vcffile, trans=true)# 65.772 ms, 104.59 MiB, 1070396 alloc

# test if eof(reader) is working
out = Matrix{Union{Float64, Missing}}(undef, 191, 1400)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_gt!(out, reader)
@test all(ismissing.(out[:, 1357:end]))

out = Matrix{Union{Float64, Missing}}(undef, 1400, 191)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_gt_trans!(out, reader)
@test all(ismissing.(out[1357:end, :]))


# test impute, center, scale
A = Matrix{Union{Float64, Missing}}(undef, 191, 1400)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_gt!(A, reader, impute=true, center=true, scale=true)
@test isapprox(mean(A[:, 5]), 0, atol=10)
@test isapprox(var(A[:, 5]), 1, atol=10)

At = Matrix{Union{Float64, Missing}}(undef, 1400, 191)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_gt_trans!(At, reader, impute=true, center=true, scale=true)
@test isapprox(mean(At[5, :]), 0, atol=10)
@test isapprox(var(At[5, :]), 1, atol=10)






# test convert_ht transpose

using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools
using Test

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"

H  = convert_ht(Float64, vcffile)
Ht = convert_ht(Float64, vcffile, trans=true)
@test all(Ht .== H')
# size(H) = (382, 1356)

@benchmark convert_ht(Float64, vcffile) # 40.711 ms, 58.89 MiB, 552403 alloc
@benchmark convert_ht(Float64, vcffile, trans=true) # 43.491 ms, 58.89 MiB, 552404 alloc

# BELOW NOT DONE YET!!

# test if eof(reader) is working
out = Matrix{Float64}(undef, 382, 1400)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_ht!(out, reader)
@test size(out) == (382, 1357)

out = Matrix{Float64}(undef, 1400, 382)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_gt_trans!(out, reader)
@test size(out) == (382, 1357)







# test impute, center, scale
A = Matrix{Union{Float64, Missing}}(undef, 191, 1400)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_gt!(A, reader, impute=true, center=true, scale=true)
@test isapprox(mean(A[:, 5]), 0, atol=10)
@test isapprox(var(A[:, 5]), 1, atol=10)

At = Matrix{Union{Float64, Missing}}(undef, 1400, 191)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_gt_trans!(At, reader, impute=true, center=true, scale=true)
@test isapprox(mean(At[5, :]), 0, atol=10)
@test isapprox(var(At[5, :]), 1, atol=10)

