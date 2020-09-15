#trying out VCFTools.jl and GeneticVariation.jl
using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
reader = VCF.Reader(openvcf(vcffile, "r"))
sampleID = VCF.header(reader).sampleID
record = read(reader)

record = VCF.Record("20\t14370\trs6054257\tG\tAAAA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
record = VCF.Record("20\t14370\t.\tG\tA\t29\tPASS;fdsa\t.\tDS\t1.99\t1.01")
record = VCF.Record("20\t14370\t.\tG\tA\t29\tPASS\t.\tGT\t1/0\t0/0")
record.genotype

#test convert_ht
using Revise
using GeneticVariation
using VCFTools

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"
H = convert_ht(Float32, vcffile, has_missing=true)

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
# old = 671.380 ms, 326.66 MiB, 3631762 alloc
# new = 573.432 ms, 257.75 MiB, 3352787 alloc

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
@benchmark mask_gt(vcffile, masks) #40.840 ms, 40.09 MiB, 521257 alloc


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

# cd("/Users/biona001/.julia/dev/VCFTools/test")
# vcffile = "test.08Jun17.d8b.vcf"
vcffile = "target.vcf"

H  = convert_ht(Float64, vcffile)
Ht = convert_ht(Float64, vcffile, trans=true)
@test all(Ht .== H')
# size(H) = (382, 1356)

@benchmark convert_ht(Float64, vcffile) # 44.523 ms, 58.89 MiB, 552403 alloc
@benchmark convert_ht(Float64, vcffile, trans=true) # 44.900 ms, 58.89 MiB, 552404 alloc

# test if eof(reader) is working
out = Matrix{Float64}(undef, 382, 1400)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_ht!(out, reader)
@test size(out) == (382, 1357)
@test all(out[:, 1358:end] .== 0.0)

out = Matrix{Float64}(undef, 1400, 382)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_ht_trans!(out, reader)
@test size(out) == (1400, 382)
@test all(out[1358:end, :] .== 0.0)



# test convert_ds transpose

using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools
using Test

cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf"

D  = convert_ds(Float64, vcffile)
Dt = convert_ds(Float64, vcffile, trans=true)
@test all(Dt .== D')
# size(D) = (191, 1356)

@benchmark convert_ds(Float64, vcffile) # 102.551 ms, 183.63 MiB, 2106379 alloc
@benchmark convert_ds(Float64, vcffile, trans=true) # 114.470 ms ms, 183.63 MiB, 2106380 alloc

# test if eof(reader) is working
out = Matrix{Union{Missing, Float64}}(undef, 191, 1400)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_ds!(out, reader)
@test size(out) == (191, 1400)
@test all(ismissing.(out[:, 1358:end]))

out = Matrix{Union{Missing, Float64}}(undef, 1400, 191)
reader = VCF.Reader(openvcf(vcffile, "r"))
copy_ds_trans!(out, reader)
@test size(out) == (1400, 191)
@test all(ismissing.(out[1358:end, :]))






using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
cd("/Users/biona001/.julia/dev/MendelImpute/simulation")

# impute 
tgtfile = "./compare6/target_masked.vcf.gz"
reffile = "./compare6/haplo_ref.vcf"
outfile = "./compare6/imputed_target.vcf.gz"
width   = 800

H  = convert_ht(Float64, reffile, trans=true)
Hb = convert_ht(Bool, reffile, trans=true)

Hb_c = convert(Matrix{Float64}, Hb)
all(Hb_c .== Hb)










using Revise
using LoopVectorization
using Random
using LinearAlgebra
using BenchmarkTools

function gemv_naive!(c, A, b)
    @inbounds for j in 1:size(A, 2)
        @simd for i in 1:size(A, 1)
            c[i] += A[i, j] * b[j]
        end
    end
end

function gemv_avx!(c, A, b)
    @avx for j in 1:size(A, 2), i in 1:size(A, 1)
        c[i] += A[i, j] * b[j]
    end
end

function gemv_avx2!(c, A, b)
    @avx for i in 1:size(A, 1)
        ci = zero(eltype(c))
        for j in 1:size(A, 2)
            ci += A[i, j] * b[j]
        end
        c[i] = ci
    end
end

Random.seed!(2020)
n = 10000
A = bitrand(n, n)
b = rand(n)
c = zeros(n)
c_avx = zeros(n)
c_avx2 = zeros(n)
c_true = zeros(n)

# check correctness
mul!(c_true, A, b)
gemv_avx!(c_avx, A, b)
gemv_avx2!(c_avx2, A, b)
gemv_naive!(c, A, b)
[c c_avx c_avx2 c_true]

# efficiency (A = bitmatrix, 1000x1000)
@benchmark mul!(c_true, A, b)      # 3.411 ms     (Julia's default)
@benchmark gemv_naive!(c, A, b)    # 4.230 ms     (Looping using @simd and @inbounds)
@benchmark gemv_avx!(c_avx, A, b)  # 566.309 μs   (using LoopVectorization)
@benchmark gemv_avx2!(c_avx, A, b) # 572.952 μs   (using LoopVectorization with accumulator)

# efficiency (A = bitmatrix, 10000x10000)
@benchmark mul!(c_true, A, b)      # 341.411 ms  (Julia's default)
@benchmark gemv_naive!(c, A, b)    # 424.198 ms  (Looping using @simd and @inbounds)
@benchmark gemv_avx!(c_avx, A, b)  # 77.201 ms   (using LoopVectorization)
@benchmark gemv_avx2!(c_avx, A, b) # 73.227 ms   (using LoopVectorization with accumulator)

# efficiency (A = bitmatrix, 50000x50000)
@benchmark mul!(c_true, A, b)      # 8.999 s   (Julia's default)
@benchmark gemv_naive!(c, A, b)    # 10.207 s  (Looping using @simd and @inbounds)
@benchmark gemv_avx!(c_avx, A, b)  # 2.685 s   (using LoopVectorization)
@benchmark gemv_avx2!(c_avx, A, b) # 2.197 s   (using LoopVectorization with accumulator)

# efficiency (A = bitmatrix, 100000x100000)
@time mul!(c_true, A, b)      # 37.167032 s   (Julia's default)
@time gemv_naive!(c, A, b)    # 42.665357 s   (Looping using @simd and @inbounds)
@time gemv_avx!(c_avx, A, b)  # 17.452804 s   (using LoopVectorization)
@time gemv_avx2!(c_avx, A, b) # 17.881693 s   (using LoopVectorization with accumulator)

# efficiency (A = Matrix{Float64}, 1000x1000)
BLAS.set_num_threads(1)
@benchmark mul!(c_true, A, b)      # 137.602 μs   (Julia's default: calls BLAS)
@benchmark gemv_naive!(c, A, b)    # 155.327 μs   (Looping using @simd and @inbounds)
@benchmark gemv_avx!(c_avx, A, b)  # 174.087 μs   (using LoopVectorization)
@benchmark gemv_avx2!(c_avx, A, b) # 189.796 μs   (using LoopVectorization with accumulator)

# efficiency (A = Matrix{Float64}, 10000x10000)
BLAS.set_num_threads(1)
@benchmark mul!(c_true, A, b)      # 41.293 ms   (Julia's default: calls BLAS)
@benchmark gemv_naive!(c, A, b)    # 47.362 ms   (Looping using @simd and @inbounds)
@benchmark gemv_avx!(c_avx, A, b)  # 97.696 ms   (using LoopVectorization)
@benchmark gemv_avx2!(c_avx, A, b) # 99.377 ms   (using LoopVectorization with accumulator)






using VectorizationBase: gesp, stridedpointer
using Random
using LoopVectorization
using BenchmarkTools
using LinearAlgebra

function gemv_avx!(c, A, b)
    @avx for j in 1:size(A, 2), i in 1:size(A, 1)
        c[i] += A[i, j] * b[j]
    end
end
function gemv_avx_512by512!(c, A, b)
    @avx for j in 1:512, i in 1:512
        c[i] += A[i, j] * b[j]
    end
end
function gemv_tile!(c, A, b)
    M, N = size(A)
    Miter = M >>> 9 # fast div(M, 512)
    Mrem = M & 511 # fast rem(M, 512)
    Niter = N >>> 9
    Nrem = N & 511
    GC.@preserve c A b for n in 0:Niter-1
        for m in 0:Miter-1
            gemv_avx_512by512!(
                gesp(stridedpointer(c), (512m,)),
                gesp(stridedpointer(A), (8m, 8n)),
                gesp(stridedpointer(b), (512n,))
            )
        end
        if mrem > 0
            
        end
    end
    # TODO: handle nrem
end

n = 1 << 14 # 2^14, multiple of 512, because we're not handling remainders
A = bitrand(n, n);
b = rand(n); c1 = zero(b); c2 = zero(b);
gemv_tile!(c1, A, b);
gemv_avx!(c2, A, b);
c1 ≈ c2

@benchmark gemv_tile!($c1, $A, $b) # 200.051 ms  (cache aware avx)
@benchmark gemv_avx!($c2, $A, $b)  # 449.613 ms  (naive avx)
@benchmark mul!($c1, $A, $b)       # 1.133 s     (Julia built-in)







using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools

#import data
cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf.gz"
snps = nrecords(vcffile)
chunks = ceil(Int, snps / 100)
snps_per_chunk = ceil(Int, snps / chunks)
(chunk - 1)*snps_per_chunk+1:(chunk * snps_per_chunk) #ranges

full_X = convert_gt(Float32, vcffile, trans=true)

Xreader = VCF.Reader(openvcf(vcffile, "r"))
X = Matrix{Union{Float32, Missing}}(undef, snps_per_chunk, nsamples(vcffile))
for chunk in 1:chunks
    copy_gt_trans!(X, Xreader)
    println(all(X .== full_X[(chunk - 1)*snps_per_chunk+1:(chunk * snps_per_chunk), :]))
end




# try data import strategy outlined in https://github.com/JuliaLang/julia/issues/34195#issuecomment-569343440

using Mmap
using VCFTools
using CodecZlib

using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools

#import data
cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "test.08Jun17.d8b.vcf.gz"

H, sampleID, record_chr, record_pos, record_ids, record_ref, record_alt = convert_ht(Float32, vcffile, save_snpinfo(), trans=true)
@code_warntype convert_ht(Float32, vcffile, save_snpinfo(), trans=true)
H2 = convert_ht(Float32, vcffile, trans=true)
all(H .== H2)
@code_warntype convert_ht(Float32, vcffile, trans=true)


X, sampleID, record_chr, record_pos, record_ids, record_ref, record_alt = convert_gt(Float32, vcffile, save_snpinfo(), trans=true)
@code_warntype convert_gt(Float32, vcffile, save_snpinfo(), trans=true)
X2 = convert_gt(Float32, vcffile, trans=true)
all(X .== X2)
@code_warntype convert_gt(Float32, vcffile, trans=true)


# test aim selection
using Revise
using GeneticVariation
using Random
using VCFTools
using BenchmarkTools

# get each sample's population origin
cd("/Users/biona001/.julia/dev/MendelImpute/data/1000_genome_phase3_v5/country_origin")
sampleID_to_population = Dict{String, String}()
for population in readdir("data/")
    for sample in readdir("data/" * population)
        sample == ".DS_Store" && continue
        sampleID_to_population[sample] = population
    end
end

# chr22 in 1000 genomes data
cd("/Users/biona001/.julia/dev/VCFTools/test")
vcffile = "chr22.1kg.phase3.v5a.vcf.gz"
pvals = aim_select(vcffile, sampleID_to_population)
