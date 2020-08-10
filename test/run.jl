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

cd("/Users/biona001/.julia/dev/MendelImpute/simulation")
tgtfile = "unphase.target_masked.vcf"

Xmm = Mmap.mmap(open(tgtfile, "r"))
String(Xmm[1:10]) # this is ##fileform









# try jiahao's convert function

using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using BenchmarkTools

vcffile = "target.chr18.typedOnly.maf0.1.masked.vcf.gz"
Is, Js, Vs = unsafe_convert_gt(vcffile)




# try convert based on splitting
using Revise
using VCFTools
using MendelImpute
using GeneticVariation
using Random
using StatsBase
using CodecZlib
using ProgressMeter
using BenchmarkTools

# check answer 
vcffile = "target.chr18.typedOnly.maf0.1.masked.vcf.gz"
@time X_true = convert_gt(UInt8, vcffile, trans=true); # 4.669417 seconds (42.34 M allocations: 2.970 GiB, 8.78% gc time)
@time X = unsafe_convert_gt2(vcffile);                 # 2.756179 seconds (39.72 M allocations: 1.707 GiB, 9.44% gc time)
@time X2 = unsafe_convert_gt3(vcffile);                # 2.737134 seconds (40.43 M allocations: 1.722 GiB, 4.26% gc time)
all(skipmissing(X_true .== X))
all(skipmissing(X_true .== X2))


# 1 thread
@btime unsafe_convert_gt2(vcffile); #2.369 s (39085435 allocations: 1.68 GiB)
@btime unsafe_convert_gt3(vcffile); #2.445 s (39680389 allocations: 1.69 GiB)
@btime convert_gt(UInt8, vcffile, trans=true); #2.150 s (33984323 allocations: 2.56 GiB)


# 4 threads
@btime unsafe_convert_gt2(vcffile); # 2.289 s (39085435 allocations: 1.68 GiB)
@btime unsafe_convert_gt3(vcffile); # 1.313 s (39684259 allocations: 1.69 GiB)
@btime convert_gt(UInt8, vcffile, trans=true); # 2.142 s (33984323 allocations: 2.56 GiB)


# 8 threads
@btime unsafe_convert_gt2(vcffile); # 2.327 s (39085435 allocations: 1.68 GiB)
@btime unsafe_convert_gt3(vcffile); # 1.103 s (39689223 allocations: 1.69 GiB)
@btime convert_gt(UInt8, vcffile, trans=true); # 2.234 s (33984323 allocations: 2.56 GiB)




using ProfileView
@profview unsafe_convert_gt3(vcffile)
@profview unsafe_convert_gt3(vcffile)


#test on data with large p
vcffile = "target.chr18.typedOnly.maf0.0005.masked.vcf.gz"
@time unsafe_convert_gt2(vcffile); #11.903016 seconds (196.12 M allocations: 8.412 GiB, 9.62% gc time)
@time unsafe_convert_gt3(vcffile); #6.058362 seconds (199.16 M allocations: 8.461 GiB, 18.60% gc time)
@time convert_gt(UInt8, vcffile, trans=true); #12.555977 seconds (170.52 M allocations: 12.864 GiB, 15.13% gc time)


#test on data with large n
vcffile = "ref.chr18.excludeTarget.vcf.gz"
@time unsafe_convert_gt2(vcffile); #282.177077 seconds (4.13 G allocations: 197.932 GiB, 5.47% gc time)
@time unsafe_convert_gt3(vcffile); #115.613233 seconds (4.13 G allocations: 197.981 GiB, 16.93% gc time)
@time convert_gt(UInt8, vcffile, trans=true); #249.568465 seconds (4.10 G allocations: 309.245 GiB, 4.93% gc time)




stream = GzipDecompressorStream(open(vcffile))
line = readline(stream)
line = readline(stream)
line = readline(stream)
line = readline(stream)
line = readline(stream)
line = readline(stream)


sample_data .= split(line, "\t")

# split by field. Assume GT field comes before all other field
fn1 = x -> split(x, ":")
split_words = map(fn1, datas[10:end])

fn2 = strs-> begin 
    println(strs)
end
map(fn2, split_words)

# profile convert_gt
using Revise
using ProfileView
using VCFTools
using GeneticVariation
using Random
using BenchmarkTools
cd("/Users/biona001/.julia/dev/VCFTools/test")
f = "test.08Jun17.d8b.vcf.gz"


@profview convert_gt(UInt8, f)
@profview convert_gt(UInt8, f)
