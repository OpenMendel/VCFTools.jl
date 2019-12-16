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
