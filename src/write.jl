"""
    write_vcf(outfile, x1, x2, [chr], [pos], [ref], [alt], [sampleID], [SNPids])

Writes haplotypes in `x1` and `x2` into VCF file, tracking phase information. 

# Inputs
- `outfile`: Output file name (ends in `.vcf` or `.vcf.gz`)
- `x1`: `n × p` matrix of the 1st haplotype for each sample. Each row is a haplotype, 
    and each entry is 0 or 1
- `x2`: `n × p` matrix of the 2nd haplotype for each sample. `x = x1 + x2`. 
    Each entry is 0 or 1

# Optional Inputs
- `chr`: Vector of String. `chr[i]` is chromosome number of SNP `i`. 
- `pos`: Vector of Int. `pos[i]` is position of SNP `i`. 
- `ref`: Vector of String. `ref[i]` is reference allele of SNP `i`. 
- `alt`: Vector of String. `alt[i]` is alternate allele of SNP `i`. 
- `sampleID`: Vector of String. `sampleID[i]` is reference allele of SNP `i`. 
- `SNPids`: Vector of String. `SNPids[i]` is reference allele of SNP `i`. 

# Reference
For multithreaded write, see
https://github.com/OpenMendel/MendelImpute.jl/blob/master/src/impute.jl#L23
"""
function write_vcf(
    outfile::AbstractString,
    x1::AbstractMatrix,
    x2::AbstractMatrix;
    chr::AbstractVector{String} = ["1"     for i in 1:size(x1, 2)],
    pos::AbstractVector{Int}    = [100i    for i in 1:size(x1, 2)],
    ref::AbstractVector{String} = ["A"     for i in 1:size(x1, 2)],
    alt::AbstractVector{String} = ["T"     for i in 1:size(x1, 2)],
    sampleID::AbstractVector{String} = [string(i) for i in 1:size(x1, 1)],
    SNPids::AbstractVector{String} = ["snp$i" for i in 1:size(x1, 2)],
    )
    n, p = size(x1)
    (n, p) == size(x2) || error("x1 and x2 have different dimensions!")
    io = openvcf(outfile, "w")
    pb = PipeBuffer()
    
    # minimum meta info
    print(pb, "##fileformat=VCFv4.2\n")
    print(pb, "##source=VCFTools.jl\n")
    print(pb, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # header line (includes sample ID)
    print(pb, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for id in sampleID
        print(pb, "\t", id)
    end
    print(pb, "\n")
    bytesavailable(pb) > 1048576 && write(io, take!(pb))
    
    # loop over genotypes
    pmeter = Progress(p, 1, "Writing VCF...")
    @inbounds for i in 1:p
        # write meta info (chrom/pos/snpid/ref/alt/imputation-quality)
        print(pb, chr[i], "\t", string(pos[i]), "\t", SNPids[i], "\t", 
            ref[i], "\t", alt[i], "\t.\tPASS\t.\tGT")
        # print ith record
        write_snp!(pb, x1, x2, i) 
        bytesavailable(pb) > 1048576 && write(io, take!(pb))
        next!(pmeter)
    end
    write(io, take!(pb))
    close(io); close(pb) # close io and buffer
end

"""
Helper function for saving the `i`th record (SNP), tracking phase information.
Here `X = X1 + X2`.
"""
function write_snp!(pb::IOBuffer, X1::AbstractMatrix, X2::AbstractMatrix, i::Int)
    x1 = @view(X1[:, i])
    x2 = @view(X2[:, i])
    n = length(x1)
    @assert n == length(x2)
    @inbounds for j in 1:n
        if x1[j] == x2[j] == 0
            print(pb, "\t0|0")
        elseif x1[j] == 0 && x2[j] == 1
            print(pb, "\t0|1")
        elseif x1[j] == 1 && x2[j] == 0
            print(pb, "\t1|0")
        elseif x1[j] == 1 && x2[j] == 1
            print(pb, "\t1|1")
        else
            error("phased genotypes can only be 0|0, 0|1, 1|0 or 1|1 but
                got $(x1[j])|$(x2[j])")
        end
    end
    print(pb, "\n")
    return nothing
end

"""
    write_vcf(outfile, x, [chr], [pos], [ref], [alt], [sampleID], [SNPids])

Writes genotypes in `x` into VCF file, not tracking phase information. 

# Inputs
- `outfile`: Output file name (ends in `.vcf` or `.vcf.gz`)
- `x`: `n × p` genotype matrix. Each row is a sample, and each entry is 0, 1, or 2. 

# Optional Inputs
- `chr`: Vector of String. `chr[i]` is chromosome number of SNP `i`. 
- `pos`: Vector of Int. `pos[i]` is position of SNP `i`. 
- `ref`: Vector of String. `ref[i]` is reference allele of SNP `i`. 
- `alt`: Vector of String. `alt[i]` is alternate allele of SNP `i`. 
- `sampleID`: Vector of String. `sampleID[i]` is reference allele of SNP `i`. 
- `SNPids`: Vector of String. `SNPids[i]` is reference allele of SNP `i`. 

# Reference
For multithreaded write, see
https://github.com/OpenMendel/MendelImpute.jl/blob/master/src/impute.jl#L23
"""
function write_vcf(
    outfile::AbstractString,
    x::AbstractMatrix;
    chr::AbstractVector{String} = ["1"     for i in 1:size(x, 2)],
    pos::AbstractVector{Int}    = [100i    for i in 1:size(x, 2)],
    ref::AbstractVector{String} = ["A"     for i in 1:size(x, 2)],
    alt::AbstractVector{String} = ["T"     for i in 1:size(x, 2)],
    sampleID::AbstractVector{String} = [string(i) for i in 1:size(x, 1)],
    SNPids::AbstractVector{String} = ["snp$i" for i in 1:size(x, 2)],
    )
    n, p = size(x)
    io = openvcf(outfile, "w")
    pb = PipeBuffer()
    
    # minimum meta info
    print(pb, "##fileformat=VCFv4.2\n")
    print(pb, "##source=VCFTools.jl\n")
    print(pb, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

    # header line (includes sample ID)
    print(pb, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for id in sampleID
        print(pb, "\t", id)
    end
    print(pb, "\n")
    bytesavailable(pb) > 1048576 && write(io, take!(pb))
    
    # loop over genotypes
    pmeter = Progress(p, 1, "Writing VCF...")
    @inbounds for i in 1:p
        # write meta info (chrom/pos/snpid/ref/alt/imputation-quality)
        print(pb, chr[i], "\t", string(pos[i]), "\t", SNPids[i], "\t", 
            ref[i], "\t", alt[i], "\t.\tPASS\t.\tGT")
        # print ith record
        write_snp!(pb, x, i) 
        bytesavailable(pb) > 1048576 && write(io, take!(pb))
        next!(pmeter)
    end
    write(io, take!(pb))
    close(io); close(pb) # close io and buffer
end


"""
Helper function for saving the `i`th record (SNP), not tracking phase information.
Here genotypes will be written as
0 => 0/0
1 => 1/0
2 => 1/1
"""
function write_snp!(pb::IOBuffer, X::AbstractMatrix, i::Int)
    x = @view(X[:, i])
    n = length(x)
    @inbounds for j in 1:n
        if x[j] == 0
            print(pb, "\t0/0")
        elseif x[j] == 1
            print(pb, "\t1/0")
        elseif x[j] == 2
            print(pb, "\t1/1")
        else
            error("unphased genotypes can only be 0, 1, or 2 but got $(x[j])")
        end
    end
    print(pb, "\n")
    return nothing
end
