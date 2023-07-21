# Test the iteration behavior of the VCFIterator
function test_iteration()
    # Create a temporary VCF file for testing
    vcf_file = "test.vcf"
    open(vcf_file, "w") do io
        write(io, "##fileformat=VCFv4.3\n")
        write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
        write(io, "chr1\t100\t.\tA\tC\t10.0\n")
        write(io, "chr2\t200\t.\tG\tT\t.\n")
        write(io, "chr3\t300\t.\tC\tG\t20.5\n")
    end

    data = read(vcf_file, String)
    println(data)
    
    # Create the VCFIterator
    iterator = vcfiterator(vcf_file)
    
    # Iterate over the variants and collect the rows
    rows = []
    next_state = 1
    while next_state !== nothing
        row, next_state = iterate(iterator, next_state)
        push!(rows, row)
    end
    
    # Expected result
    expected_rows = [
        VCFRow("chr1", 100, ".", "A", "C", 10.0),
        VCFRow("chr2", 200, ".", "G", "T", NaN),
        VCFRow("chr3", 300, ".", "C", "G", 20.5)
    ]
    
    # Test that the iterated rows match the expected result
    @test rows == expected_rows
    
    # Test that the length of the iterator matches the number of variants
    @test length(iterator) == 3
    
    # Clean up the temporary VCF file
    rm(vcf_file)
end

# Run the tests
@testset "VCFIterator tests" begin
    @testset "Iteration behavior" begin
        test_iteration()
    end
end