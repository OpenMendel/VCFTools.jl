# Test the iteration behavior of the VCFIterator

# vcf file 
# read the data and print 
# create the VCFIterator 
# expected rows vs rows 

function test_iteration()
    
    # Create a temporary VCF file for testing
    isfile("test.08Jun17.d8b.vcf.gz") || Downloads.download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
    joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz"))

    vcf_file = joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz")
        
    # Create the VCFIterator
    iterator = vcfiterator(vcf_file)
    
    # Iterate over the variants and collect the rows
    rows = []
    next_state = 1
    final_state = 3
    while next_state !== final_state
        row, next_state = iterate(iterator, next_state)
        push!(rows, row)
    end
    
    # Expected result
    expected_rows = []
    for n in 1:final_state
        push!(expected_rows, vcf_file[i])
    end 

    # Test that the iterated rows match the expected result
    @test rows == expected_rows
        
    # Clean up the temporary VCF file
    rm(vcf_file)
end

@testset "iterator(vcf)" begin
    test_iteration()
end 