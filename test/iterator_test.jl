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
    next_state = 1
    final_state = 2
    row = []
    while next_state !== final_state
        row, next_state = iterate(iterator, next_state)
    end

    @test row == Any[("22", 20000086, ["rs138720731"], "T", ["C"], 100.0)]
    @test next_state == 2

    next_state = 1
    final_state = 3
    row = []
    while next_state !== final_state
        if next_state == 2
            row = []
        end
        row, next_state = iterate(iterator, next_state)
    end

    @test row == Any[("22", 20000146, ["rs73387790"], "G", ["A"], 100.0)]
    @test next_state == 3

    # Clean up the temporary VCF file
    rm(vcf_file)
end

@testset "iterator(vcf)" begin
    test_iteration()
end 