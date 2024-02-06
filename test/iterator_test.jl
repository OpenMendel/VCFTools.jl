# Test the iteration behavior of the VCFIterator

# vcf file 
# read the data and print 
# create the VCFIterator 
# expected rows vs rows 

"""
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
    row = VCFRow("", 0, [""], "", [""], 0)
    while next_state !== final_state
        row, next_state = iterate(iterator, next_state)
    end

   #@test row == VCFRow("22", 20000086, ["rs138720731"], "T", ["C"], 100.0)
   @test string(row.CHROM) == string("22")
   @test Int(row.POS) == Int(20000086)
   @test string(row.ID) == string(["rs138720731"])
   @test string(row.REF) == string("T")
   @test string(row.ALT) == string(["C"])
   @test float(row.QUAL) == float(100.0)
   @test next_state == 2

   next_state = 1
    final_state = 3
    row = VCFRow("", 0, [""], "", [""], 0)
    while next_state !== final_state
        if next_state == 2
            row = VCFRow("", 0, [""], "", [""], 0)
        end
        row, next_state = iterate(iterator, next_state)
    end
   
    # return a tuple not an array 

    #@test row == VCFRow("22", 20000146, ["rs73387790"], "G", ["A"], 100.0)
    @test string(row.CHROM) == string("22")
    @test Int(row.POS) == Int(20000146)
    @test string(row.ID) == string(["rs73387790"])
    @test string(row.REF) == string("G")
    @test string(row.ALT) == string(["A"])
    @test float(row.QUAL) == float(100.0)
    @test next_state == 3

    # Clean up the temporary VCF file
    rm(vcf_file)
end

@testset "iterator(vcf)" begin
    test_iteration()
end 

# put genetic variant base tests here 
"""

# Create a temporary VCF file for testing
isfile("test.08Jun17.d8b.vcf.gz") || Downloads.download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz"))

vcf_file = joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz")
vcf_iterator = VCFIterator(vcf_file)

for i in 1:50
    vcf_row, _ = iterate(vcf_iterator, i)

    # Print the genotype information for the current record
    println("Genotypes for record $i: ", vcf_row.GENOTYPE)
end

row = VCFRow("1", 9001061, [], "T", ["G"], 32.0, [1.0, 1.0])



"""
vcf_row = VCFRow("22", 20000086, ["rs138720731"], "T", ["C"], 100.0, Union{Missing, Float64}[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
data = VCFData(vcf_file, vcf)
@test maf(data, vcf_row) == 0.0
"""
