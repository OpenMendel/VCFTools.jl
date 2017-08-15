using Distributions, GeneticVariation, GZip

export gtstats

"""
    gtstats(vcffile, [out])

Calculate genotype statistics for each marker in a VCF file with GT field data.
It is similar to the gsstats.jar utility  <https://faculty.washington.edu/browning/beagle_utilities/utilities.html#gtstats>.

# Input
- `vcffile`: VCF file, ending with .vcf or .vcf.gz
- `out`: IOStream for output. Defulat is the `STDOUT`
One line with 15 tab-delimited fiels is written per marker:
    - 1-8)  VCF fixed fields (CHROM, POS, ID, REF, ALT, QUAL, FILT, INFO)
    -   9)  Missing genotype count
    -  10)  Missing genotype frequency
    -  11)  Non-REF allele count
    -  12)  Non-REF allele frequency
    -  13)  Minor allele count             (REF allele vs Non-REF alleles)
    -  14)  Minor allele frequency         (REF allele vs Non-REF alleles)
    -  15)  HWE P-value                    (REF allele vs Non-REF alleles)

# Output
- `lines`: number of lines written to `out`
"""
function gtstats(vcffile::AbstractString, out::IO=STDOUT)
    # open VCF file
    if vcffile[(end - 3):end] == ".vcf"
        reader = VCF.Reader(open(vcffile, "r"))
    elseif vcffile[(end - 6):end] == ".vcf.gz"
        reader = VCF.Reader(GZip.open(vcffile, "r"))
    else
        throw(ArgumentError("VCF filename should end with vcf or vcf.gz"))
    end
    # loop over records
    samples = length(VCF.header(reader).sampleID)
    lines = 0
    for record in reader
        # if no "GT" field, skip this record
        VCF.findgenokey(record, "GT") == 0 && continue
        # calcuate summary statistics
        lines += 1
        missings = refhomzygs = althomzygs = 0
        for i in 1:samples
            gt = VCF.genotype(record, i, "GT")
            if gt == "."
                missings += 1
            elseif gt == "0/0" || gt == "0|0"
                althomzygs += 1
            elseif gt == "1/1" || gt == "1|1"
                refhomzygs += 1
            end
        end
        heterozygs   = samples - missings - refhomzygs - althomzygs
        missfreq     = missings / samples
        altalleles   = heterozygs + 2althomzygs
        refalleles   = heterozygs + 2refhomzygs
        altfreq      = altalleles / (refalleles + altalleles)
        reffreq      = 1 - altfreq
        minoralleles = altfreq < 0.5? altalleles : refalleles
        maf          = altfreq < 0.5? altfreq : reffreq
        # Hardy-Weinburg equilibrium test
        if minoralleles == 0
            hwepval = 1.0
        else
            refhomzygs_e = (samples - missings) * reffreq * reffreq
            althomzygs_e = (samples - missings) * altfreq * altfreq
            heterozygs_e = samples - missings - refhomzygs_e - althomzygs_e
            hwestat = (refhomzygs - refhomzygs_e)^2 / refhomzygs_e +
                      (althomzygs - althomzygs_e)^2 / althomzygs_e +
                      (heterozygs - heterozygs_e)^2 / heterozygs_e
            hwepval = ccdf(Chi(1), hwestat)
        end
        # output
        write(out, record.data[1:record.format[1][1]-2], '\t')
        print(out, missings, '\t', missfreq, '\t', altalleles, '\t',
        altfreq, '\t', minoralleles, '\t', maf, '\t', hwepval, '\n')
    end
    return lines
end

"""
    gtstats(vcffile, out)

Calculate genotype statistics for each marker in a VCF file with GT field data.
Output is written to ta file specified by `out`.
"""
function gtstats(vcffile::AbstractString, out::AbstractString)
    ofile = open(out, "w")
    gtstats(vcffile, ofile)
    close(ofile)
end
