"""
    conformgt_by_id(reffile, tgtfile, outfile, chrom, posrange, checkfreq)

Match the VCF records in `tgtfile` to those in `reffile` according to ID.
The function will:
0. Find corresponding VCF records in the target and reference files
0. Exclude target VCF records whose ID cannot be matched to any reference VCF record
0. Exclude target VCF records whose test of equal allele frequency is rejected at significance level `checkfreq`
0. Adjust target VCF records so that chromosome strand and allele order match the VCF reference file
0. The matched VCF records are written into files `outfile.tgt.vcf.gz` and
`outfile.ref.vcf.gz`, both with only "GT" data

# Input
- `reffile`: VCF file with reference genotype (GT) data
- `tgtfile`: VCF file with target genotype (GT) data
- `outfile`: the prefix for output filenames
- `chrom`: chromosome name, must be identical in target and reference files
- `posrange`: position range in the reference file
- `checkfreq`: significance level for testing equal alelle frequencies between mached target and reference records. If the test pvalue is  `≤ checkfreq`, the records are not output. Setting `checkfreq=0` or `checkfreq=false` (default) implies not checking allele frequencies. Setting `checkfreq=1` effectively rejects all tests and no matched records are output

# Output
- `lines`: number of matched VCF records
"""
function conformgt_by_id(
    reffile::AbstractString,
    tgtfile::AbstractString,
    outfile::AbstractString,
    chrom::AbstractString,
    posrange::Range,
    checkfreq::Number = false
    )
    # open reference and target VCF files
    reader_ref = VCF.Reader(openvcf(reffile, "r"))
    reader_tgt = VCF.Reader(openvcf(tgtfile, "r"))
    records_ref, records_tgt = nrecords(reffile), nrecords(tgtfile)
    # create output files
    writer_ref = VCF.Writer(openvcf(join([outfile, ".ref.vcf.gz"]), "w"),
        VCF.header(reader_ref))
    writer_tgt = VCF.Writer(openvcf(join([outfile, ".tgt.vcf.gz"]), "w"),
        VCF.header(reader_tgt))
    # collect IDs in the reference panel
    info("Scan IDs in reference panel"; prefix = "conformgt_by_id: ")
    pbar = ProgressMeter.Progress(records_ref, 1)
    id_list_ref = Vector{String}[]
    record_counter = 0
    for record_ref in reader_ref
        # update progress bar
        record_counter += 1
        ProgressMeter.update!(pbar, record_counter)
        # chromosome not matched
        VCF.chrom(record_ref) == chrom || continue
        # if past search range, break
        pos_ref = VCF.pos(record_ref)
        if pos_ref > last(posrange)
            ProgressMeter.update!(pbar, records_ref)
            break
        end
        # if not in search range, skip this record
        pos_ref < first(posrange) && continue
        # if no "GT" field, skip this record
        VCF.findgenokey(record_ref, "GT") == 0 && continue
        # if no ID, skip this record
        VCF.hasid(record_ref) || continue
        # add ID to list
        id_ref = VCF.id(record_ref)
        if length(id_list_ref) > 0 && id_ref == id_list_ref[end]
            print("Duplicate ID in reference panel: $(id_ref); ")
            print("will use the first occurence\n")
            continue
        end
        push!(id_list_ref, id_ref)
    end
    close(reader_ref)
    # match target IDs to reference IDs
    info("Match target IDs to reference IDs"; prefix = "conformgt_by_id: ")
    # re-set reference reader
    reader_ref = VCF.Reader(openvcf(reffile, "r"))
    record_ref = VCF.Record()
    state_ref = start(reader_ref)
    # loop over target records
    pbar = ProgressMeter.Progress(records_tgt, 1)
    lines = record_counter = 0
    for record_tgt in reader_tgt
        # update progress bar
        record_counter += 1
        ProgressMeter.update!(pbar, record_counter)
        # if no "GT" field, skip this record
        VCF.findgenokey(record_tgt, "GT") == 0 && continue
        # if no ID, skip
        VCF.hasid(record_tgt) || continue
        # if not in reference panel, skip
        # TODO can be optimized here, shrink id_list_ref
        id_tgt = VCF.id(record_tgt)
        id_tgt ∈ id_list_ref || continue
        # search for matching reference record
        while !done(reader_ref, state_ref)
            record_ref, state_ref = next(reader_ref, state_ref)
            VCF.hasid(record_ref) && VCF.id(record_ref) == id_tgt && break
        end
        # check REF/ALT label match
        reflabel_tgt = VCF.ref(record_tgt)
        altlabel_tgt = VCF.alt(record_tgt)
        reflabel_ref = VCF.ref(record_ref)
        reflabel_tgt == reflabel_ref || reflabel_ref ∈ altlabel_tgt || continue
        record_gt_ref, record_gt_tgt = match_gt_allele(record_ref, record_tgt)
        # allele frequency test if requested
        if checkfreq > 0
            _, _, _, n0_ref, n1_ref, = gtstats(record_gt_ref, nothing)
            _, _, _, n0_tgt, n1_tgt, = gtstats(record_gt_tgt, nothing)
            pval = binomial_proportion_test(n0_ref, n1_ref, n0_tgt, n1_tgt)
            pval ≤ checkfreq && continue
        end
        # write to matched target and reference file
        lines += 1
        VCF.write(writer_ref, record_gt_ref)
        VCF.write(writer_tgt, record_gt_tgt)
    end
    close(reader_ref); close(reader_tgt)
    flush(writer_ref); flush(writer_tgt)
    close(VCF.BioCore.IO.stream(writer_ref))
    close(VCF.BioCore.IO.stream(writer_tgt))
    # return
    info("$(lines) records are matched"; prefix = "conformgt_by_id: ")
    return lines
end


"""
    conformgt_by_pos(reffile, tgtfile, outfile, chrom, posrange, checkfreq)

Match the VCF GT records in `tgtfile` to those in `reffile` according to
chromosome position. The function will:
0. Find corresponding VCF records in the target and reference files
0. Exclude target VCF records whose position cannot be matched to any reference VCF record
0. Exclude target VCF records whose test of equal allele frequency is rejected at significance level `checkfreq`
0. Adjust target VCF records so that chromosome strand and allele order match the VCF reference file
0. The matched VCF records are written into files `outfile.tgt.vcf.gz` and
`outfile.ref.vcf.gz`, both with only "GT" data

# Input
- `reffile`: VCF file with reference genotype (GT) data
- `tgtfile`: VCF file with target genotype (GT) data
- `outfile`: the prefix for output filenames
- `chrom`: chromosome name, must be identical in target and reference files
- `posrange`: position range in the reference file
- `checkfreq`: significance level for testing equal alelle frequencies between mached target and reference records. If the test pvalue is  `≤ checkfreq`, the records are not output. Setting `checkfreq=0` or `checkfreq=false` (default) implies not checking allele frequencies. Setting `checkfreq=1` effectively rejects all tests and no matched records are output

# Output
- `lines`: number of matched VCF records
"""
function conformgt_by_pos(
    reffile::AbstractString,
    tgtfile::AbstractString,
    outfile::AbstractString,
    chrom::AbstractString,
    posrange::Range,
    checkfreq::Number = false
    )
    # open reference and target VCF files
    reader_ref = VCF.Reader(openvcf(reffile, "r"))
    reader_tgt = VCF.Reader(openvcf(tgtfile, "r"))
    records_ref, records_tgt = nrecords(reffile), nrecords(tgtfile)
    # create output files
    writer_ref = VCF.Writer(openvcf(join([outfile, ".ref.vcf.gz"]), "w"),
        VCF.header(reader_ref))
    writer_tgt = VCF.Writer(openvcf(join([outfile, ".tgt.vcf.gz"]), "w"),
        VCF.header(reader_tgt))
    # initialize reference reader
    record_ref, state_ref, pos_ref = VCF.Record(), start(reader_ref), 0
    # loop over target records
    info("Match target POS to reference POS"; prefix = "conformgt_by_pos: ")
    pbar = ProgressMeter.Progress(records_tgt, 1)
    lines = record_counter = 0
    for record_tgt in reader_tgt
        # update progress bar
        record_counter += 1
        ProgressMeter.update!(pbar, record_counter)
        # if no "GT" field, skip this record
        VCF.findgenokey(record_tgt, "GT") == 0 && continue
        # chromosome not matched, skip this record
        VCF.chrom(record_tgt) == chrom || continue
        # if past search range, break
        pos_tgt = VCF.pos(record_tgt)
        if pos_tgt > last(posrange)
            ProgressMeter.update!(pbar, records_tgt)
            break
        end
        # if before search range, skip this record
        pos_tgt < first(posrange) && continue
        # if before the current ref position, skip this record
        pos_tgt < pos_ref && continue
        if pos_tgt == pos_ref
            matched = true
        elseif pos_tgt > pos_ref
            matched = false
            # search for matching reference record
            while !done(reader_ref, state_ref)
                record_ref, state_ref = next(reader_ref, state_ref)
                # if no "GT" field, skip this record
                VCF.findgenokey(record_tgt, "GT") == 0 && continue
                # if chromosome not matched, skip this record
                VCF.chrom(record_ref) == chrom || continue
                pos_ref = VCF.pos(record_ref)
                if pos_ref < pos_tgt
                    continue
                elseif pos_ref == pos_tgt
                    matched = true; break
                elseif pos_ref > pos_tgt
                    break
                end
            end
        end
        # if no matching reference record, skip this record
        matched || continue
        # check REF/ALT label match
        reflabel_tgt = VCF.ref(record_tgt)
        altlabel_tgt = VCF.alt(record_tgt)
        reflabel_ref = VCF.ref(record_ref)
        reflabel_tgt == reflabel_ref || reflabel_ref ∈ altlabel_tgt || continue
        record_gt_ref, record_gt_tgt = match_gt_allele(record_ref, record_tgt)
        # allele frequency test if requested
        if checkfreq > 0
            _, _, _, n0_ref, n1_ref, = gtstats(record_gt_ref, nothing)
            _, _, _, n0_tgt, n1_tgt, = gtstats(record_gt_tgt, nothing)
            pval = binomial_proportion_test(n0_ref, n1_ref, n0_tgt, n1_tgt)
            pval ≤ checkfreq && continue
        end
        # write to matched target and reference file
        lines += 1
        VCF.write(writer_ref, record_gt_ref)
        VCF.write(writer_tgt, record_gt_tgt)
    end
    close(reader_ref); close(reader_tgt)
    flush(writer_ref); flush(writer_tgt)
    close(VCF.BioCore.IO.stream(writer_ref))
    close(VCF.BioCore.IO.stream(writer_tgt))
    # return
    info("$(lines) records are matched"; prefix = "conformgt_by_pos: ")
    return lines
end


"""
    match_gt_allele(record1, record2)

Match the REF allele label of record 2 to that of record 1. Return the possibly
REF/ALT flipped record 2 with only "GT" genotype data.
"""
function match_gt_allele(record1::VCF.Record, record2::VCF.Record)
    @assert VCF.findgenokey(record1, "GT") > 0 "record1 has no GT field"
    @assert VCF.findgenokey(record2, "GT") > 0 "record2 has no GT field"
    ref1 = VCF.ref(record1)
    ref2 = VCF.ref(record2)
    alt2 = VCF.alt(record2)
    record1_out = filter_genotype(record1, ["GT"])
    if ref1 == ref2
        # matching REF labels
        record2_out = filter_genotype(record2, ["GT"])
    elseif ref1 ∈ alt2
        # switch REF/ALT label of record2
        record2_out = flip_gt_allele(record2, ref1)
    else
        throw(ArgumentError("$(ref1) not found in record 2"))
    end
    record1_out, record2_out
end

"""
    filter_genotype(record, [genokey=["GT"]])

Filter a VCF record according to `genokey` and output a VCF record with
genotype formats only in `genokey`.
"""
function filter_genotype(
    record::VCF.Record,
    genokey::Vector{T} = ["GT"]
    ) where T <: AbstractString
    format_out = genokey ∩ VCF.format(record)
    gt_out = [Dict(gk => VCF.genotype(record, i, gk) for gk in format_out)
        for i in 1:length(record.genotype)]
    return VCF.Record(record; genotype = gt_out)
end

"""
    filter_genotype(vcffile, outfile, [genokey=["GT"]])

Filter a VCF file according to `genokey` (default `["GT"]`). Output a VCF
file with genotype formats only in `genokey`. Record (markers) with no fields
in `genokey` are skipped.
"""
function filter_genotype(
    vcffile::T,
    outfile::T = "filtered.vcf.gz",
    genokey::Vector{T} = ["GT"]
    ) where T <: AbstractString
    reader = VCF.Reader(openvcf(vcffile, "r"))
    writer = VCF.Writer(openvcf(outfile, "w"), VCF.header(reader))
    for record in reader
        VCF.write(writer, filter_genotype(record, genokey))
    end
    flush(writer); close(reader); close(writer)
    nothing
end

"""
    flip_gt_allele(record[, altlabel])

Flip the REF and ALT allele lable in a VCF record and output a VCF record with
only "GT" genotype format. By default the first label in ALT field is used as
REF label in the flipped record.
"""
function flip_gt_allele(record::VCF.Record, altlabel::String = VCF.alt(record)[1])
    altlabel ∈ VCF.alt(record) || throw(ArgumentError("ALT label not found"))
    VCF.findgenokey(record, "GT") == 0 && throw(ArgumentError("GT format not found"))
    gt_out = [Dict("GT" => VCF.genotype(record, i, "GT"))
        for i in 1:length(record.genotype)]
    record_out = VCF.Record(record; ref = VCF.alt(record)[1],
        alt = [VCF.ref(record)], genotype = gt_out)
    flip_01!(record_out.data, record_out.genotype[1][1][1]:record_out.genotype[end][end][end])
    return record_out
end

"""
    flip_01!(s[, r])
Flip the digits 0 and 1 in a UInt8 vector `s` in range `r`.
"""
function flip_01!(s::Vector{UInt8}, r::Range = 1:length(s))
    for i in r
        if s[i] == 0x30
            s[i] = 0x31
        elseif s[i] == 0x31
            s[i] = 0x30
        end
    end
    s
end

"""
    binomial_proportion_test(x1, x2, y1, y2)

Perform hypothesis test for the equality of two binomial proportions. The
contingency table is structured as
        | Group 1 | Group 2
Count 1 |    x1   |   x2
Count 2 |    y1   |   y2
It tests whether `x1 / (x1 + y1)` is equal to `x2 / (x2 + y2)`. When any of the
counts is less than 5, it resorts to the Fisher's exact test.
"""
function binomial_proportion_test(x1::T, y1::T, x2::T, y2::T) where T <: Integer
    if min(x1, y1, x2, y2) < 5
        # Fisher exact test
        pval = pvalue(FisherExactTest(x1, x2, y1, y2))
    else
        # large sample test
        n1, n2 = x1 + y1, x2 + y2
        p̂1, p̂2, p̂ = x1 / n1, x2 / n2, (x1 + x2) / (n1 + n2)
        ẑ = (p̂1 - p̂2) / √(p̂ * (1 - p̂) * (1 / n1 + 1 / n2))
        pval = ẑ > 0 ? 2ccdf(Normal(), ẑ) : 2cdf(Normal(), ẑ)
    end
    pval
end
