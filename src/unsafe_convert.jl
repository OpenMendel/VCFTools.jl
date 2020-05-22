#### This file is adapted from https://github.com/jiahao/VCF.jl


#Eat the second argument for non-GZip IOstreams
Base.position(b::Base.GenericIOBuffer, ::Bool) = position(b)

"Initial guess for how many nonzero entries will be read in"
const INITIAL_SIZE_GUESS = 1_000_000

"Fudge factor used to reallocate storage vectors"
const SLACK_FACTOR = 1.25

"Update progress bar each time this many lines are read"
const UPDATE_LINE_INTERVAL = 500

"""
Read VCF files into a [0,1,2]-genotype matrix

Implementation note:

    This is a very simplistic VCF file reader that skips a lot of the data
    stored in the file.

    All header lines are ignored.

    Only field 2 (the chromosomal position) and fields 10+ (the variants) are
    read in.

    The phases of the variants are ignored.

Input:

    filename: Name of VCF file. Can be `.gz` compressed

Output:

    Is: Row indexes of variants
    Js: Column indexes of variants
    Vs: Number of variants at position I for sample J

Side effects:

    Also displays progress meter

Reference:

    Based on the VCF Version 4.2 Specification at

    http://samtools.github.io/hts-specs/VCFv4.2.pdf
"""
function unsafe_convert_gt(vcffile)
    Is = Int64[]
    Js = Int64[]
    Vs = Int8[]
    sizehint!(Is, INITIAL_SIZE_GUESS)
    sizehint!(Js, INITIAL_SIZE_GUESS)
    sizehint!(Vs, INITIAL_SIZE_GUESS)

    csize = filesize(vcffile)
    p = Progress(csize, 2, "", 50) #Progress bar

    mode = :startline
    lineno = pos = fieldidx = skiplines = 0
    thisentry = 0
    rowid = colid = 1
    tmpbuf = IOBuffer()
    n = 0

    #Now actually try to read file
    if csize==0
        @warn("File $vcffile is empty")
        @goto done
    end

    stream = if endswith(vcffile, ".gz")
        GzipDecompressorStream(open(vcffile))
    else
        IOBuffer(Mmap.mmap(vcffile))
    end

    #Iterate through file
    while true
        c = read(stream, UInt8)
        # print(String([c]))
        n += 1
        if mode == :startline
            lineno += 1
            fieldidx = 1
            mode = c == 0x23 ? :skipline : :seekpos #Skip comment (0x23 = '#' in ASCII)

        elseif mode == :skipline
            if c == 0x0a #'\n'
                skiplines += 1
                mode = :startline
            end

        elseif mode == :seekpos
            startline = false

            if c == 0x09 #'\t'
                fieldidx += 1

                if fieldidx==2
                    mode = :readpos
                    seek(tmpbuf, 0)
                end
            end

        elseif mode == :readpos
            if c == 0x09 #'\t'
                pos = parse(Int64, take!(tmpbuf))
                fieldidx = 3
                mode = :seekdata
            else
                write(tmpbuf, c)
            end

        elseif mode == :seekdata
            if c == 0x09 #'\t'
                fieldidx += 1
                if fieldidx == 10
                    mode = :readdata
                    colid = thisentry = 0
                end
            end

        elseif mode == :readdata
            if c == 0x30 && position(tmpbuf)==0 # '0'
                continue

            elseif 0x30 ≤ c ≤ 0x39 || c == 0x2e # '0' ≤ c ≤ '9' || c == '.'
                #Horrible hack to force tmpbuf pointer forward without storing a value
                tmpbuf.ptr += 1
                #write(tmpbuf, c)

            elseif c == 0x0a #'\n'
                rowid += 1
                mode = :savedata

            elseif c == 0x09 || c == 0x0d #c == '\t' || c == '\r')
                #Reset IO buffer
                if position(tmpbuf) > 0
                    thisentry += 1
                    seek(tmpbuf, 0)
                end
                colid += 1

                mode = :savedata

            elseif c == 0x7c || c == 0x2f # '|' or '/'
                #Reset IO buffer
                if position(tmpbuf) > 0
                    thisentry += 1
                    seek(tmpbuf, 0)
                end

            elseif c == 0x3a #':' #Extra genome-level data I don't understand
                mode = :skipfield

            else
                error("Unknown c = '$(String([c]))' on line $lineno, file offset $n")
            end

        elseif mode == :skipfield
            if c == 0x09 #'\t'
                mode = :readdata
            elseif c == 0x0a #'\n'
                mode = :startline
            end
        end

        if mode == :savedata
            if lineno % UPDATE_LINE_INTERVAL == 0 #Update progress bar
                coffset = position(stream)
                update!(p, coffset)
            end

            if thisentry != 0#position(tmpbuf) > 0
                push!(Is, pos)
                push!(Js, colid)
                push!(Vs, thisentry) #parse(Int8, take!(tmpbuf)))
                thisentry = 0

                if length(Vs) == INITIAL_SIZE_GUESS
                    #More data than initially guessed
                    #Estimate memory consumption and resize in-memory storage
                    coffset = position(stream)
                    estnnz = ceil(Int, length(Vs)*csize/coffset)
                    # @info("Number of nonzero entries read so far:", length(Vs))
                    # @info("Estimated number of nonzero entries:", estnnz)
                    estnnz = ceil(Int, SLACK_FACTOR*estnnz)
                    sizehint!(Is, estnnz)
                    sizehint!(Js, estnnz)
                    sizehint!(Vs, estnnz)
                end
            end

            mode = c == 0x0a ? :startline : :readdata #'\n' = 0x0a
        end

        eof(stream) && break
    end

    @label done

    finish!(p)
    println("Actual number of header lines:", skiplines)
    println("Actual number of record lines:", lineno-skiplines)
    println("Actual number of nonzeros:", length(Vs))

    if length(Vs) > 0
        nrows = maximum(Is)
        ncols = maximum(Js)
        println("Number of record rows: ", rowid)
        println("Number of logical rows: ", nrows)
        println("Number of logical cols: ", ncols)
        println("Average density: ", length(Vs)/(Int64(nrows)*Int64(ncols)))
    end

    close(stream)

    Is, Js, Vs
end

function Base.parse(::Type{T}, s::AbstractArray{UInt8}) where T<:Integer
    n = zero(T)
    @fastmath for c in s
        n::T = 10n + (c - T(0x30))
        if n<0 || !(0x30 ≤ c ≤ 0x39 ) #!('0'<=c<='9')
            error(string(Char(c),':',s,':',n))
        end
    end
    n
end

"""
Naive convert function that splits a vcf file line by line, assuming data only contains GT field. 
"""
function unsafe_convert_gt2(vcffile::String)
    stream = if endswith(vcffile, ".gz")
        GzipDecompressorStream(open(vcffile))
    else
        IOBuffer(Mmap.mmap(vcffile))
    end

    p, n = nrecords(vcffile), nsamples(vcffile)
    out = Matrix{Union{Missing, UInt8}}(undef, p, n)
    sample_data = Vector{String}(undef, n + 9)
    GT = Vector{Union{Missing, UInt8}}(undef, n)

    l = 1
    for line in eachline(stream)
        if !startswith(line, "#")
            # split by sample
            sample_data .= split(line, "\t")

            # split by field
            # split_words = map(x -> split(x, ":"), sample_data[10:end])

            # Assume GT field comes before all other field
            map!(parse_gt, GT, @view(sample_data[10:end]))
            copyto!(@view(out[l, :]), GT)
            l += 1
        end
    end

    return out
end

"""
Same as unsafe_convert_gt2 but tries to read lines in bulk and process them in parallel
"""
function unsafe_convert_gt3(vcffile::String)
    vcf_buffer_records = 1024

    stream = if endswith(vcffile, ".gz")
        GzipDecompressorStream(open(vcffile))
    else
        IOBuffer(Mmap.mmap(vcffile))
    end

    p, n = nrecords(vcffile), nsamples(vcffile)
    out = Matrix{Union{Missing, UInt8}}(undef, p, n)
    sample_data = [Vector{String}(undef, n + 9) for i in 1:Threads.nthreads()]
    GT = [Vector{Union{Missing, UInt8}}(undef, n) for i in 1:Threads.nthreads()]

    chunks = floor(Int, p / vcf_buffer_records)
    remaining_records = p - vcf_buffer_records * chunks
    records = Vector{String}(undef, vcf_buffer_records)

    # process header
    headerlines = 0
    while Base.peek(stream) == 0x23 # starts with '#'
        headerlines += 1
        records[headerlines] = readline(stream)
        # TODO process header line
    end

    # process chunk by chunk
    for c in 1:chunks
        # load a bunch of records
        for i in 1:vcf_buffer_records
            records[i] = readline(stream)
        end

        # process records in parallel
        l = (c - 1) * vcf_buffer_records # beginning of current chunk
        Threads.@threads for i in 1:length(records)
            id = Threads.threadid()

            # split by sample
            sample_data[id] .= split(records[i], "\t")

            # split by field
            # split_words = map(x -> split(x, ":"), sample_data[10:end])

            # Assume GT field comes before all other field
            map!(parse_gt, GT[id], @view(sample_data[id][10:end]))
            copyto!(@view(out[l + i, :]), GT[id])
        end
    end

    # process remaining records
    resize!(records, remaining_records)
    for i in 1:remaining_records
        records[i] = readline(stream)
    end
    l = chunks * vcf_buffer_records # beginning of last chunk
    Threads.@threads for i in 1:length(records)
        id = Threads.threadid()

        # split by sample
        sample_data[id] .= split(records[i], "\t")

        # split by field
        # split_words = map(x -> split(x, ":"), sample_data[10:end])

        # Assume GT field comes before all other field
        map!(parse_gt, GT[id], @view(sample_data[id][10:end]))
        copyto!(@view(out[l + i, :]), GT[id])
    end

    return out
end

function parse_gt(GT::AbstractString)
    if ishomozygous_zero(GT)
        return 0x00
    elseif isheterozygous(GT)
        return 0x01
    elseif ishomozygous_one(GT)
        return 0x02
    else
        return missing
    end
end
ishomozygous_zero(vcfentry::AbstractString) = vcfentry == "0/0" || vcfentry == "0|0"
isheterozygous(vcfentry::AbstractString) = vcfentry == "1/0" || vcfentry == "0/1" || vcfentry == "0|1" || vcfentry == "1|0"
ishomozygous_one(vcfentry::AbstractString) = vcfentry == "1/1" || vcfentry == "1|1"
Base.ismissing(vcfentry::AbstractString) = vcfentry == "./." || vcfentry == ".|."
