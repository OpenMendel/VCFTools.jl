using Documenter, VCFTools

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = Documenter.HTML(),
    sitename = "VCFTools",
    modules = [VCFTools],
    clean = true,
    pages = [
        "Home" => "index.md",
        "VCF Summary" => "man/summaryinfo.md",
        "Filter" => "man/filter.md",
        "Convert" => "man/convert.md",
        "Match markers in two VCF" => "man/conformgt.md",
        "API"  => "man/api.md"
    ]
)

deploydocs(
    repo   = "github.com/biona001/VCFTools.jl.git",
    target = "build"
)
