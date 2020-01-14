using Documenter, VCFTools

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = Documenter.HTML(),
    sitename = "VCFTools",
    modules = [VCFTools],
    clean = true,
    pages = [
        "Home" => "index.md",
        "VCF Summary" => "src/summaryinfo.md",
        "Filter" => "src/filter.md",
        "Convert" => "src/convert.md",
        "Match markers in two VCF" => "src/conformgt.md",
        "API"  => "src/api.md"
    ]
)

deploydocs(
    repo   = "github.com/biona001/VCFTools.jl.git",
    target = "build"
)
