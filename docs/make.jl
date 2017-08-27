using Documenter, VCFTools

ENV["DOCUMENTER_DEBUG"] = "true"
makedocs()
deploydocs(
  deps   = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
  repo   = "github.com:OpenMendel/VCFTools.jl.git",
  julia  = "0.6",
  osname = "osx"
  )

# ## for local build
# makedocs(
# #      doctest   = false,
#       format    = :html,
#       clean     = true,
#       sitename  = "VCFTools.jl",
#       modules   = [VCFTools],
#       pages     = [
#           "Home"  => "index.md",
# 	  "VCF summary" => "summaryinfo.md",
#      "Match markers in two VCF" => "conformgt.md",
#      "API" => "api.md"
#       ]
# )
#
# deploydocs(
#       repo    = "github.com:OpenMendel/VCFTools.jl.git",
#       target  = "build",
#       julia   = "0.6",
#       deps    = nothing,
#       make    = nothing
# )
