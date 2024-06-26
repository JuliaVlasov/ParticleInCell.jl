push!(LOAD_PATH, "../src/")
using Documenter
using Plots
using ParticleInCell

ENV["GKSwstype"] = "100"

makedocs(
    modules = [ParticleInCell],
    sitename = "ParticleInCell.jl",
    authors = "Pierre Navaro",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliavlasov.github.io/ParticleInCell.jl",
        assets = String[],
    ),
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Maxwell solver" => "maxwell.md", #"Landau damping" => "landau_damping.md",
        "Two-stream instability" => "tsi.md",
        "Vlasov-Maxwell 2D" => "vlasov-maxwell.md",
        "Vlasov-Poisson 2D" => "vlasov-poisson.md",
        "Weighting schemes" => "weighting_schemes.md",
        "Types" => "types.md",
        "Functions" => "functions.md",
        "Contents" => "contents.md",
    ],
)

deploydocs(;
    branch = "gh-pages",
    devbranch = "master",
    repo = "github.com/JuliaVlasov/ParticleInCell.jl",
)
