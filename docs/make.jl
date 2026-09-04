using Documenter
using FVMFramework # Replace with your actual package name

makedocs(
    sitename = "FVMFramework",
    modules = [FVMFramework],
    authors = "Will Martin <willemartin@icloud.com>",
    pages = [
        "FVMFramework.jl: A Framework for solving Finite Volume Method Problems in Julia" => "index.md",
        "Setting up a Simulation" => "quickstart/sim_config_setup.md"
    ]
)

deploydocs(
    repo = "https://github.com/Jolt8/FVMFramework.jl",
)