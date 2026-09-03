using Documenter
using FVMFramework # Replace with your actual package name

makedocs(
    sitename = "FVMFramework",
    modules = [FVMFramework],
    authors = "Will Martin <willemartin@icloud.com>",
    pages = [
        "Home" => "index.md",
    ]
)

deploydocs(
    repo = "https://github.com/Jolt8/FVMFramework.jl",
)