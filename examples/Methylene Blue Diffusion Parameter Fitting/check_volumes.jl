
using Pkg; Pkg.activate("c:/Users/wille/Desktop/FVMFramework")
using FVMFramework
using Ferrite
using FerriteGmsh

mesh_path = "c:/Users/wille/Desktop/FVMFramework/examples/Methylene Blue Diffusion Parameter Fitting/cone_end/dialysis_tubing_cone_output.msh"
grid = togrid(mesh_path)

u_proto = (
    mass_fractions = (methylene_blue = zeros(length(grid.cells)), water = zeros(length(grid.cells))),
)
config = create_fvm_config(grid, u_proto)

let
    tube_vol = 0.0
    for cell_id in grid.cellsets["dialysis_tubing_interior"]
        tube_vol += config.geo.cell_volumes[cell_id]
    end

    reservoir_vol = 0.0
    for cell_id in grid.cellsets["surrounding_fluid"]
        reservoir_vol += config.geo.cell_volumes[cell_id]
    end

    println("Tube Volume: ", tube_vol * 1e6, " ml")
    println("Reservoir Volume: ", reservoir_vol * 1e6, " ml")
    println("Total Volume: ", (tube_vol + reservoir_vol) * 1e6, " ml")
    println("Ratio (Tube / Total): ", tube_vol / (tube_vol + reservoir_vol))
    println("Equilibrium Methylene Blue Concentration (assuming 0.0004 initial in tube): ", 0.0004 * tube_vol / (tube_vol + reservoir_vol))
end
