TODO:
    - Create a universal vector that passes 
    
    - Implement Multithreading
        - Implement GPU acceleration later

    - Split internal physics, sources, and capacity parameters instead of one big MethanolReformerPhysics

    - Implement FEM eventually

    - Implement a parameter optimization pipeline for fitting and finding best parameters (probably two separate files/folders/methods)

    - Implement PreAllocationTools.jl

    - Improve BC definition to allow stuff like a constant mass flux in of a certain mixture

    - Implement VLE (Vapor-Liquid-Equilibrium)

    - Implement Navier-Stokes

    - Create methods for interfacing with other FVM or DG solvers like Trixi.jl or Oceanigans.jl

    - Create methods that run multiple physics functions inside them like fluid_fluid_flux!(...) that does continuity_and_momentum_darcy, diffusion_temp_exchange!, species_advection!, enthalpy_advection!, etc.