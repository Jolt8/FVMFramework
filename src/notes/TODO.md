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

    - Create methods to more easily diagnose problems through logging. Graphing a certain variable over time would be very useful
        - this really needs to be worked on 

    - Create a dedicated method for distributing a BC over a set of cell volumes or set of cell areas 

    - Improve access of mass fractions 

    - Add a method to finish_fvm_config(...) to take all cells that are not under a cell set and add them to a default set with either a user defined default region or a default default region with no physics or variables 
        - We would also have to exclude it from the connections map