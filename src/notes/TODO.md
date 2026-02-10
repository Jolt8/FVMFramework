TODO:
    - Implement GPU acceleration later

    - Split internal physics, sources, and capacity parameters instead of one big MethanolReformerPhysics
        - This is a maybe, it would be more modular but would require the redefnition of properties and would not be great for connections
        - I'm going to ignore this one, but it might be implemented later if there's a need for it

    - Implement FEM eventually for things like stress analysis

    - Implement a parameter optimization pipeline for fitting (method 1) and finding (method 2) best parameters (probably two separate files/folders/methods)

    - Implement VLE (Vapor-Liquid-Equilibrium)

    - Implement Navier-Stokes

    - Create methods for interfacing with other FVM or DG solvers like Trixi.jl or Oceanigans.jl

    - Create methods that run multiple physics functions inside them like fluid_fluid_flux!(...) that does continuity_and_momentum_darcy, diffusion_temp_exchange!, species_advection!, enthalpy_advection!, etc.
        - Done! I'm still keeping this one open because we may or may not revoke this
    
    - Create methods to more easily diagnose problems through logging. Graphing a certain variable over time would be very useful
        - this really needs to be worked on 
        - Just some brainstorming:
            - It would be nice if every single variable was tracked
            - This could probably be achieevd through a debug version of each flux and internal physics function where if problems emerge such as values like NaN, Inf, -Inf, 0.0, etc.

    - Improve the way mass fractions are accessed
        - This could probably be achieved through another ComponentVector where it would be something like u.mass_fractions.methanol[cell_id] or u.mass_fractions[cell_id].methanol

    - Add a method to finish_fvm_config(...) to take all cells that are not under a cell set and add them to a default set with either a user defined default region or a default default region with no physics or variables 
        - We would also have to exclude it from the connections map

    - Create a method for doing distributing computing by recruiting my home PC as a worker to prevent having to do simulation on my laptop in class

    - Implement a method within the connecitons constructor that always puts a certain type of connection first. For example, always putting fluids first instead of solids in fluid_solid connections to make sure that the rho of idx_a can be found via cell_rho_ideal!() and the rho of idx_b can be accessed by phys[phys_b].rho. This would prevent dynamic dispatch which is nice and is significantly more scalable as long as the user remembers what order they put the priority list in.
        - the priority list would look something like: AbstractPhysics[MethanolReformingArea, FluidPhysics, SolidPhysics, etc...]