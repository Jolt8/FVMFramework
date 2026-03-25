TODO:
    - For second-order spatial discretization, I've been recommended to use Venkatakrishnan limiters, Van Albada Limiters, or maybe even Learned Neural Flux Limiters
    
    - I should probably implement second order meshes at some point. In my Biology IA that contained a lot of cylinders, I had to make a very fine mesh to capture the curve of the cylinders. This was especially bad for very thin cylinders.

    - The other thing that would be great to implement is a method to handle moving meshes, I'm sure this would be kinda challenging, but I've already implement deform box mesh moving in the past, so I don't think that Arbitrary Langrangian-Eulerian would be too far out of reach

    - We have to choose between implementing second-order spatial discretization, second-order meshes, moving meshes, and adaptive mesh refinement as the next project

    - Implement a method to check units right after finish_fvm_config 
        - This would also mean that unitful values would have to be stripped when finish_fvm_config is called
        - Note that you can turn any quantity into base SI units by doing upreferred(1.0u"cm")
        - You can also set your preferred units with Unitful.preferunits(u"s", u"m", u"kg", u"K")
        - This has kinda always bugged me, but there's no way to simplify units like in 1.0u"kg*m^2/s^2"|> u"J", but I guess that's for the best because it might convert something like Pa * s to Poise (who the fuck uses Poise)
        - Also, did you know that m_fix = Unitful.FixedUnits(u"m") exists where you can return errors when something like 1.0m_fix + 1.0u"cm" is called, that would probably be useful in checking for unit consistency
    
    - Implement GPU acceleration later

    - Split internal physics, sources, and capacity parameters instead of all of them being lumped in one big region
        - This is a maybe, it would be more modular but would require the redefnition of properties and would not be great for connections
        - I'm going to ignore this one, but it might be implemented later if there's a need for it

    - Implement FEM eventually for things like stress analysis (long way away from this one)

    - Implement a parameter optimization pipeline for fitting (method 1) and finding (method 2) best parameters (probably two separate files/folders/methods)
        - Partially done, I don't think it's going to be worth it to try to create methods that try to cover the many different way optimization is done
        - instead, I think it's much better we just show how users can do something and then let them implement their own methods on a case-by-case basis 

    - Implement VLE (Vapor-Liquid-Equilibrium)

    - Implement Navier-Stokes

    - Create methods for interfacing with other FVM or DG solvers like Trixi.jl or Oceanigans.jl

    - Create methods to more easily diagnose problems through logging. Graphing a certain variable over time would be very useful
        - this really needs to be worked on 
        - Just some brainstorming:
            - It would be nice if every single variable was tracked
            - This could probably be achieevd through a debug version of each flux and internal physics function where if problems emerge such as values like NaN, Inf, -Inf, 0.0, etc.
        - We should probably look into the package JuliaDebugger that provides breakpoints that would be insanely useful

    - Create a method for doing distributing computing by recruiting my home PC as a worker to prevent having to do simulation on my laptop in class

    - Implement MPI