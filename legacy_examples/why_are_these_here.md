Why are these files here any why are they so long?
    - These files represent previous workflows, projects, attempts and other stuff previously done in my other repository
    - They mainly represent methods that have yet to be implemented in this repository
    - My plan going forward is to eventually implement the methods demonstrated in these separate files into the main branch through a more modular interface
    
To give more info on the purpose of each file in no particular order:
    - Check Mesh Function For Later.jl
        - This one is a helper function to check the quality of the mesh similar to what OpenFOAM does while meshing
        - Right now, the approach it takes to checking the mesh is pretty naive so it will probably need a rewrite in the future
        - For now, it's a good example of something to implement in the future

    - Ea and A Parameter Estimation All Temps FR.jl
        - This project was part of my Chemistry IA to determine the activation energy and pre-exponential factor of the esterification of acetic acid with ethanol. The project was highly successful despite the questionable quality of the data. 
        - In the future, I'm hoping to automate many of the methods used in this to make parameter fitting for this repo very easy as that's its primary purpose
        - This includes methods for 
            - Using XLSX to pull data for the simulation
            - Fixing a variable in the solver to measured values in the real world
                - In that file, the temperature of the reaction was measured and in the simulation the temperature was fixed to that value
            - And a lot of other minor things including setting up the optimization itself
    
    - Generalized DG.jl
        - This was a partially finished version of implementing Discontinuous Galerkin (DG) in Julia, although I never got far because I just couldn't wrap my head around DG
        - There was plans to automatically switch between DG for heat transfer, very linear problems/areas and stuff like that to FVM automatically for VLE stuff, shocks, and heavy discontinuities, although I never made it that far
        - In the future, I might try to implement DG again if it can be rather simply merged with the current methods, but I doubt it
        - Instead, an interface with Trixi.jl will probably be created 

    - Generalized FVM Best with Adjoint Shape Optimization.jl
        - This was my first ever attempt at geometry optimizaiton that I quickly moved away from once I figured out how to make AD pass through a function to move the mesh with a deform box
        - That being said, it will probably be a lot easier to get Reverse Diff working with it, so it might be useful in the future
        - The way this one works is that each cell has an "alpha" value that describes how "solid" it is. Basically, you try to force the solver to make alpha either 1 or 0 and not 0.5 while minimizing an objective function.
        - However, compared to actual geometry optimziation, getting the solver to push the alpha of each cell to 1 or 0 is very difficult and often more difficult than actually trying to optimize what you're trying to optimize

    - Geometry Optimization Navier Stokes.jl
        - Includes methods to perform Navier Stokes in Julia
        - This was not finished, but it was very close and probably needs a little bit of polishing to bring it to a suitable state for this repository

    - Geometry Optimization New Features.jl
        - This is the most up-to-date version of doing geometry optimization for FVM in Julia
        - It's still not compatible with any reverse differentation codebase, but I'm waiting for Enzyme to be updated to Julia 1.12 before I do this
            - To summarize why I don't like ReverseDiff.jl and Zygote.jl:
                - ReverseDiff is generally less performant than Enzyme, so when it comes to putting in the effort to make any type of reverse differentiation work, immediately spending effort to move to Enzyme seems more worth it 
                - Zygote is completely off the table because I love mutation and Zygote hates mutation. Also, Enzyme is becoming more and more popular in the future
                - If this repository is ever seriously used by Chemical Engineers or myself in the future (which I'm planning to do), I will probably pay someone with a lot of knowledge of Enzyme to get Enzyme working with this repository, but that's a future problem. 
        - This script works by creating a boundary box that captures the the nodes of the solver's mesh. Then, the xyz position of each of the boundary box's nodes is adjusted using AD and SciMLSensitivity to minimize an objective/loss function.
    
    - Geometry Optimization Notes.txt
        - This contains some of the lessons learned while doing geometry optimization in Julia, so I'll just summarize it
        - Structs for organizing caches doesn't work (unless PreallocationTools.jl can be present inside a struct) because they are so strictly typed
        - Use PreallocationsTools.jl for mutating caches within the solver function to avoid GC
        - We might need to pass in the fixed values passed into the f_closure as const before passing them in
