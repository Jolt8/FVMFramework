Collection of Notes for the geometry optimization solver I'm working on

NOTE TO FUTURE SELF ON PURE MUTATION:
    - I tried to pass a struct that acted as a cache for all the variables to prevent GC when recreating the vectors in rebuild_fvm_geometry 
    - I also tried to pass a cache of node coordinates to mutate every iteration 
    - Anyways, that didn't work because of this error message when doing the Optimization.solve(...)
        MethodError: no method matching Float64(::ForwardDiff.Dual{ForwardDiff.Tag{SciMLBase.ParamJacobianWrapper{…}, Float64}, Float64, 12})
        The type Float64 exists, but no method is defined for this combination of argument types when trying to construct it.
    - Based on this, it seems that my implementation located in (Geometry Optimization Only Mutation.jl) was enforcing Float64 within the struct because I can't switch between Float64 and the weird ForwardDiff.Dual Float64 type (or at least I can't think of a way how). 
    - Thus, it seems like I have to stick with the functional versions of apply_ffd_motion and rebuild_fvm_geometry 

NOTE: For using caches
    - We need to use PreAllocationTools.jl when the GC from using caches updated in the solver itself becomes a problem 

NOTE: 
    - We might have to declare the unchanging values passed into the closure of the function as const before passing them in to prevent GC
    - ex. 
        - const cell_neighbor_map = cell_neighbor_map 
        - const const_neighbor_map_respective_node_ids = neighbor_map_respective_node_ids


NOTE: What I have learned after a lot of debugging 
    - We should probably just give up on fixing zygote and ReverseDiff
    - Instead, we should focus on getting Enzyme to work by preventing any dynamic dispatch because it supports array mutation
    - While we try to get that to work, we should just use AutoForwardDiff as the deform box 