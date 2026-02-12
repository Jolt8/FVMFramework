Merging component Vectors
    - We can merge u and caches and then only pass u to each physics function and then I'll always be able to call u.rho even if rho is a state variable (inputted u value) or a cached variable the only issue with this approach is that it's a little less clear which parameters are state variables (u_true) and caches

Geomery data in u_merged?
    - Perhaps even the geometry data could be merged with u this would be great because we could make geometry optimizaiton functions mutating of geometry data caches if we did this, we couldn't store cell_volumes or cell centroids within a component vector that represented a single cell

Another great idea for optimization:
    - We could create a method that creates pointers to the p in this function for optimizaiton. For example, u.reactions.reforming_reactions.kf_Ea = optimized!(p_proto) then, while running the funciton, the value that represents u.reactions.reforming_reactions.kf_Ea could be replaced by p.reaction.reforming_reactions.kf_Ea 
    - the p.reactions.reforming_reactions.kf_Ea would be passed in as a flat vector and then converted into a component array wit ComponentVector(p, p_merged_axes). The input p values would be constructed over the initial setup through optimized!(p_proto)
    - I'm sure we'll figure something out
    - The bottom line is that a real value should be able to be replaced by simply encasing the original value with optimized!() if we encased the original value with optimized!(), we could use that as our initial p_guess
    - I think the way we could do all of this is by intercepting pointers 
    - I think the way we could do this is by keeping track of a tracked_u and tracked_du and then appending the field that is trying to be accessed as the fake function runs if we were to do this, to actually get the physics to run, we would have to supply a value that can almost always be processed, like probably 1.0 or random numbers and then whenever a field is accessed, we intercept the field attempting to be accessed, 

Julia package complexity and direction this framework should be heading
    - Sometimes I wonder if the increasingly underlying complexity of this framework is worth it one of the things that I don't really like about some julia libraries is when they have an overwhelming preference towards simplicity and representing physics as pure math while sacrificing customizability
       - A great example of this is Gridap.jl which I remember seeing a while back and it just screamed "not friendly to user understanding or user defined custom methods" Despite probably being very flexible
    - That being said, I feel like there a point where a package reaches a critical mass like Enzyme.jl or ForwardDiff.jl where it covers so many edge cases that the lack of being not so friendly to users doesn't really matter
    - I would say that for doing the finite volume method users should be able to understand how the framework works
    - I would also like the distinction that a framework feeling flexibility != transparent (kinda obvious, but whatever)
    - I would define flexibility as the framework's ability do do various things easil, while customizability refers to how transparent the lowest level methods of that package are. That being said, when flexible packages are done well like ForwardDiff or Enzyme, it's really nice

Distinguishing between u_state, u_fixed, u_cached to create u_merged
    - Although this is a very far off goal, what we could do is run a fake version of methanol_reformer_op_different_connections once and then based on the data collected we see which variables are state variables, optimized_variables, fixed variables, cached/intermediate variables etc.. 
    - This could even be great because it could help us decide whether or not caching certain values is worth it or not
    - And then we construct u_state, u_fixed, u_cached and the user only has to interact with defining variables for u_fixed 
    - This would be great because we wouldn't have to updates properties within each region function like mw_avg!() or rho_ideal!()
    - This would also be great for debugging (probably) and would ensure that cached variables are always updated correctly
    - One thing I'm a little skeptical of is whether or now thinking about physics funcitons while not even being entirely sure what your u_state variables are would be too confusing
    - Like right now, this makes a lot of sense for organization, but it might just be that I have the gift of knowing what variables belong to which before developing this method, and then it would become pretty confusing when working with a new problems where I don't know which are u_state variables, cached, 
