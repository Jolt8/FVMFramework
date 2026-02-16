General notes on using a tracer:
    - (This is by far my most important question) I'm wondering what kind of workflow will be preferrable whenever the tracer methods are finished. Would it be better to have the user only define the physics and not the properties of each region, run the tracer and then request the user to fill out each property that is required to solve given the physics they defined. In contrast, the workflow that many other CFD solvers use is one where the solved-for variables are defined first, vel_x, vel_y, vel_z, rho, pressure, temp, mass_fractions, etc. and then the correct physics are used to describe how these parameters change with time.
    - One of the disadvantages of the tracer would be the difficulty in dealing with if statements, for example let's say we were simulating laminar flow through a pipe and we wanted to use physics that were less expensive than full darcy-weisbach or navier-stokes. We could use friction factor and have an if statement like:
        if Re < 2300 
            u.friction_factor[cell_id] = 64
        else
            u.friction_factor[cell_id] = some_function(du, u, cell_id)
        end
    This is problematic because if the some_function referenced a different cache variable that would not be encountered in the Re < 2300 case, the tracer would fail. I'm pretty sure that TracerLocalSparsityDetector can already handle this, so we'd have to find out how they deal with this (hopefully it's not too difficult).
    - Something to note, using a tracer doesn't make it so that we can't use generic references in physics functions, we still can achieve this through merging u_state, u_caches, and u_properties vectors passed into the function. 
        - However, for stuff like rho where it can be a state, cache or fixed property, I feel like reconstructing the u_merged from these input variables would become quite complex, probably less complex than using the tracer, but I feel like the small amount of extra development time developing the tracer for potentially very big benefits would probably be worth it. 
            - This becomes especially true when we eventually implement the optimzie!() method to create a pointer to a p[idx] for optimization. 
    - Another thing I wonder is if the Tracer actually has to do any math within the function. I'm not sure how much benefit we net from doing math with the Tracer, but actually returning real values has been a struggle, although sample() helps a ton. Also, I really hate my current method because of the similar_sets and subsets definitions which defeat the point of the tracer in the first place. I'm guessing I'm going to have to spend a lot of time creating custom methods for returning data when a ComponentVector is accessed, but it's going to be a struggle. 
    - I think we're going to have to create custom methods for all the mathematical operators (not fun), to actually not have to do any math within the function. 
    - I think one thing that's going to be pretty hard to tackle when designing the tracer is getting how variables are related to each other right. For example, for molar_concentrations, it's almost always going to be found from mass_fractions while mass_fractions can also be found from molar_concentrations, so I feel like these two variables are obviously linked, but I don't know if making this distinguishment would actually be of any use. 
    - Oh, another thing that I realized is that our current method for doing capacities is actually pretty confusing for solver organization. For example with heat flux, right now we keep track of du.temp in units (J/s) and then divide by the heat capacity of the cell (J/K) to get the actual du.temp (K/s). What we could do is store du.heat as a cache and then divide by the heat capacity of each cell to get du.temp. While I don't think we should do this, it's just a little inconsistent with how every single variable's name usually tries to be accurate to the actual physical quantity it's defining. Another way we could achieve this is by just dividing the flux by the capacity and then taking it away from du.field.
    - One of the main goals of this framework is to make it as customizable as possible to allow anyone with a basic understanding of Julia to use its methods. I'm always worried about how "magic" the Tracer will feel compared to how grounded the rest of this framework is in pretty basic Julia. I mean, maybe this is just something users will have to accept just like how I've never really tried to understand how TracerLocalSparsityDetector works and just accepted it as "magic". 
    - This is probably the biggest question: will the tracer build a merged u vector index map or will it generate the code for the f function itself
        - DO NOT DO CODE GENERATION!!!, I just realized that option is very close to MTK and I hate working with MTK

Notes on actually writing it
    - State Variables:
        YES:
            - Write to du.field
            - Read from u.field
        NO:
            - write to u.field
        
        Additional Rules:
            - Never read again after being written to

        This would also include du.integrator_error, but du.integrator_error would always be set to 0.0 originally
    - Caches
        YES/NO:
            - Write to du.field
            - Write to u.field
        YES:
            - Read from u.field
        Additional Rules:
            - Always read again after being written to 
    - Fixed Parameters
        YES:
            - Read from u.field
        NO:
            - Write to du.field
            - write to u.fiel
    - Global parameters
        - stuff like optimize!(u.field)
        pointer to p[idx]
    - dealing with loops
        - We need the user to provide minimal context to the tracer, however even these could be requested and written to a permanent JSON file when the tracer is first run
        Ex:
            context(
                species = [:methanol, :water, :co2]
                phases = [:gas, :liquid]
            )
        - Example: trace(f_closure)
            at the end it would say something like "oh, there's a few loops at lines X, Y, Z and I have no idea what we're looping through. Then it could write what the first part of the loop looks like to the json 
                Ex:
                u.mass_fractions needs species names (e.g. [:methanol, :water])
                u.reactions needs reaction names (e.g. [:MSR, :WGS])
            and then the user could provide simply [:methanol, :water, :co2] 
            
            Furthermore, instead of using a json to actually use variables you could define species_names = [:methanol, :water, :co2] at the top and then reference species_names for all other loops

    - Perhaps to distinguish between some types of variables a macro like @flux could be used like: 
        @flux du.temp[cell] = du.heat[cell] / u.rho[cell] * u.cp[cell] * u.vol[cell]
    to distinguish between certain types of variables
    - Another method would be to just provide a list of possible Fixed Variables and Caches and then let the user decide which variables they use. 
        - Caches would only show up in this list if they're written at du.field[cell]

    Other features to potentially add in the future:
        - Mermaid diagram of dependencies like how mass_fractions relates to reaction_rate, heat_source, and temp
    
    Steps to implement
        - Struct that overloads getproperty, setproperty!, getindex, setindex! and propertynames. 
        - Math operator overloads that return valid-ish numbers 
        - Variable classification that accumulatees read write logs and sorts variables into states, fixed parameters, caches and other parameters
        - Loop discovery on propertynames and a fake set of keys
        - Loop context and template script generator 