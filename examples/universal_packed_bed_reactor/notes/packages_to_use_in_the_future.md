Out of the packages that I've seen so far, the three that seem by far the most useful are DiffEqFlux.jl, Surrogates.jl, and SymbolicRegression.jl 

The reason why I like these three so much is that they align very well with the goals of this project

For context, the discipline that I will be going into in the future is chemical engineering. So far, I've determined that parameter estimation is *probably* the most valuable thing that I can use Julia for, although I might be proved wrong since I don't know if there's something that I haven't encountered yet in my current engineering journey that is actually super valuable. Like for example, I've never had a case where I wanted to make a microcontroller maximize hydrogen yield for example. However I do have experience with how valuable parameter estimation because I've found that academic sources, even stuff that's widely cited as very accurate such as the Peppley Amphlett, Mann methanol reforming model to not be trustable. Thus, I've determined that every time I will design a chemical reactor, I will probably have to come up with my own reaction kinetics model and find the two (or more) arrhenius parameters for that system validated against data in my own pilot plant where the goodness of fit can be validated directly. If I get good at this, it has the advantage that I will be able to find the reaction kinetics (or basically any kind of physics ever like finding the overall heat transfer coefficient for a packed bed reactor) for basically any reaction or catalyst extremely quickly without having to rely on academic sources. Also, the advantage of my this approach is that since I can model the solution transiently, any variations in basically any variable like temperature that would have previously had to been carefully controlled to allow a steady state assumption to work. In fact, my current project is going to be finding the 2 arrhenius parameters of the 3 main reaction in methanol reforming using just 5 thermocouples and a reactor that I built for under 500$ (no specialty equipment required!, except a GC-MS to see if this method is practical). However, if you think that there's a niche in chemical engineering where Julia and its associated packages would be more useful, tell me as I'm very flexible. Also, if you think there's any packages that I've missed, tell me. However, I've really enjoyed doing parameters because of how intuitive it is (at least for me) to reason about and design an experiment to get the data you need. If you want to see my intuition at work check "preliminary notes on doing thermal fits.md"

DiffEqFlux.jl 
    - I think this one would be absolutely amazing for reducing the stiffness (and thus performance) of certain elements of the simulations that I run
    - For example, vaporization models or chemical kinetics models are famously stiff while stuff like heat transfer can be simulated normally

Surrogates.jl
    - This one seems riskier and I have a feeling that it would perform significantly worse (at least in terms of training) compared to DiffEqFlux due to the fact that we basically have to account for variations in all of the parameters (or at least I think we do)
    - I sometimes wonder if a fully trained Surrogates.jl model that simulated stuff like my packed_bed_reactor simulation would be faster than just using DiffEqFlux.jl and then solving it with a very fast solver like KrylovJL_GMRES with a iluzero or algebraicmultigrid preconditioner
    - This is the one that I'm probably the least enthusiastic about 

    - My skepticism was warrented, this probably isn't that useful
    - I think this would only be useful for creating digitial twins once the actual reactor is running

General Note:
    - When it comes to training these models I wonder if it would be faster to run it with just FBDF which supports ForwardDiff, FBDF with KrylovJL_GMRES with some kind of reverse differentiation (although reverse differentiators have been inconsistent), or just pure FBDF and KrylovJL_GMRES without any automatic differentiation 
    - Note that KrylovJL_GMRES or basically any GMRES solver offered for Julia does not support forward mode automatic differentiation 
    - I think that FBDF with KrylovJL_GMRES and either Mooncake or Enzyme will be the best
    

SymbolicRegression.jl
    - This one would work great with DiffEqFlux.jl because within a finite volume method simulation that might use DiffEqFlux for fluid dynamics or chemical kinetics, once we have the DiffEqFlux.jl properly calibrated to either simulated or real-world data, we can use that model to do symbolic regression
    - I think with this we could even do something like independently discovering the Ergun Equation or the countless heat exchanger empirical models out there
    - This would be insanely useful as we could probably arrive at something like the Ergun Equation with only a few weeks to a few months of work for cheaply as well
    - Furthermore, we might even be able to develop a DiffEqFlux model and then use symbolic regression to generate an equation that approximates Navier-Stokes CFD for a lower computational cost by taking in data from somehting like OpenFOAM or Ansys Fluent
    - I think being able to pull from industry solvers that have been validated against hundreds of test-cases will be very useful
    - We could maybe even create a modification of the Ergun Equation that accounts for the pressure drop due to rough walls
    - I actually wonder if doing this would be so easy that instead of creating a model that tries to encapsulate all bed void fractions and all media particle sizes, you would just build an apparatus to create a correlation just for that specific geometry of packing media
        - This would be great for gyroids 
    - We could even use this one to generate an equation for an arduino/microcontroller to control a reactor for example based on several outside measurements 

NeuralPDE.jl
    - This one is not useful at all


Other packages to look at:
    - GlobalSensitivty.jl 
        - WE NEED THIS 
    - DataDrivenDiffEq.jl
        - Apparently better than SymbolicRegression for finding physics formulas 
    - Catalyst.jl
        - I should probably look into this one just to make sure I don't need it
        - I wonder how expensive it is compared to the chemical kinetics equations that my simulators uses 
    - JuMP.jl
