Notes on things we're probably going to have to create a correlation for:
    1. The fraction of the heat expelled by the electric heater that each cell received represented as a polynomial
    2. The overall heat transfer coefficient between the pipe and the surroundings using the standard formula 
    3. The adjusted heat capacity of the entire reactor (including the pipe itself, ceramic insulation, packing, etc.) which is probably doable because we know the specific heat capacity and density of air and this would probably become a lumped parameter model
    4. Thermal gradients being smoothed out due to the packing of the reactor and the reactor walls
    5. The lag each thermocouple experiences in response to a change in temperature
    6. The difference in temperature between the center of the reactor and the 
    7. The heating each thermocouple experiences due to variable levels of insulation between the thermocouple and the heating wire surrounding it

All of these parameters are going to be pretty hard to decouple and probably harder to decouple than the reaction kinetics themselves

I wonder if we could use something with similar thermal properties to methanol like propylene glycol (or maybe just use pure water) to make constructing a model easier by allowing for more parameters to be lumped together

The thing I'm most concerned about is the fact that the signal that each theromcouple received is not only skewed by the varing levels of insulation between the thermocouple and the heating wire surrounding it, but also by the varying levels of insulation between the heating wire and the pipe itself, and in addition the varying levels of insulation between the thermocouples and the pipe

I think we can decouple the varying levels of insulation between the thermocouple and the pipe by turning the heater off and then flow water through the reactor at a known temperature that's either higher or lower than room temperature for several minutes to allow it to reach steady state 

If we know the time between when the water is introduced and when it comes out the top, we can roughly estimate how long it takes for the water 
to reach the top

If we also measure the inlet and outlet temperature of the water over time, we can roughly estimate the temperature gradient inside the reactor
then, we can use this data to hopefully account for the varying levels of insulation between the thermocouple and the pipe

Also, if we flowed in the water faster or slower we could then decouple the heat transfer along the length of the reactor 

I don't think we could do this by collecting data when the reactor is half full for example due to the fact that liquid filling the void in the bed has a synergistic with increasing heat transfer in the catalyst packing

Wow, I really cannot come up with any other tests to decouple the other parameters, I think we're going to have to create some *very* empirical correlations for them

Yes it can tell the difference. If I have water slowly flow up the walls of the reactor, the difference in temperature between the thermocouple and the water inside the pipe (using the assumed thermal gradient based on the difference between the inlet and outlet temperatures) can find the heat transfer coefficient between the thermocouple and the pipe. Then, by allowing the hot water to flow through the reactor for a long time, we can find the bias of each of the thermocouples.
