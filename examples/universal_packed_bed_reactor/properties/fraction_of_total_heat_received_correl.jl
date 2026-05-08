function calculate_fraction_of_total_heat_received!(du, u, geo)

    #this would probably be some kind of polynomial function that takes in the cell id and returns the fraction of total heat received
    #what we'll probably do is have a bunch of different versions of this function and see which one fits the data the best and is also the simplest
    #we don't want the one that fits the data the best 

    #actually instead of using cell_id, we'll probably use the function that calculates the fraction of the total length of the pipe that the cell is at
    #because it makes using various formulas easier

    expected_fraction = 1.0 / length(geo.cell_volumes) #placeholder for now
    middle_cell_length_along_reactor = (u.cell_lengths_along_pipe[1] + u.cell_lengths_along_pipe[end]) / 2
    #cell_length_fraction_from_center = 2 * (cell_id - middle_cell) / (length(cell_volumes) - 1)

    #for 100: 2 * ((100 - 50.5) / (100 - 1)) = 1
    #for 50: 2 * ((50 - 50.5) / (100 - 1)) = -0.010101...
    #for 1: 2 * ((1 - 50.5) / (100 - 1)) = -1
    #perfect!

    for cell_id in eachindex(geo.cell_volumes)
        cell_length_fraction = 2 * (u.cell_lengths_along_pipe[cell_id] - middle_cell_length_along_reactor) / (u.cell_lengths_along_pipe[end] - u.cell_lengths_along_pipe[1])
        u.per_cell_fraction_of_total_heat_received[cell_id] = u.heater_poly_p2[1] * cell_length_fraction^2 + u.heater_poly_p1[1] * cell_length_fraction + u.heater_poly_p0[1]
    end 
end