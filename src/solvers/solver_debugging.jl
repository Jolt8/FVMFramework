function check_state(state, region, cell_id, suspicious_values)
    for field in keys(state)
        state_field = state[field]
        if state_field isa NamedTuple
            for (property_name, property_field) in pairs(state_field)
                if any(in(suspicious_values), property_field) || any(isnan, property_field)
                    println("$(field).$(property_name) is $(property_field) in this region $(region.name)")
                end
            end
        elseif state_field isa AbstractArray && !(state_field isa SubArray)
            if state_field[cell_id] in suspicious_values || isnan(state_field[cell_id])
                println("$(field) is $(state_field[cell_id]) in this region $(region.name)")
            end
        else
            if state_field[1] in suspicious_values || isnan(state_field[1])
                println("$(field) is $(state_field[1]) in this region $(region.name)")
            end
        end
    end
end

function debug_region!(du, u, region)
    cell_id = region.region_cells[1]
    suspicious_du_values = [-Inf, Inf, NaN]
    check_state(du, region, cell_id, suspicious_du_values)

    suspicious_u_values = [-Inf, Inf, NaN]
    check_state(u, region, cell_id, suspicious_u_values)
end