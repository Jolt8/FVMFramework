function rebuild_u_named(u_flat, u_proto)
    u_axes = getaxes(u_proto)

    u_named = ComponentVector[]

    for step_id in eachindex(u_flat)
        push!(u_named, ComponentVector(u_flat[step_id], u_axes))
    end

    return u_named
end

function rebuild_u_named_vel(u_flat, u_proto)
    u_axes = getaxes(u_proto)

    pre_u_named = ComponentVector[]

    for step in eachindex(u_flat)
        push!(pre_u_named, ComponentVector(u_flat[step], u_axes))
    end

    u_named = ComponentVector[]

    for step in eachindex(pre_u_named)
        # Create a velocity matrix of size (3, n_cells) for ParaView
        # vcat of adjoints (vx', vy', vz') creates a single 3xN Matrix
        vx = pre_u_named[step].vel_x
        vy = pre_u_named[step].vel_y
        vz = pre_u_named[step].vel_z
        vel_matrix = vcat(vx', vy', vz')

        curr_u_named = ComponentVector(
            velocity = vel_matrix,
            pressure = pre_u_named[step].pressure, 
            temp = pre_u_named[step].temp,
            mass_fractions = pre_u_named[step].mass_fractions
        )

        push!(u_named, curr_u_named)
    end

    return u_named
end