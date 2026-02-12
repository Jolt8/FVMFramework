#not sure if this should be called heat or temp cap, 
#it's turning heat input into a temp change, maybe heat to temp then?
function cap_heat_flux_to_temp_change!( 
        du, u, 
        cell_id,
        vol
    )
    # J/s /= m^3 * kg*m^3 * J/(kg*K)
    # = K/s
    du.temp[cell_id] /= vol * u.rho[cell_id] * u.cp[cell_id]
end

#turns rho changes to pressure changes
function cap_mass_flux_to_pressure_change!(
        du, u, 
        cell_id,
        vol
    )
    # kg/s /= (m^3 / (J/(mol*K) * K))
    #remember, J = Pa*m^3
    # = Pa/s
    du.pressure[cell_id] /= (vol / (R_gas * u.temp[cell_id]))
end


