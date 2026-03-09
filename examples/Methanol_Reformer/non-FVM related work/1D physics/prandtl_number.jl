function prandtl_number!(du, u, cell_id, vol)
    dynamic_viscosity = u.mu[cell_id] #u.viscosity(u.temp, u.pressure)[1] #not sure why this returns a matrix 
    k = u.k[cell_id] #u.k(u.temp, u.pressure)[1]
    cp = u.cp[cell_id]  #u.cp(u.temp, u.pressure)[1]

    u.Pr[cell_id] = (dynamic_viscosity * cp) / (k)
end

#u = (temp = [270.0], pressure = [100000.0], Pr = [0.0], viscosity = dynamic_viscosity_interp_TP, k = k_interp_TP, cp = cp_interp_TP)
#du = (temp = [0.0], pressure = [0.0], Pr = [0.0])

#prandtl_number!(du, u, 1, 1)
#@benchmark prandtl_number!($du, $u, 1, 1)