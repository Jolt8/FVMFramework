Solver Organization
    - Physics types and properties 
    - Sorting of different physics types into different types of connections
        - (ex. Navier-Stokes on fluid_fluid interfaces)
        - (ex. Conjugate heat transfer on fluid_solid interfaces)
    - Region functions containing the following in the following order of operations:
        - Internal Physics (ex. chemical reactions)
        - Sources (ex. volumetric heating)
        - Boundary conditions 
        - Capacities (ex. dividing heat flux (J/s) by (rho * vol * cp) (kg*m^3 * m^3 * J/(kg*K)) to get K/s)