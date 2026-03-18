using ComponentArrays
using Polyester
using BenchmarkTools

# 1. ComponentArray Compile-Time Unroller
# Extracts the top-level property keys inside a ComponentVector at compile time
# and emit fully unrolled code that evaluates the closure `f` for each property.
@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{ComponentVector, N}) where {N}
    # For a ComponentVector, the axes are in the 4th type parameter
    AxTuple = groups[1].parameters[4]
    
    # Extract the first Axis
    Ax = AxTuple.parameters[1]
    
    # The first parameter of Axis is the NamedTuple containing the layout
    Props = keys(Ax.parameters[1])
    
    exprs = []
    for n in Props
        # Emit an expression to get the sub-view for each block
        view_exprs = [:(getproperty(groups[$i], $(QuoteNode(n)))) for i in 1:N]
        # Emit the call to the closure: f(cell_id, :field_name, view1, view2...)
        push!(exprs, :(f(cell_id, $(QuoteNode(n)), $(view_exprs...))))
    end
    
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

# 2. Virtual Component Array
# A lightweight struct holding multiple ComponentArrays and a compile-time map.
# E.g. (mass_fractions = 1, pressure = 1, rho = 2, net_rates = 2, viscosity = 3)
struct VirtualComponentArray{Data <: Tuple, Mapping <: NamedTuple}
    data::Data
    mapping::Mapping
end

# Transparently route nested field access to the appropriate underlying buffer
@inline function Base.getproperty(V::VirtualComponentArray, s::Symbol)
    # At compile time, standard Julia infers which index to use from mapping
    map_idx = getfield(getfield(V, :mapping), s)
    # Then accesses the field from the specific ComponentArray natively
    return getproperty(getfield(V, :data)[map_idx], s)
end

@inline Base.getindex(V::VirtualComponentArray, s::Symbol) = getproperty(V, s)
@inline Base.keys(V::VirtualComponentArray) = keys(getfield(V, :mapping))


# ==============================================================================
# Testing & Benchmarking 
# ==============================================================================

n_cells = 100000

# Group 1: State vars (u, du)
u_proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)
u_arr = ComponentArray(u_proto)
du_arr = ComponentArray(u_proto)

# Group 2: Intermediates using DiffCache (Caches)
cache_proto = (
    mass = zeros(n_cells),
    rho = ones(n_cells),
    net_rates = (
        reforming_reactions = (WGS_rxn = zeros(n_cells), MD_rxn = zeros(n_cells)),
    )
)
cache_tmp = ComponentArray(cache_proto)

# Group 3: Global or constant properties
prop_proto = (
    viscosity = ones(n_cells),
)
prop_arr = ComponentArray(prop_proto)

# The routing map. We route dot keys to the corresponding array index
# 1 -> State variables (u/du)
# 2 -> Intermediate calculations (caches)
# 3 -> Properties
vmap = (
    mass_fractions = 1,
    pressure = 1,
    mass = 2,
    rho = 2,
    net_rates = 2,
    viscosity = 3
)

# Build the abstracted Virtual Components
u = VirtualComponentArray((u_arr, cache_tmp, prop_arr), vmap)
du = VirtualComponentArray((du_arr, cache_tmp, prop_arr), vmap)

# Core physics solver
function test_virtual_component_loop!(du, u, n_cells)
    @batch for cell_id in 1:n_cells
        # Direct dot syntax routes without allocation
        du.mass_fractions.methylene_blue[cell_id] += 1.0 
        
        # Unrolled iteration over mass_fractions (Routes to tuple index 1)
        foreach_field_at!(cell_id, du.mass_fractions) do cid, species_name, mass_fraction
            mass_fraction[cid] += 1.0 
        end
        
        # Unrolled iteration over sub-sub caches (Routes to tuple index 2)
        foreach_field_at!(cell_id, du.net_rates.reforming_reactions) do cid, rxn_name, rxn_rate
            rxn_rate[cid] += 1.0 
        end
        
        # Mixes everything (Routes to 1, 2, and 3 internally)
        du.mass[cell_id] += du.rho[cell_id] * u.viscosity[cell_id]
    end
end

println("="^50)
println("Running preliminary unallocated check...")
test_virtual_component_loop!(du, u, n_cells)
println("Did u.mass_fractions.methylene_blue update? ", du_arr.mass_fractions.methylene_blue[1] == 2.0)

println("\nChecking @btime performance for virtualizing...")
VSCodeServer.@profview test_virtual_component_loop!(du, u, n_cells)
println("="^50)
