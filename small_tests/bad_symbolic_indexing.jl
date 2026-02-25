using ComponentArrays
using Polyester
using BenchmarkTools

struct RegionU{S1, S2}
    state::ComponentVector
    cache::ComponentVector
    fixed::NamedTuple
end

struct RegionTest{S1, S2}
    properties::NamedTuple
end

function FVMSystemAccessor(u::ComponentVector, cache::ComponentVector, region::RegionTest{S1, S2}) where {S1, S2}
    return RegionU{S1, S2}(
        u,
        cache,
        region.properties
    )
end

@generated function Base.getproperty(u::RegionU{S1, S2}, sym::Symbol) where {S1, S2}
    return quote
        if sym in $S1
            return getproperty(getfield(u, :state), sym)
        elseif sym in $S2
            return getproperty(getfield(u, :cache), sym)
        elseif sym === :state || sym === :cache || sym === :fixed
            return getfield(u, sym)
        else
            return getproperty(getfield(u, :fixed), sym)
        end
    end
end

@generated function Base.setproperty!(u::RegionU{S1, S2}, sym::Symbol, val) where {S1, S2}
    return quote
        if sym in $S1
            setproperty!(getfield(u, :state), sym, val)
        elseif sym in $S2
            setproperty!(getfield(u, :cache), sym, val)
        else
            setproperty!(getfield(u, :fixed), sym, val)
        end
    end
end

n_cells = 10000

du_proto = ComponentArray(
    mass_fractions=ComponentVector(methanol=zeros(n_cells), water=zeros(n_cells), carbon_monoxide=zeros(n_cells), hydrogen=zeros(n_cells), carbon_dioxide=zeros(n_cells)),
    temp=zeros(n_cells)
)

du_proto_nt = (; du_proto...)

u_proto = ComponentArray(
    mass_fractions=ComponentVector(methanol=zeros(n_cells), water=zeros(n_cells), carbon_monoxide=zeros(n_cells), hydrogen=zeros(n_cells), carbon_dioxide=zeros(n_cells)),
    temp=zeros(n_cells)
)

u_proto_nt = (; u_proto...)

du_cache = ComponentArray(
    rho=zeros(n_cells),)

du_cache_nt = (; du_cache...)

u_cache = ComponentArray(
    rho=zeros(n_cells)
)

u_cache_nt = (; u_cache...)

properties = (cp=zeros(n_cells), k=zeros(n_cells))

reforming_area_region = RegionTest{(:mass_fractions, :temp),(:rho,)}(
    properties
)

du = ComponentVector(du_proto; NamedTuple(du_cache)..., properties...)
#FVMSystemAccessor(du_proto, du_cache, reforming_area_region)
u = ComponentVector(u_proto; NamedTuple(u_cache)..., properties...)

du_axes = getaxes(du)[1]
u_axes = getaxes(u)[1]

du = Vector(du)
u = Vector(u)

#ahhh, so turning them all into a component vector is what's causing multithreading to fail

function test_accessors(du_proto_nt, du_cache_nt, u_proto_nt, u_cache_nt, properties, n_cells)
    du = (; du_proto_nt..., du_cache_nt...)
    u = (; u_proto_nt..., u_cache_nt..., properties...)
    for cell_id in 1:n_cells
        du.mass_fractions
        du.mass_fractions.methanol
        du.mass_fractions.methanol[cell_id]
        du.mass_fractions.methanol[cell_id] = 1.0
        u.mass_fractions
        u.mass_fractions.methanol
        u.mass_fractions.methanol[cell_id]
        u.mass_fractions.methanol[cell_id] = 1.0

        du.temp
        du.temp[cell_id]
        du.temp[cell_id] = 270.0
        u.temp
        u.temp[cell_id]
        u.temp[cell_id] = 270.0

        du.rho
        du.rho[cell_id]
        du.rho[cell_id] = 1.0
        u.rho
        u.rho[cell_id]
        u.rho[cell_id] = 1.0

        u.cp
        u.cp[cell_id]
        u.cp[cell_id] = 1.0
    end
end

function test_multiple_samples(du_proto, du_cache, u_proto, u_cache, properties, n_cells, n_samples)
    for _ in 1:n_samples
        test_accessors(du_proto, du_cache, u_proto, u_cache, properties, n_cells)
    end
end

VSCodeServer.@profview test_multiple_samples(du_proto_nt, du_cache_nt, u_proto_nt, u_cache_nt, properties, n_cells, 1000)
@btime test_accessors($du_proto_nt, $du_cache_nt, $u_proto_nt, $u_cache_nt, $properties, $n_cells)