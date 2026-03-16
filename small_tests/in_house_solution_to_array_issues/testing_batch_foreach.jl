include("FVMArray.jl")

u_caches = (
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), #it would also be nice if commans were not required for fields with 1 field within them
    rho = zeros(n_cells)
)

u_caches_axes, u_caches_offset = create_axes(u_caches)
u_caches_vec = zeros(u_caches_offset)

u_cache_test = FVMArray(u_caches_vec, u_caches_axes)

function test_foreach_field(u_cache_test)
    @batch for i in 1:100
        foreach_field_at!(1, u_cache_test.net_rates.reforming_reactions) do reaction, net_rate
            net_rate[reaction] += 1.0
        end
    end
end

u_cache_test .= 0.0
@btime test_foreach_field(u_cache_test) #275.415 ns, (0 allocations, 0 bytes)

function test_batch_foreach_field(u_cache_test)
    @batch for i in 1:100
        batch_foreach_field(u_cache_test.net_rates.reforming_reactions) do reaction, net_rate #this allocates a lot
            net_rate .+= 1.0
            #u_cache_test.net_rates.reforming_reactions[reaction] .+= 1.0
        end
    end
end

u_cache_test .= 0.0
@btime test_batch_foreach_field(u_cache_test) #83.1 μs (1902 allocations, 68.84 KiB)
