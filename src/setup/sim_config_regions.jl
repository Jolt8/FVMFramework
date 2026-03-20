function add_region!(
    config, name;
    type,
    initial_conditions,
    properties,
    optimized_syms,
    cache_syms_and_units,
    region_function
)

    region_cells = collect(getcellset(config.grid, name))

    for cell_id in region_cells
        for property_name in propertynames(initial_conditions)
            var = getproperty(config.u_proto, property_name)
            initial_condition = getproperty(initial_conditions, property_name)
            if var isa ComponentVector
                for sub_name in propertynames(var)
                    getproperty(var, sub_name)[cell_id] = getproperty(initial_condition, sub_name)
                end
            else
                var[cell_id] = initial_condition
            end
        end
    end

    for field in optimized_syms
        config.optimized_parameters[field] = properties[field]
    end

    region = RegionSetupInfo(name, type, initial_conditions, properties, cache_syms_and_units, region_function, region_cells)

    all_region_names = [region.name for region in config.regions]

    if name in all_region_names
        existing_region_idx = findfirst(x -> x == name, all_region_names)

        config.regions[existing_region_idx] = region
    else   
        push!(config.regions, region)

        config.regions[1].cache_syms_and_units = merge(config.regions[1].cache_syms_and_units, optimized_syms) 
    end
    return 
end