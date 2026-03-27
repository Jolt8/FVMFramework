
function _build_blank_dict!(current_dict, properties, n_cells)
    for (property_name, value) in pairs(properties)
        if value isa NamedTuple
            if !haskey(current_dict, property_name)
                current_dict[property_name] = Dict{Symbol,Any}()
            end
            _build_blank_dict!(current_dict[property_name], value, n_cells)
        else
            if !haskey(current_dict, property_name)
                current_dict[property_name] = zeros(n_cells) .* upreferred(unit(value))
            end
        end
    end
end

function _dict_to_namedtuple(d::Dict)
    # Recursively convert Dicts to NamedTuples
    named_tuple_pairs = Pair{Symbol, Any}[]
    for (property_name, value) in pairs(d)
        if value isa Dict
            push!(named_tuple_pairs, property_name => _dict_to_namedtuple(value))
        else
            push!(named_tuple_pairs, property_name => value)
        end
    end
    return (; named_tuple_pairs...)
end


function _drill_down_and_fill_properties!(merged_properties, properties, cells)
    for (property_name, value) in pairs(properties)
        if merged_properties[property_name] isa NamedTuple
            _drill_down_and_fill_properties!(merged_properties[property_name], properties[property_name], cells)
        elseif merged_properties[property_name] isa AbstractArray
            for cell_id in cells
                merged_properties[property_name][cell_id] = properties[property_name]
            end
        elseif merged_properties[property_name] isa Number
            merged_properties[property_name] = properties[property_name]
        else
            error("The property $property_name was not handled")
        end
    end
end

function _drill_down_and_fill_patch_properties!(merged_properties, properties, cells)
    for (property_name, value) in pairs(properties)
        if merged_properties[property_name] isa NamedTuple
            _drill_down_and_fill_patch_properties!(merged_properties[property_name], properties[property_name], cells)
        elseif merged_properties[property_name] isa AbstractArray
            merged_properties[property_name][1] = properties[property_name]
        elseif merged_properties[property_name] isa Number
            merged_properties[property_name] = properties[property_name]
        else
            error("The property $property_name was not handled")
        end
    end
end

function merge_region_properties(config)
    n_cells = length(config.geo.cell_volumes)

    prop_dict = Dict{Symbol,Any}()

    #convert component arrays to named tuples
    for region in config.regions
        _build_blank_dict!(prop_dict, NamedTuple(region.properties), n_cells)
    end

    for patch in config.patches
        _build_blank_dict!(prop_dict, NamedTuple(patch.properties), n_cells)
    end

    merged_properties = ComponentArray(_dict_to_namedtuple(prop_dict))

    for region in config.regions
        _drill_down_and_fill_properties!(merged_properties, NamedTuple(region.properties), region.region_cells)
    end

    for patch in config.patches
        _drill_down_and_fill_patch_properties!(merged_properties, NamedTuple(patch.properties), patch.cell_neighbors)
    end

    return merged_properties
end


function _drill_down_and_fill_caches!(merged_caches, region_cache_syms_and_units, special_caches, merged_properties)
    for (property_name, property_unit) in pairs(region_cache_syms_and_units)
        if property_name in keys(special_caches)

        elseif merged_caches[property_name] isa NamedTuple
            _drill_down_and_fill_caches!(merged_caches[property_name], region_cache_syms_and_units, (_ = 0.0,), merged_properties)
        elseif merged_caches[property_name] isa AbstractArray
            #since we use the properties as the initial value for the cache, we might not even need properties and could just preemptively merge them
            if hasproperty(merged_properties, property_name) && merged_properties[property_name] isa AbstractArray
                merged_caches[property_name] .= merged_properties[property_name]
            else
                merged_caches[property_name] .= 0.0 * upreferred(property_unit)
            end
        else
            error("The cache $property_name was not handled")
        end
    end
end

function merge_region_caches(config, special_caches, merged_properties)
    n_cells = length(config.geo.cell_volumes)

    cache_dict = Dict{Symbol, Any}()

    for region in config.regions
        for (cache_name, cache_unit) in pairs(region.cache_syms_and_units)
            if !haskey(cache_dict, cache_name)
                if cache_name in keys(special_caches)
                    cache_dict[cache_name] = special_caches[cache_name]
                else
                    cache_dict[cache_name] = zeros(n_cells) .* upreferred(cache_unit)
                end
            end
        end
    end

    merged_caches = _dict_to_namedtuple(cache_dict)

    for region in config.regions
        _drill_down_and_fill_caches!(merged_caches, region.cache_syms_and_units, NamedTuple(special_caches), merged_properties)
    end

    return ComponentArray(merged_caches)
end