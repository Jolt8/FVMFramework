function mw_avg!(stream)
    inv_mw = 0.0

    for_fields!(stream.mass_fractions, stream.molecular_weights) do species, mass_fractions, molecular_weights
        inv_mw += mass_fractions[species] / molecular_weights[species]
    end
    
    stream.mw_avg = 1.0 / inv_mw
end

function rho_ideal!(stream)
    stream.rho = (stream.pressure * stream.mw_avg) / (stream.R_gas * stream.temp)
end

function cp_avg!(stream)
    stream.cp_avg *= 0.0

    for_fields!(stream.mass_fractions, stream.species_cps) do species, mass_fractions, species_cps
        stream.cp_avg += mass_fractions[species] * species_cps[species]
    end
end

function cv_avg!(stream)
    stream.cv_avg *= 0.0

    for_fields!(stream.mass_fractions, stream.species_cvs) do species, mass_fractions, species_cvs
        stream.cv_avg += mass_fractions[species] * species_cvs[species]
    end
end

function specific_heat_ratio!(stream)
    stream.specific_heat_ratio = stream.cp_avg / stream.cv_avg
end

function mass_flows!(stream)
    for_fields!(stream.mass_flows, stream.mass_fractions) do species, stream_mass_flows, stream_mass_fractions
        stream_mass_flows[species] = stream.mass_flow * stream_mass_fractions[species]
    end
end

function molar_flow!(stream)
    stream.molar_flow = stream.mass_flow / stream.mw_avg
end

function molar_flows!(stream)
    for_fields!(stream.molar_flows, stream.molar_fractions) do species, stream_molar_flows, stream_molar_fractions
        stream_molar_flows[species] = stream.molar_flow * stream_molar_fractions[species]
    end
end

function volumetric_flow!(stream)
    stream.volumetric_flow = stream.mass_flow / stream.rho
end

function volumetric_flows!(stream)
    for_fields!(stream.volumetric_flows, stream.mass_fractions) do species, stream_volumetric_flows, stream_mass_fractions
        stream_volumetric_flows[species] = (stream.mass_flow * stream_mass_fractions[species]) / stream.rho
    end
end

function molar_concentrations!(stream)
    for_fields!(stream.mass_fractions, stream.molecular_weights, stream.molar_concentrations) do species, mass_fractions, molecular_weights, molar_concentrations
        molar_concentrations[species] = (stream.rho * mass_fractions[species]) / molecular_weights[species]
    end
end

function molar_fractions!(stream)
    for_fields!(stream.mass_fractions, stream.molecular_weights, stream.molar_fractions) do species, mass_fractions, molecular_weights, molar_fractions
        molar_fractions[species] = (mass_fractions[species] / molecular_weights[species]) * stream.mw_avg
    end
end

function update_stream_properties!(stream)
    mw_avg!(stream)
    rho_ideal!(stream)
    molar_concentrations!(stream)
    molar_fractions!(stream)
    mass_flows!(stream)
    molar_flow!(stream)
    molar_flows!(stream)
    volumetric_flow!(stream)
    volumetric_flows!(stream)
    cp_avg!(stream)
    cv_avg!(stream)
    specific_heat_ratio!(stream)
end

