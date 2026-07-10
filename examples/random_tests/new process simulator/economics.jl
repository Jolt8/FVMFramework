function component_cost_with_interest(
    streams, stream_id, system,
    desired_capacity, reference_capacity, reference_cost, 
    scaling_exponent, component_lifetime, interest_rate,
    component_operating_costs, component_maintenence_costs
)
    cost_of_desired_capacity_plant = reference_cost * (desired_capacity / reference_capacity)^scaling_exponent

    crf = (interest_rate * (1.0 + interest_rate)^component_lifetime) / ((1.0 + interest_rate)^component_lifetime - 1.0)

    annual_capacity = ethylene_capacity * 365.0 * 24.0 * 60.0 * 60.0u"s"
    annualized_capex = cost_of_desired_capacity_plant * crf

    annual_opex = component_operating_costs + component_maintenence_costs

    #cost_per_capacity = (annualized_capex + annual_opex) / annual_capacity |> Unitful.Unit(1.0u"USD" / (desired_capacity * 0.0))
    #cost_per_total_output = cost_per_capacity * component_lifetime |> Unitful.Unit(1.0u"USD" / (desired_capacity * 0.0u"s"))
    #I think the cost per capacity and cost per total output should be evaluated later once the entire plant has been laid out

    system.total_plant_cost += (annualized_capex + annualized_opex)
end

function component_cost_no_interest(
    streams, stream_id, system,
    desired_capacity, reference_capacity, reference_cost, 
    scaling_exponent, component_lifetime, interest_rate,
    component_operating_costs, component_maintenence_costs
)
    cost_of_desired_capacity_plant = reference_cost * (desired_capacity / reference_capacity)^scaling_exponent

    annual_capacity = ethylene_capacity * 365.0 * 24.0 * 60.0 * 60.0u"s"
    annualized_capex = cost_of_desired_capacity_plant

    annual_opex = component_operating_costs + component_maintenence_costs

    #cost_per_capacity = (annualized_capex + annual_opex) / annual_capacity |> Unitful.Unit(1.0u"USD" / (desired_capacity * 0.0))
    #cost_per_total_output = cost_per_capacity * component_lifetime |> Unitful.Unit(1.0u"USD" / (desired_capacity * 0.0u"s"))
    #I think the cost per capacity and cost per total output should be evaluated later once the entire plant has been laid out

    system.total_plant_cost += (annualized_capex + annualized_opex)
end