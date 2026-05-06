using Unitful, UnitfulAssets
using MonteCarloMeasurements
using Plots

#works:
#test = Particles(collect(1u"kg":1u"kg":10u"kg")) #this does give a display error, however
Particles(2000, Normal(0, 1)) * u"s"

#doesn't work:
#1.0u"kg" ± 1.0u"kg"

#anyways, back to the script
H2_kg_per_kg_plastic = 0.3u"kg/kg" #0.2-0.4
H2_mol_per_kg_plastic = H2_kg_per_kg_plastic / 2.016u"g/mol"

CO_kg_per_kg_plastic = 0.75 #0.5-1.0
CO_mol_per_kg_plastic = CO_kg_per_kg_plastic / 28.01u"g/mol"

test_1 = ustrip(H2_mol_per_kg_plastic / H2_mol_per_kg_plastic)
test_2 = ustrip(CO_mol_per_kg_plastic / H2_mol_per_kg_plastic)

H2_to_CO_ratio = ustrip(H2_mol_per_kg_plastic) / ustrip(CO_mol_per_kg_plastic)

#eh, we'll just say that a 2.5 ratio can be achieved for ZEO OX

plastic_cost_per_kg = 0.07u"USD/kg" #0.02 - 0.15 #INPUT

kg_syngas_per_kg_plastic = 2.0u"kg/kg" #1.5-2.5 

kg_ethylene_prod_per_kg_syngas = 0.31u"kg/kg" #0.333-0.286

ethylene_cost_per_kg = 0.9u"USD/kg" #1.21 - 0.53 #OUPUT

ethylene_prod_per_kg_plastic = kg_syngas_per_kg_plastic * kg_ethylene_prod_per_kg_syngas 

#Using Methanol-to-Olefins (MTO) as a closer analog to OX-ZEO than Fischer-Tropsch
function cost_per_kg_of_ethylene_based_on_kgs_of_capacity_of_plant(ethylene_capacity)
    #cost_b = 1.5e9u"USD"
    #capacity_b = 19.0u"kg/s"
    cost_b = 190u"MUSD"
    capacity_b = 4.75u"kg/s"

    cost_of_desired_capacity_plant = cost_b * (ethylene_capacity / capacity_b)^0.65

    plant_lifetime_years = 25 

    interest_rate = 0.08
    crf = (interest_rate * (1.0 + interest_rate)^plant_lifetime_years) / ((1.0 + interest_rate)^plant_lifetime_years - 1.0)

    annual_capacity = ethylene_capacity * 365.0 * 24.0 * 60.0 * 60.0u"s"
    annualized_capex = cost_of_desired_capacity_plant * crf

    annual_opex = cost_of_desired_capacity_plant * 0.04 # 4% for maintenance and labor

    return (annualized_capex + annual_opex) / annual_capacity
end

capacity_range = collect(range(0.46293629362936295u"kg/s", 25.0u"kg/s", 100))
ox_zeo_plant_added_cost_per_kg_ethylene = cost_per_kg_of_ethylene_based_on_kgs_of_capacity_of_plant.(capacity_range)
plot(capacity_range, cost_per_kg_of_ethylene_based_on_kgs_of_capacity_of_plant.(capacity_range))

ethylene_cost_per_kg_vec = [ethylene_cost_per_kg - cost_per_kg_of_ethylene_based_on_kgs_of_capacity_of_plant(cap) for cap in capacity_range]
#min_profitable_idx = argmin(abs.(ethylene_cost_per_kg_vec .- 0.0u"USD/kg"))
#capacity_range[min_profitable_idx]

plot(capacity_range, ethylene_cost_per_kg_vec)
#10kg/s seems to be optimal

ethylene_sell_per_kg = ethylene_cost_per_kg - cost_per_kg_of_ethylene_based_on_kgs_of_capacity_of_plant(10.0u"kg/s")

plastic_cost_per_kg
ethylene_sell_per_kg
ethylene_prod_per_kg_plastic

profit_per_kg_plastic = ethylene_prod_per_kg_plastic * ethylene_sell_per_kg - plastic_cost_per_kg
#oof, 0.087 USD/kg 

function cost_per_kg_of_plastic_processed_based_on_kgs_of_capacity_of_plant(plastic_capacity)
    cost_b = 500u"MUSD" |> u"USD"
    capacity_b = 94.49841u"kg/s"

    cost_of_desired_capacity_plant = cost_b * (plastic_capacity / capacity_b)^0.65 #also 0.65

    plant_lifetime_years = 25 #also 25 years
    
    interest_rate = 0.08
    crf = (interest_rate * (1.0 + interest_rate)^plant_lifetime_years) / ((1.0 + interest_rate)^plant_lifetime_years - 1.0)

    annual_capacity = plastic_capacity * 365.0 * 24.0 * 60.0 * 60.0u"s"
    annualized_capex = cost_of_desired_capacity_plant * crf

    annual_opex = cost_of_desired_capacity_plant * 0.04 # 4% for maintenance and labor

    plant_cost_per_kg_plastic_processed = (annualized_capex + annual_opex) / annual_capacity

    #energy usage
    solar_cost_per_j_over_lifetime = 1.25e-8u"USD/J" #1.1-1.4
    #not including batteries 

    #OX-ZEO is exothermic, but to assume a worse case scenario, we won't use its heat of reaction to preheat the incoming streams

    plastic_energy_usage_per_kg = 3.0u"MJ/kg" |> u"J/kg" #this is on the high end, usually 1.3-3.6 MJ/kg

    cost_of_energy_per_kg_plastic = solar_cost_per_j_over_lifetime * plastic_energy_usage_per_kg
    #only accounts for 0.0375 USD/kg

    steam_cost_per_kg_plastic = 0.02u"USD/kg" # Cost of steam for gasification

    return plant_cost_per_kg_plastic_processed + cost_of_energy_per_kg_plastic + steam_cost_per_kg_plastic
end

plastic_capacity_range = collect(range(0.46293629362936295u"kg/s", 25.0u"kg/s", 100))
corresponding_ethylene_capacity_range = plastic_capacity_range * ethylene_prod_per_kg_plastic

function profit_per_kg_plastic_1(plastic_capacity)
    plastic_processing_plant_cost_per_kg = cost_per_kg_of_plastic_processed_based_on_kgs_of_capacity_of_plant(plastic_capacity)
    corresponding_ethylene_capacity = ethylene_prod_per_kg_plastic * plastic_capacity
    ethylene_plant_cost_per_kg = cost_per_kg_of_ethylene_based_on_kgs_of_capacity_of_plant(corresponding_ethylene_capacity)

    profit_per_kg_plastic = ethylene_prod_per_kg_plastic * (ethylene_cost_per_kg - ethylene_plant_cost_per_kg) - plastic_cost_per_kg - plastic_processing_plant_cost_per_kg
    return profit_per_kg_plastic
end

profit_per_kg_plastic_range = profit_per_kg_plastic_1.(plastic_capacity_range)
plot(plastic_capacity_range, profit_per_kg_plastic_range)
#around a 0.2$ profit at 10kg/s of plastic processed
#around a 0.3$ profit at 25kg/s of plastic processed
#pretty good!
#I'm sure that other costs would eat into that heavily, but economy of scale estimate for the ox-zeo plant 
#was probably way too high since I used the fischer tropsch process as a reference and the 
#ethylene plant is taking up most of the cost per kg

min_profitable_plastic_capacity = argmin(abs.(profit_per_kg_plastic_range .- 0.0u"USD/kg"))
plastic_capacity_range[min_profitable_plastic_capacity] #1.7 kg/s