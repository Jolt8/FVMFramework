using Unitful, UnitfulAssets
using Plots

# --- Constants ---
syngas_input = 10.0u"kg/s"  # Let's look at your 10kg/s scale
interest_rate = 0.08

# Energy Content (LHV)
# 1kg of syngas (H2+CO) has ~18 MJ of chemical energy
# 1kg of Methanol has ~20 MJ of chemical energy
syngas_lhv = 18.0u"MJ/kg"
methanol_lhv = 20.0u"MJ/kg"

# --- 1. SOFC Electrical Generation ---
function sofc_cost_per_gj(power_production_mw)
    # SOFCs are expensive: ~$2000/kW for the system
    cost_b = 100u"MUSD" |> u"USD" 
    capacity_b = 50.0u"MW" 
    
    # Scaling factor for fuel cells is often higher (0.8) because they are modular
    cost_of_plant = cost_b * (power_production_mw / capacity_b)^0.8
    
    # CRITICAL: SOFC stacks degrade. 7 year stack life is standard.
    lifetime = 7 
    crf = (interest_rate * (1.0 + interest_rate)^lifetime) / ((1.0 + interest_rate)^lifetime - 1.0)
    
    annual_energy = power_production_mw * 365.25 * 24 * 3600u"s"
    annual_cost = (cost_of_plant * crf) + (cost_of_plant * 0.05) # 5% OpEx
    
    return annual_cost / annual_energy |> u"USD/GJ"
end

# --- 2. Methanol Synthesis ---
function methanol_cost_per_gj(syngas_input)
    # Methanol plants scale better (0.65)
    # $200M for a ~200MW thermal input plant
    cost_b = 2.2u"GUSD" |> u"USD"
    capacity_b = 114u"kg/s"

    syngas_kg_to_kg_methanol = 1.0u"kg/kg"

    methanol_capacity = syngas_input * syngas_kg_to_kg_methanol
    
    cost_of_plant = cost_b * (methanol_capacity / capacity_b)^0.65
    
    lifetime = 25 # Chemical plants last much longer than fuel cells
    crf = (interest_rate * (1.0 + interest_rate)^lifetime) / ((1.0 + interest_rate)^lifetime - 1.0)
    
    annual_mass = methanol_capacity * 365.25 * 24 * 3600u"s"
    annual_cost = (cost_of_plant * crf) + (cost_of_plant * 0.04)
    
    return annual_cost / annual_mass |> u"USD/kg"
end

function fischer_tropsch_economy(syngas_input)
    cost_b = 20u"GUSD" |> u"USD"
    SAF_capacity_b = 220u"kg/s"

    SAF_capacity_a = syngas_input * 0.3 + syngas_input * 0.1
    
    cost_of_plant = cost_b * (SAF_capacity_a / SAF_capacity_b)^0.7
    
    lifetime = 25 # Chemical plants last much longer than fuel cells
    crf = (interest_rate * (1.0 + interest_rate)^lifetime) / ((1.0 + interest_rate)^lifetime - 1.0)
    
    annual_mass = SAF_capacity_a * 365.25 * 24 * 3600u"s"
    annual_cost = (cost_of_plant * crf) + (cost_of_plant * 0.04)
    
    return annual_cost / annual_mass |> u"USD/kg"
end

function haber_bosch_process(syngas_input)
    cost_b = 675u"MUSD" |> u"USD"
    NH3_capacity_b = 23.8u"kg/s"

    NH3_capacity_a = syngas_input * 0.126u"kg/kg" #for a H2:CO ratio of 2:1
    
    cost_of_plant = cost_b * (NH3_capacity_a / NH3_capacity_b)^0.65
    
    lifetime = 25 # Chemical plants last much longer than fuel cells
    crf = (interest_rate * (1.0 + interest_rate)^lifetime) / ((1.0 + interest_rate)^lifetime - 1.0)
    
    annual_mass = NH3_capacity_a * 365.25 * 24 * 3600u"s"
    annual_cost = (cost_of_plant * crf) + (cost_of_plant * 0.04)
    
    return annual_cost / annual_mass |> u"USD/kg"
end


# --- Performance Calcs ---
sofc_efficiency = 0.60

power_out_sofc = (syngas_input * syngas_lhv * sofc_efficiency) |> u"MW"

# Revenue
elec_price = 39.0u"USD/GJ"       # ~$0.14/kWh (High/Green price)
meth_price = 500u"USD/1000kg" |> u"USD/kg"
SAF_price = 2.00u"USD/kg"

naphtha_yield_per_kg_syngas = 0.10u"kg/kg" 
naphtha_price = 0.70u"USD/kg"

NH3_price = 700u"USD/1000kg" |> u"USD/kg"

# Profitability

function calc_sofc_profit(syngas_input)
    power_out_sofc = syngas_input * syngas_lhv * sofc_efficiency
    sofc_profit = (elec_price - sofc_cost_per_gj(power_out_sofc)) * power_out_sofc
    return sofc_profit |> u"USD/s"
end

function calc_methanol_profit(syngas_input)
    meth_profit = (meth_price - methanol_cost_per_gj(syngas_input)) * syngas_input
    return meth_profit |> u"USD/s"
end

function calc_fischer_tropsch_profit(syngas_input)
    ft_unit_cost_saf = fischer_tropsch_economy(syngas_input)
    ft_saf_revenue = (SAF_price - ft_unit_cost_saf) * (syngas_input * 0.3)
    ft_naphtha_revenue = naphtha_price * (syngas_input * naphtha_yield_per_kg_syngas)
    return (ft_saf_revenue + ft_naphtha_revenue) |> u"USD/s"
end

function calc_haber_bosch_profit(syngas_input)
    hb_unit_cost_nh3 = haber_bosch_process(syngas_input)
    hb_nh3_revenue = (NH3_price - hb_unit_cost_nh3) * (syngas_input * 0.126u"kg/kg")
    return hb_nh3_revenue |> u"USD/s"
end

println("SOFC Profit: ", calc_sofc_profit(syngas_input) |> u"USD/hr")
println("Methanol Profit: ", calc_methanol_profit(syngas_input) |> u"USD/hr")
println("Fischer Tropsch Profit: ", calc_fischer_tropsch_profit(syngas_input) |> u"USD/hr")
println("Haber Bosch Profit: ", calc_haber_bosch_profit(syngas_input) |> u"USD/hr")

#hmm, it appears Methanol is the best option up until around 18 kg/s
#SOFC is never worth it
#how is the Haber-Bosch process ever profitable

syngas_input_values = collect(range(0.1u"kg/s", 10u"kg/s", 100))
plot(syngas_input_values, calc_sofc_profit.(syngas_input_values), label = "SOFC", xlabel = "Syngas Input (kg/s)", ylabel = "Profit (USD/s)", legend = :topleft)
plot!(syngas_input_values, calc_methanol_profit.(syngas_input_values), label = "Methanol")
plot!(syngas_input_values, calc_fischer_tropsch_profit.(syngas_input_values), label = "Fischer Tropsch")
plot!(syngas_input_values, calc_haber_bosch_profit.(syngas_input_values), label = "Haber Bosch")
