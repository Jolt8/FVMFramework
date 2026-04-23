using Unitful
using Plots

#this is rough calculation of the energy required for direct reduction of iron ore using 
#methane reforming inside the furnace

heat_energy_per_kg_of_iron = 10u"MJ/kg"

typical_furnace_processed_kg_per_day = 10000u"1000kg/24hr"

typical_furnace_heat_rate = heat_energy_per_kg_of_iron * typical_furnace_processed_kg_per_day |> u"MW" #1157 MW !!!

concentrated_solar_power_wattage_per_area = 100u"MW/3.03514km^2"

typical_furnace_area = typical_furnace_heat_rate / concentrated_solar_power_wattage_per_area |> u"km^2" #35 km^2




#comparing solar to methane combustion

#methane combustion
wattage_required = typical_furnace_heat_rate |> u"GW"

CO2_per_GJ = 53u"kg/GJ"
carbon_tax = 100u"1/1000kg"
dollars_per_GJ = ustrip(CO2_per_GJ * carbon_tax) * 1/u"GJ" #5.3 dollars per GJ
methane_cost_per_energy = 2.5u"1/GJ"

methane_cost_per_second = (wattage_required * (methane_cost_per_energy + dollars_per_GJ)) #2.31 CAD/s

methane_cost_per_day = ustrip(methane_cost_per_second * (60 * 60 * 24)) #200,000 CAD/day

#solar
typical_solar_cost_per_GW = 2_000_000_000u"1/GW" #cad

total_solar_cost = ustrip(wattage_required * typical_solar_cost_per_GW) #1.023 * 890,000,000 = 910,470,000

hours_to_payback_solar = total_solar_cost / methane_cost_per_day * u"24hr"
year_to_payback_solar = hours_to_payback_solar |> u"yr"
#25.3 years with 0$ carbon tax

#8.125 years with 100$ carbon tax

#solar is not worth it :(
#even with concentrated solar power which would probably be cheaper than solar panels, 
#the fact that it's intermittend makes it not worth it when constant upkeep is basically required for economies of scale

#however, this is not taking into account any carbon taxes or the cost of doing something with the CO2 emissions
#which would probably make solar much more attractive
#yeah, if a 100$ / ton carbon tax was used, natural gas would be around 9 dollars per GJ, which would make solar more than worth it
#I can see why a carbon tax is commonly cited as the best way to fight climate change


function carbon_tax_vs_solar_payback(carbon_tax)
    CO2_per_GJ = 53u"kg/GJ"
    dollars_per_GJ = ustrip(CO2_per_GJ * carbon_tax) * 1/u"GJ" #5.3 dollars per GJ
    methane_cost_per_energy = 2.5u"1/GJ"

    methane_cost_per_second = (wattage_required * (methane_cost_per_energy + dollars_per_GJ)) #2.31 CAD/s

    methane_cost_per_day = ustrip(methane_cost_per_second * (60 * 60 * 24)) #200,000 CAD/day

    #solar
    typical_solar_cost_per_GW = 3_500_000_000u"1/GW" #cad

    total_solar_cost = ustrip(wattage_required * typical_solar_cost_per_GW) #1.023 * 890,000,000 = 910,470,000

    hours_to_payback_solar = total_solar_cost / methane_cost_per_day * u"24hr"
    years_to_payback_solar = hours_to_payback_solar |> u"yr"

    return years_to_payback_solar
end

carbon_tax_vs_solar_payback(37u"1/1000kg")

carbon_taxes = collect(0:10:300) .* u"1/1000kg"

payback_times = [carbon_tax_vs_solar_payback(carbon_tax) for carbon_tax in carbon_taxes]

plot(carbon_taxes, payback_times, xlabel = "Carbon Tax (\$/ton)", ylabel = "Payback Time (years)", title = "Carbon Tax vs Solar Payback")

solar_lifespan = 25u"yr"
plot!(carbon_taxes, ones(length(carbon_taxes)) * solar_lifespan, label = "25 year payback")

#hooray, anything above 0.7 dollars per ton of CO2 makes solar more worth it!
#this is assuming that the solar panels and batteries last around 25 years
#however, I think the batteries would have to be replaced much sooner than that, and I don't know what percentage of the cost they make up
#yeesh, they only last about 10-15 years and they make up around 30-50% of the capital cost
#still, basically any carbon tax above probably around 5 dollars per ton of CO2 makes solar more worth it
#oh, the 2 billion dollars per accounts for all costs like battery replacement and maintanence
#even if we go to the high end of 3.5 billion dollars per GW, the carbon tax only has to be 37$/ton of CO2 for solar to be worth it
#also, solar panels can usually last more than 25 years, with some lasting 30-40 years
#also, if concentrated solar power can ever reach a level that allows for 24/7 operation, it would likely be much cheaper than solar panels
#or natural gas