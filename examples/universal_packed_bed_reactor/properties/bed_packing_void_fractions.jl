using Unitful

filled_container_volume = 25u"ml"
mass_of_empty_silicon_carbide_packing = 42.027u"g"
mass_of_water_filled_silicon_carbide_packing = (54.222u"g" + 54.978u"g") / 2

water_density = 998u"kg/m^3"

silicon_carbide_bed_void_fraction = ((mass_of_water_filled_silicon_carbide_packing - mass_of_empty_silicon_carbide_packing) / water_density) / filled_container_volume |> u"m^3/m^3"

silicon_carbide_packing_density = mass_of_empty_silicon_carbide_packing / filled_container_volume |> u"kg/m^3"
silicon_carbide_comparison_to_actual_density = mass_of_empty_silicon_carbide_packing / (filled_container_volume - (silicon_carbide_bed_void_fraction * filled_container_volume)) |> u"kg/m^3"

silicon_carbide_density = 3.21u"g/cm^3"
extra_check_packing_density = silicon_carbide_density * (1 - silicon_carbide_bed_void_fraction) |> u"kg/m^3"

silicon_carbide_volume = pi * (16.0u"mm" / 2)^2 * (3.0u"inch") |> u"ml"

bed_void_fraction_second_check = (filled_container_volume - (mass_of_empty_silicon_carbide_packing / silicon_carbide_density)) / filled_container_volume |> u"m^3/m^3"
#0.476, that's pretty close to 0.496, so it's a pretty good approximation
#I think the first is probably going to be more accurate because it accounts for the packing media being porous 

final_silicon_carbide_bed_void_fraction = 0.504u"m^3/m^3"
#0.504
final_silicon_carbide_packing_density = (1 - final_silicon_carbide_bed_void_fraction) * silicon_carbide_density |> u"kg/m^3"
#1592.16 kg m^-3


#COPPER MESH
mass_of_copper_mesh_used = 27.0u"g" #this is not what I got on my scale, but this is the approximate weight of 9 x 9.5 inch piece of 100 mesh copper mesh, which I trust more
copper_density = 8.96u"g/cm^3"

copper_mesh_section_avaliable_volume = pi * (16.0u"mm" / 2)^2 * (9.1u"inch") |> u"ml"
copper_mesh_occupied_volume = mass_of_copper_mesh_used / copper_density |> u"ml"

copper_mesh_bed_void_fraction = (copper_mesh_section_avaliable_volume - copper_mesh_occupied_volume) / copper_mesh_section_avaliable_volume |> u"m^3/m^3"
#0.9351588365422343

total_reactor_volume = pi * (16.0u"mm" / 2)^2 * (12.1u"inch") |> u"ml"
total_reactor_volume_without_silicon_carbide = total_reactor_volume - final_silicon_carbide_bed_void_fraction * silicon_carbide_volume

reactor_avaliable_volume = 38.22u"ml"

copper_mesh_volume = total_reactor_volume_without_silicon_carbide - reactor_avaliable_volume
other_copper_mesh_bed_void_fraction = (copper_mesh_section_avaliable_volume - copper_mesh_volume) / copper_mesh_section_avaliable_volume |> u"m^3/m^3"
#0.656

#I guess, 0.70 is a decent approximation 

final_copper_mesh_bed_void_fraction = 0.65u"m^3/m^3"
#0.65
final_copper_mesh_packing_density = (1 - final_copper_mesh_bed_void_fraction) * copper_density |> u"kg/m^3" 
#3136.0 kg m^-3
