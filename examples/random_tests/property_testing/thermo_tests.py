import numpy as np
import math 
import scipy
import matplotlib.pyplot as plt

from thermo import ChemicalConstantsPackage, CEOSGas, PRMIX, CEOSLiquid, FlashVL, equilibrium, EquilibriumStream, EquilibriumState, eos_mix, Flash


constants, correlations = ChemicalConstantsPackage.from_IDs(['Benzene', 'Toluene'])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
mole_fractions = [0.5, 0.5]

state = flasher.flash(zs=mole_fractions, T=368.1294444943463, P=100000)
#state = flasher.flash(zs = mole_fractions, P = 100000, VF = 0.5)

print(state.Ks(state, state))
print(state.Ks(state, state.gas))
print(state.Ks(state.gas, state.liquid0))
print(state.Ks(state, state.liquid0))
#state.gas.HV

#print(state.mu_l())
print(state.gas.mu())
print(state.liquid0.mu())
print(state.quality)


#flasher.flash(zs=mole_fractions, T=300, P=100000)

temp_range = range(250, 570, 10)
pressure_range = range(10000, 1000000, 10)

resulting_liquid_mu = {} 
resulting_gas_mu = {} 

for temp in temp_range:
    for pressure in pressure_range:
        state = flasher.flash(zs=mole_fractions, T=temp, P=pressure)
        try:
            resulting_liquid_mu[(temp, pressure)] = state.liquid0.mu()
        except:
            resulting_liquid_mu[(temp, pressure)] = None
            
        try:
            resulting_gas_mu[(temp, pressure)] = state.gas.mu()
        except:
            resulting_gas_mu[(temp, pressure)] = None


def create_thermo_interpolation(flasher, temp_range, pressure_range):
    resulting_property_dict = {}

    for temp in temp_range:
        for pressure in pressure_range:
            state = flasher.flash(zs=mole_fractions, T=temp, P=pressure)
            try:
                resulting_property_dict[(temp, pressure)] = state.liquid0.mu()
            except:
                resulting_property_dict[(temp, pressure)] = None            

    return resulting_property_dict

        