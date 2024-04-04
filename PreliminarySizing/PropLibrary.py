import numpy as np
import rocketcea
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant



# Hydrogen Peroxide  
peroxide98 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[98,2]) 
peroxide96 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[96,4]) 
peroxide95 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[95,5]) 
peroxide90 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[90,10]) 
peroxide85 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[85,15]) 


# Ethaol water mixtures
ethanol90 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[90,10]) 
ethanol85 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[85,15]) 
ethanol80 = rocketcea.blends.newFuelBlend(fuelL=['C2H5OH', 'H2O'], fuelPcentL=[80,20]) 


# Methanol water mixes
# Ethaol water mixtures
methanol90 = rocketcea.blends.newFuelBlend(fuelL=['CH3OH', 'H2O'], fuelPcentL=[90,10]) 
methanol85 = rocketcea.blends.newFuelBlend(fuelL=['CH3OH', 'H2O'], fuelPcentL=[85,15]) 
methanol80 = rocketcea.blends.newFuelBlend(fuelL=['CH3OH', 'H2O'], fuelPcentL=[80,20]) 

# aniline
card_str = """
fuel C6H7N(L)  C 6.0   H 7.0    N 1.0  wt%=100
h, kj/mol=31.3    t(k)=298.15   rho=1.03 
"""
add_new_fuel( 'aniline', card_str )
aniline = rocketcea.blends.newFuelBlend(fuelL=['aniline'], fuelPcentL=[100]) 


# furfurylalcohol
card_str = """
fuel C5H6O2(L)  C 5.0   H 6.0    O 2.0  wt%=100
h, kj/mol=-276.2    t(k)=298.15   rho=1.13 
"""
add_new_fuel( 'furfurylalcohol', card_str )
furfurylalcohol = rocketcea.blends.newFuelBlend(fuelL=['furfurylalcohol'], fuelPcentL=[100]) 

# EMIM SCN 
card_str = """
fuel C7H11N3S(L)  C 7.0   H 11.0   N 3.0   S 1.0  wt%=100
h, kj/mol=52.8    t(k)=298.15   rho=1.11 
"""
add_new_fuel( 'EMIMSCN', card_str )
EMIMSCN = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN'], fuelPcentL=[100]) 

# TMPDA 
card_str = """
fuel C7H18N2(L)  C 7.0   H 18.0   N 2.0   wt%=100
h, kj/mol=-52.75    t(k)=298.15   rho=0.774
"""
add_new_fuel( 'TMPDA', card_str )
TMPDA = rocketcea.blends.newFuelBlend(fuelL=['TMPDA'], fuelPcentL=[100]) 


