import matplotlib.pyplot as plt
import random as rd

import pandas as pd 

def FuelComp(M):
    """
    FuelComp draws a graph comparing the Energy output (mWh), Pollution toll (kg of CO2),
    and Cost ($) of various fuel sources. These sources are: D-T fusion, Uranium-235 Fission,
    Coal, and Crude Oil.\
    ___
    
    Parameter: M = Mass of fuel (kg)
    (an ideal Deuterium Tritium mix in the case of DT fusion)
    ___
    
    Example: print(FuelComp(1))
    
    """
    # Multiply the various values of the base information matrix to account for mass selection
    data=[["Energy(kWh)", M* 10**8, M*24*(10**6), M*8, M*12], 
          ["CO2(kg)", M*0, M*0, M*3.67, M*3.15], 
          ["Cost($)",M*30 *(10**6), M*15 *(10**6), M*0.3, M*0.68] 
         ]

    # Plot the adjusted information matrix and add legend, etc
    df = pd.DataFrame(data, columns=["Factors", "DT Fusion", "U-235 Fission",
                                     "Coal", "Crude Oil"])

    fig = plt

    df.plot(title="Comparison of " + str(M) + "kg of various fuels",
            x="Factors", y=["DT Fusion", "U-235 Fission", "Coal", "Crude Oil"],
            kind="barh",logx="sym")

    return fig
    
def EnergyDistr(nD,nT):
    """
    EnergyDistr produces a line graph of the energy distribution throughout a multi 
    stage fusion process. The energy of each of the reactants (Deuterium and Tritium)
    and products (Hydrogen and neutron) is graphed against minute fractions of a second.
    (the time inaccuracy is simply non resolvable given our current level)
    ___
    
    Parameters: nD, nT = number of Deuterium and Tritium particles respectively.
    ___
    
    Example: print(EnergyDistr(75,75)) 
    
    """
    # Setting up parameters and converting to energy values.
    nH = 0

    D = nD * 2.2
    T = nT * 2.75 
    H = 0
    E = 0

    # Create arrays to store energy values
    YnD = []
    YnT = []
    YnH = []
    YE = []

    # Begin iterative process of until there are no more of at least one reactant
    while nD * nT > 0:
        
        # Use random numbers to simulate cross section reactivity 
        rnd = rd.random()
        
        if rnd < 0.1:
            
            nD -= 1
            nT -= 1
            
            D = D - 2.2
            T = T - 2.75 
            H += 3.5
            nH += 1
            E = 14.1 * nH
            
            # Remove #s for text stream describing the reaction
            #print("nT = " + str(nT) +", nD = " + str(nD) 
                  #+ ", nH = " + str(nH) + ", E = " + str(E) + "MeV")
            
            # Append arrays to hold energy values 
            YnD.append(D)
            YnT.append(T)
            YnH.append(H)
            YE.append(E)
            
        else:
            
            #print("No reaction")
            
            # Append arrays to hold energy values 
            YnD.append(D)
            YnT.append(T)
            YnH.append(H)
            YE.append(E)
       
    # Plot with appropriate additions
    fig = plt
    plt.plot(YnD, label = "Deuterium Energy")
    plt.plot(YnT, label = "Tritium Energy")
    plt.plot(YnH, label = "Hydrogen Energy")
    plt.plot(YE, label = "Usable Neutron Energy")

    plt.legend()
    plt.title("Energy distribution simulation of D-T fusion")
    plt.ylabel("Energy (MeV)")
    plt.xlabel("Time (fractions of a second)")
    
    return fig
