### RUN SURVIVAL PROB IN A LOOP AND PLOT ##

import numpy as np
import matplotlib.pyplot as plt
import os



####### TEST ######
# Global variables
GF = 1.663788e-14 #(eV^-2)
hbar = 6.58211957e-16 #(eV.s)
c = 299792458 #(m/s)

s12 = np.sqrt(0.307)
s13 = np.sqrt(2.18e-2)
s23 = np.sqrt(0.545)
c12 = np.sqrt(1.0 - s12 * s12)
c13 = np.sqrt(1.0 - s13 * s13)
c23 = np.sqrt(1.0 - s23 * s23)

m21 = 7.53e-5 #(Delta m_21^2, in eV^2)
m32 = 2.453e-3 #(Delta m_32^2, in eV^2)
m31 = m32 + m21 #(Delta m_31^2, in eV^2)

def Matt_Osc_m21(E, L, Ne):
    '''Oscillations with matter effects, assuming \Delta m_{21}^2 dominance'''

    # convert units to eV
    E *= 1e6 #(MeV to eV)
    L *= 1e3 / (c * hbar) #(km to eV^-1)
    Ne *= c*c*c * hbar*hbar*hbar #(m^-3 to ev^3)

    # Calculate Quantities
    A_CC = -2 * np.sqrt(2) * E * GF * Ne

    s_2thetaM_2 = np.sin(np.arctan(np.tan(2 * np.arcsin(s12)) / (1 - (c13*c13 * A_CC)/(np.cos(2 * np.arcsin(s12)) * m21))))

    m_M21 = np.sqrt((m21 * np.cos(2 * np.arcsin(s12)) - c13*c13 * A_CC)**2 + (m21 * np.sin(2 * np.arcsin(s12)))**2)
    m_M = 0.5 * (m21 + c13*c13 * A_CC + m_M21)

    # Survival Probability
    P = 1 - s_2thetaM_2**2 * np.sin((m_M * L)/(4 * E))**2
    P = c13*c13*c13*c13 * P + s13*s13*s13*s13

    return P

def Matt_Osc_m31(E, L, Ne):
    '''Oscillations with matter effects, assuming \Delta m_{21}^2 dominance'''

    # convert units to eV
    E *= 1e6 #(MeV to eV)
    L *= 1e3 / (c * hbar) #(km to eV^-1)
    Ne *= c*c*c * hbar*hbar*hbar #(m^-3 to ev^3)

    # Calculate Quantities
    A_CC = -2 * np.sqrt(2) * E * GF * Ne

    s_2thetaM_2 = np.sin(np.arctan(np.tan(2 * np.arcsin(s13)) / (1 - A_CC/(np.cos(2 * np.arcsin(s13)) * m31))))

    m_M = np.sqrt((m31 * np.cos(2 * np.arcsin(s13)) - A_CC)**2 + (m31 * np.sin(2 * np.arcsin(s13)))**2)

    # Survival Probability
    P = 1 - s_2thetaM_2**2 * np.sin((m_M * L)/(4 * E))**2

    return P

############################



# Create new file (delete first, if it exists)
file_exists = os.path.exists('results.txt')
if file_exists:
    os.system('rm results.txt')
open('results.txt', 'w')

# Run simulation for range of distances (prints results to file)
N = 500
P_21 = np.zeros(N)
P_31 = np.zeros(N)
L = np.linspace(0, 100, N)
# rho_lithosphere = 2.7 #g/cm^3
rho_lithosphere = 0
Ne = rho_lithosphere / 4.821e-15 #electron density cm^-3
Ne *= 1e6 #(cm^-3 to m^3)
E = 1
for i in range(N):
    #print(i)
    os.system("./Matter_Oscillations {} {} {} 1".format(E, L[i], Ne))  #Input neutrino energy (MeV) and propagation length (km)
    #P_21[i] = Matt_Osc_m21(E, L[i], Ne)
    #P_31[i] = Matt_Osc_m31(E, L[i], Ne)

# Read in results
f = open('results.txt', 'r')
data = np.loadtxt(f,skiprows=0)

P = data[:,0]
Pvac = data[:,1]
Pglobes = data[:,2]
Pvacglobes = data[:,3]

#print(L)
#print(Pvac)
# E_eV = E * 1e6 #(MeV to eV)
# Ne *= c*c*c * hbar*hbar*hbar #(m^-3 to ev^3)
# A_CC = -2 * np.sqrt(2) * E_eV * GF * Ne

# L_21 = int((2*np.pi * E_eV * hbar * c)/(m21 * 1000))
# L_31 = int((2*np.pi * E_eV * hbar * c)/(m31 * 1000))


#plt.plot(L, Pvac, label='Vacuum Case (my calc)')
plt.plot(L, P-Pvac, label='With Matter Effects (my calc)', linestyle='dashed')
#plt.plot(L, Pvacglobes - Pvac, label='Vacuum Case (GLoBES)', linestyle='dashed')
plt.plot(L, Pglobes-Pvacglobes, label='With Matter Effects (GLoBES)', linestyle='dashed')
#plt.plot(L, P_21, label=r'$\Delta m_{21}^2(=7.53\times 10^{-5}$eV${}^2)$ dominance approximation')
#plt.plot(L, P_31, label=r'$\Delta m_{31}^2(=2.53\times 10^{-3}$eV${}^2)$  dominance approximation')
plt.xlabel("Baseline (km)")
plt.ylabel(r"$\nu_\mu \rightarrow \nu_e$ transition probability")
plt.legend(loc='best')
Title = 'Comparing Survival Probability for {}MeV neutrino,\n\
    with Matter Effects (constant density {}g/cm^3) to Vacuum Case'.format(E, rho_lithosphere)
plt.title(Title)
plt.show()

